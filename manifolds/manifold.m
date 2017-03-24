classdef (Abstract) manifold < handle & matlab.mixin.Heterogeneous
    % A Manifold.
    % This class provides all generell functions available on manifolds and
    % implements those, that are only based on these function interfaces
    %
    % Instances of manifolds have to provide these functions.
    %
    % PROPERTIES
    %   type : a string indicating the type of manifold
    %
    % FUNCTIONS
    % _Abstract, Static
    %    exp(p,v)            : exponential map at p w.r.t v.
    %    log(p,q)            : inverse exponential map of q at p
    %    dist(p,q)           : distance of p and q on the manifold
    % _Static
    %    proxad(f,lambda,w)  : proximal mappinf w.r.t. absolute differences
    %    addNoise(X,sigma)   : add Noise on the manifold to the signal/data
    %                          X with standard deviation sigma
    % _normal functions
    %    proxDist(g,f,lambda) : proximal mapping of distance terms
    %    mean(f)              : Karcher mean of the values f
    %    meadian(f)           : median of the values of f
    %    midPoint(x,z)        : Compute the mid point between x and z.
    %    geodesic(x,y,pts)    : Compute the geodesic between x and y with
    %                           pts points
    %
    % ---
    % Manifold-valued Image Restoration Toolbox 1.0
    % R. Bergmann ~ 2014-10-18 | 2016-10-07
    % see LICENSE.txt
    properties
        useMex = true; %Whether or not to use mex-files in the manifold-functions
    end
    properties (Abstract)
        type; % Type of manifold
        ItemSize; % Data Item dimension, given as size(q), q from M; might differ from manifold dimensions
        Dimension; % Dimension of the manifold
    end
    methods (Abstract)
        % proxad(f,lambda,w,<options>)
        %     perform a proximal mapping with respect to absolute differences
        %     given by w of data f and parameter lambda on the manifold. Any
        %     manifold should also provide the function as
        % proxad(problem), where problem is a struct containing all necessary
        %     parameters.
        x = proxad(this,f,lambda,w,varargin);
        % exponential map of v in TpM  at p on the manifold
        q = exp(this,p,v);
        % inverse exponential map of q at p on the manifold
        v = log(this,p,q);
        % distance between p and q on the manifold
        d = dist(this,p,q);
        % addNoise(X,sigma) add noise w.r.t. the manifold itself and a
        % standard deviation sigma to a (multidimensional) signal X of
        % manifold valued pixels.
        Y = addNoise(this,X,sigma,varargin)
    end
    methods
        function x = proxDist(this,g,f,lambda)
            % proxDist(g,f,lambda)
            % Proximal step towards f from given data g with parameter
            % lambda on an arbitrary manifold. This is the proximal map of
            % the distance function squared with fixed g.
            % INPUT
            %  f,g    : data point( sets/columns )
            %  lambda : stepsize towards f
            % OUTPUT
            %       x : result point( sets) of the proximal map
            % ---
            % Manifold-valued Image Restoration Toolbox 1.0 ~ R. Bergmann, 2014-10-19
            if all(g(:) == f(:))
                x=f;
                return
            end
             v = this.log(g,f);
             if sum(size(lambda))==2 % a numner
                 t = lambda/(1+lambda);
             else
%                t = repmat(lambda./(1+lambda),[ones(1,length(size(lambda))),this.ItemSize]);
                 t = lambda./(1+lambda); %avoid repmat
                 % manifolds first with one-dim stuff to avoid repmat
                 l = length(size(lambda));
                 t = permute(t,[l+(1:length(size(this.ItemSize))),1:l]);
             end
             x = this.exp(g, t.*v);
        end
        function [x1,x2] = proxTVSq(this,f1,f2,lambda)
            % proxTV(f,lambda)
            % Proximal steps of the discrete total variation squared of
            % f1 and f2 with parameter lambda on an arbitrary manifold.
            % This is the proximal map of
            % the piointwise distance function squared d^2(f1,f2).
            % INPUT
            %  f1,f2    : data columns
            %  lambda   : proxParameter
            % OUTPUT
            %  x1,x2    : resulting columns of the proximal map
                        % ---
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19
            if all(f1(:) == f2(:))
                x1=f1;
                x2=f2;
                return
            end
            % Calculate step length in (0,1/2]
            step = lambda/(1+2*lambda)*this.dist(f1,f2);
            step = permute(repmat(step(:),[1,this.ItemSize]),[length(this.ItemSize)+1:-1:1]);
            x1 = this.exp(f1,step.*this.log(f1,f2));
            x2 = this.exp(f2,step.*this.log(f2,f1));
        end
        function [x1,x2] = proxTV(this,f1,f2,lambda)
            % proxTV(f,lambda)
            % Proximal steps for the total variation term of f1 and f2 with
            % parameter lambda on an arbitrary manifold. This is
            % the proximal map of the piointwise distance function d(f1,f2).
            % INPUT
            %  f1,f2    : data columns
            %  lambda   : proxParameter
            % OUTPUT
            %  x1,x2    : resulting columns of the proximal map
            % ---
            % Manifold-valued Image Restoration Toolbox 1.0 ~ R. Bergmann, 2014-10-19
            if all(f1(:) == f2(:))
                x1=f1;
                x2=f2;
                return
            end

            % Calculate step length in (0,1/2]
            step = min(1/2,lambda./this.dist(f1,f2));
            step = permute(repmat(step(:),[1,this.ItemSize]),[length(this.ItemSize)+1:-1:1]);
            x1 = this.exp(f1,step.*this.log(f1,f2));
            x2 = this.exp(f2,step.*this.log(f2,f1));
        end
        function x = mean(this,varargin)
            x = this.mean_gd(varargin{:});
        end
        function x = median(this,varargin)
            % mean(f) calculates the mean of the input data with a gradient
            % descent algorithm. This implementation is based on
            %
            % B. Afsari, Riemannian Lp center of mass: Existence,
            %    uniqueness, and convexity,
            %    Proc. AMS 139(2), pp.655-673, 2011.
            % and adapted to the median defined in
            % P. T. Fletcher, S. Venkatasubramanian, and S. Joshi:
            %    The geometric median on Riemannian manifolds with
            %    application to robust atlas estimation.
            %
            % INPUT
            %    f :  m x n Data points ([this.Itemsize,m,n]) to compute
            %         m means of n points each, pp.
            % OUTPUT
            %    x :  m data points of the medians calculated
            %
            % OPTIONAL
            % 'Weights' : (1/n*ones([m,n]) 1xn or mxn weights for the mean
            %            the first case uses the same weights for all means
            % 'InitVal' : m Initial Data points for the gradient descent
            % 'MaxIterations': Maximal Number of Iterations
            % 'Epsilon'      : Maximal change before stopping
            % 'Alpha'        : Step Size in (0,2)
            %
            % Manifold-valued Image Restoration Toolbox 1.0, J. Persch 2015-07-24 | R. Bergmann 2015-07-30
            ip = inputParser;
            addRequired(ip,'f');
            addParameter(ip,'Weights',NaN);
            addParameter(ip,'InitVal',NaN);
            addParameter(ip,'Alpha',1);
            addParameter(ip,'MaxIterations',100);
            addParameter(ip,'Epsilon',10^-5);
            parse(ip, varargin{:});
            vars = ip.Results;
            f = vars.f;
            dims = size(f);
            if length(dims) ~= length(this.ItemSize)+2
                if all(dims(1:length(this.ItemSize)) == this.ItemSize) && length(dims)<length(this.ItemSize)+2
                    x = f;
                    return
                end
                error('f wrong size');
            end
            % shift manDim in first dimension
            f = reshape(f,[prod(this.ItemSize),dims(1+length(this.ItemSize):end)]);
            con_dim = size(f);
            m = con_dim(end-1);
            n = con_dim(end);
            if isnan(vars.Weights)
                w = 1/n*ones(m,n);
            elseif isvector(vars.Weights)
                w = vars.Weights(:)';
                if length(w) ~=n
                    error('length(w) does not match data points');
                end
                w = repmat(w,m,1);
            else
                w = vars.Weights;
                if any(size(w) ~= [m,n])
                    error('dim w do not match data points');
                end
                w = w./repmat(sum(w,2),1,n);
            end
            % Resize w to fit to the Manifold
            w = repmat(permute(w,[2+(1:length(this.ItemSize)),1,2]),[this.ItemSize,1,1]);
            if isnan(vars.InitVal)
                x = reshape(f,[prod(this.ItemSize),m,n]);
                x = reshape(x(:,:,1),[this.ItemSize,m]);
            else
                x = vars.InitVal;
                if (length(size(x))== length(this.ItemSize)) ...
                        && (all(size(x) == this.ItemSize))
                    x = repmat(x,[ones(1,length(this.ItemSize)),m]);
                elseif any(size(x) ~= [this.ItemSize,m])
                    error(['InitVal has to be of size [',num2str(this.ItemSize),'] or [',...
                        num2str([this.ItemSize,m]),'].']);
                end
            end
            if vars.Epsilon > 0
                epsilon = vars.Epsilon;
            else
                warning('Epsilon should be larger than zero, set Epsilon to 10^-6')
                epsilon = 10^-6;
            end
            if vars.MaxIterations > 0
                iter = vars.MaxIterations;
            else
                warning('Iterations should be larger than zero, set Iterations to 100')
                iter = 100;
            end
            f  = reshape(f, [this.ItemSize,m,n]);
            x_old = x;
            i = 0;
            while (max(this.dist(x,x_old)) > epsilon && i < iter) || i == 0
                x_old = x;
                V = this.log(repmat(x,[ones(1,length(this.ItemSize)+1),n]),f);
                % Divide by the distance
                d = this.dist(repmat(x,[ones(1,length(this.ItemSize)+1),n]),f);
                l = length(this.ItemSize);
                d = repmat(permute(d,[2+(1:l),1,2]),[this.ItemSize,1,1]);
                V(d>0) = V(d>0)./d(d>0);
                V = V.*w;
                weight = sum(d.*w,length(this.ItemSize)+2);
                V = sum(V,length(this.ItemSize)+2);
                x = this.exp(x,vars.Alpha*weight.*V);
                i= i+1;
            end
        end
        function m = midPoint(this,x,z)
            % m = midPoint(x,z)
            %   Compute the (geodesic) mid point of x and z.
            %
            % INPUT
            %    x,z : two point(sets) of manifold points
            %
            % OUTPUT
            %      m : resulting mid point( sets)
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.0 ~ R. Bergmann ~ 2014-10-19 | 2015-01-29
            m = this.exp(x, this.log(x,z)./2);
        end
        function geo = geodesic(varargin)
            % geo = geodesic(this,x,y,pts)
            % Compute the geodesic between x and y using pts-2 points to
            % interpolate
            %
            % INPUT
            %   x,y : two points of the manifold
            % OPTIONAL
            %   pts : (100) optional length of geodesic, or set to length t
            %               if t is chosen
            %     t : vector of points lead to geo = \gamma_{x,y}(t)
            %
            % OUTPUT
            %   geo : the geodesic between x,y with geo(1) = x, geo(pts)=y
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.0 ~ J. Persch ~2015-10-29
            ip = inputParser;
            addRequired(ip,'this');
            addRequired(ip,'x');
            addRequired(ip,'y');
            addParameter(ip,'pts',100);
            addParameter(ip,'t',[]);
            parse(ip, varargin{:});
            vars = ip.Results;

            this = vars.this;
            x = vars.x;
            y = vars.y;

            if size(x,length(this.ItemSize)+1)>1 ||size(y,length(this.ItemSize)+1)>1
                error('x,y should be points on the manifold');
            end
            if vars.pts > 0 || isscalar(vars.pts)
                pts = vars.pts;
            else
               error('pts should be a scalar >0');
            end
            if isvector(vars.t)
                t = vars.t;
                pts = length(t);
            elseif isempty(vars.t)
                t = 0:1/(pts-1):1;
            else
               error('t should be a vector');
            end
           geo = zeros([prod(this.ItemSize),pts]);
           v = this.log(x,y);
           for i = 1:pts
               geo(:,i) = reshape(this.exp(x,v,t(i)),prod(this.ItemSize),1);
           end
           geo = reshape(geo,[this.ItemSize,pts]);
        end
        function [v, mean_f] = var(this,f)
            % var(f) computes the empirical variance
            %       1/(numel(f)-1) * sum (f-mean(f))^2
            % of f
            % INPUT:
            % f      Manifold valued Set
            % OUTPUT:
            % v      variance of the set (scalar)
            % mean_f mean value of the set
            %
            % Manifold-valued Image Restoration Toolbox 1.2 - J. Persch, R. bergmann, 2017-01-06
            dimen = size(f);
            num_el = prod(dimen(length(this.ItemSize)+1:end));
            f = reshape(f,[this.ItemSize,1,num_el]);
            mean_f = this.mean(f);
            v = 1/(num_el-1)*sum(this.dist(repmat(mean_f,[ones(1,length(this.ItemSize)),1,num_el]),f).^2)/this.Dimension;
        end
    end
    methods (Access=protected)
        function x = mean_gd(this,varargin)
            % mean(f) calculates the mean of the input data with a gradient
            % descent algorithm. This implementation is based on
            %
            % B. Afsari, Riemannian Lp center of mass: Existence,
            %    uniqueness, and convexity,
            %    Proc. AMS 139(2), pp.655-673, 2011.
            %
            % INPUT
            %    f :  m x n Data points ([this.Itemsize,m,n]) to compute
            %         m means of n points each, pp.
            % OUTPUT
            %    x :  m data points of the means calculated
            %
            % OPTIONAL
            % 'Weights' : (1/n*ones([m,n]) 1xn or mxn weights for the mean
            %            the first case uses the same weights for all means
            % 'InitVal' : m Initial Data points for the gradient descent
            % 'MaxIterations': Maximal Number of Iterations
            % 'Epsilon'      : Maximal change before stopping
            %
            %
            % Manifold-valued Image Restoration Toolbox 1.0, J. Persch 2015-07-24 | R. Bergmann 2015-07-30
            ip = inputParser;
            addRequired(ip,'f');
            addParameter(ip,'Weights',NaN);
            addParameter(ip,'InitVal',NaN);
            addParameter(ip,'MaxIterations',50);
            addParameter(ip,'Epsilon',5*10^-7);
            parse(ip, varargin{:});
            vars = ip.Results;
            f = vars.f;
            dims = size(f);
            if length(dims) ~= length(this.ItemSize)+2
                if all(dims(1:length(this.ItemSize)) == this.ItemSize) && length(dims)<length(this.ItemSize)+2
                    x = f;
                    return
                end
                error('f wrong size');
            end
            % shift manDim in first dimension
            f = reshape(f,[prod(this.ItemSize),dims(1+length(this.ItemSize):end)]);
            con_dim = size(f);
            m = con_dim(end-1);
            n = con_dim(end);
            if isnan(vars.Weights)
                w = 1/n*ones(m,n);
            elseif isvector(vars.Weights)
                w = vars.Weights(:)';
                if length(w) ~=n
                    error('length(w) does not match data points');
                end
                w = repmat(w,m,1);
            else
                w = vars.Weights;
                if any(size(w) ~= [m,n])
                    error('dim w do not match data points');
                end
                w = w./repmat(sum(w,2),1,n);
            end
            % Resize w to fit to the Manifold
            w = repmat(permute(w,[2+(1:length(this.ItemSize)),1,2]),[this.ItemSize,1,1]);
            if isnan(vars.InitVal)
                x = reshape(f,[prod(this.ItemSize),m,n]);
                x = reshape(x(:,:,1),[this.ItemSize,m]);
            else
                x = vars.InitVal;
                if (length(size(x))== length(this.ItemSize)) ...
                        && (all(size(x) == this.ItemSize))
                    x = repmat(x,[ones(1,length(this.ItemSize)),m]);
                elseif any(size(x) ~= [this.ItemSize,m])
                    error(['InitVal has to be of size [',num2str(this.ItemSize),'] or [',...
                        num2str([this.ItemSize,m]),'].']);
                end
            end
            if vars.Epsilon > 0
                epsilon = vars.Epsilon;
            else
                warning('Epsilon should be larger than zero, set Epsilon to 10^-6')
                epsilon = 10^-6;
            end
            if vars.MaxIterations > 0
                iter = vars.MaxIterations;
            else
                warning('Iterations should be larger than zero, set Iterations to 100')
                iter = 100;
            end
            f  = reshape(f, [this.ItemSize,m,n]);
            x_old = x;
            i = 0;
            while (max(this.dist(x,x_old)) > epsilon && i < iter) || i == 0
                x_old = x;
                V = this.log(repmat(x,[ones(1,length(this.ItemSize)+1),n]),f);
                V = V.*w;
                V = sum(V,length(this.ItemSize)+2);
                x = this.exp(x,V);
                i= i+1;
            end
        end
        function x = mean_cppa(this,varargin)
            % mean(f,lambda)
            % Perform the (karcher) mean on an arbitrary manifold. Computed
            % by employing the cyclic proximal point algorithm from [1]
            %
            % INPUT
            %    f     :  m x n Data points ([this.Itemsize,m,n]) to compute
            %            m means of n points each, pp.
            %   lambda : initial value for the sequence lambda_k in the CPPA
            %
            % OUTPUT
            %   x     : manifold mean of the values in (each column of) f depending on q
            %
            % OPTIONAL PARAMETERS
            %   'InitVal'       : m Initial Data points for the cppa
            %   'Weights'       : [ones(size(f,2))] perform a karcher mean
            %                   with differen weights.
            %   'MaxIterations' : (400) Maximal number of Iterations
            %   'Epsilon'       : (10^{-6}) Lower bound for the max change in one cycle
            %
            %           one of the parameters MaxIterations/Epsilon can be dactivated
            %           by specifying them as Inf but not both.
            %
            %   For Details on the algorithm see
            % [1] M. Bacak, Computing medians and means in Hadamard spaces.
            % SIAM Journal on Optimization, 24(3), pp. 1542-1566, 2013.
            % ---
            % Manifold-valued Image Restoration Toolbox 1.0 ~ R. Bergmann ~ 2015-07-16
            ip = inputParser;
            addRequired(ip,'f');
            addRequired(ip,'lambda');
            addParameter(ip,'Weights',NaN);
            addParameter(ip,'InitVal',NaN);
            addParameter(ip,'MaxIterations',500);
            addParameter(ip,'Epsilon',5*10^-9);
            parse(ip, varargin{:});
            vars = ip.Results;
            f = vars.f;
            dims = size(f);
            if length(dims) ~= length(this.ItemSize)+2
                error('f wrong size');
            end
            f = reshape(f,[prod(this.ItemSize),dims(1+length(this.ItemSize):end)]);
            m = dims(end-1);
            n = dims(end);
            manDim = dims(1:end-2);
            if isnan(vars.Weights)
                w = 1/n*ones(m,n);
            elseif isvector(vars.Weights)
                w = vars.Weights(:)';
                if length(w) ~=n
                    error('length(w) does not match data points');
                end
                w = repmat(w,m,1);
            else
                w = vars.Weights;
                if any(size(w) ~= [m,n])
                    error('dim w do not match data points');
                end
            end
            % Resize w to fit to the Manifold
            if isnan(vars.InitVal)
                x = reshape(f,[prod(this.ItemSize),m,n]);
                x = reshape(x(:,:,1),[this.ItemSize,m]);
            else
                x = vars.InitVal;
                if (length(size(x))== length(this.ItemSize)) ...
                        && (all(size(x) == this.ItemSize))
                    x = repmat(x,[ones(1,length(this.ItemSize)),m]);
                elseif any(size(x) ~= [this.ItemSize,m])
                    error(['InitVal has to be of size [',num2str(this.ItemSize),'] or [',...
                        num2str([this.ItemSize,m]),'].']);
                end
            end
            if vars.Epsilon > 0
                epsilon = vars.Epsilon;
            else
                warning('Epsilon should be larger than zero, set Epsilon to 10^-6')
                epsilon = 10^-6;
            end
            if vars.MaxIterations > 0
                maxiter = vars.MaxIterations;
            else
                warning('Iterations should be larger than zero, set Iterations to 100')
                maxiter = 100;
            end
            itD = inf; %max Distance in last iteration
            i=0;
            while ( (max(itD(:))>=epsilon) && (i<maxiter) )
                i = i+1;
                lambdait = vars.lambda/i;
                xold = x;
                for j=1:n %cycle through all proxes
                    aj = reshape(f,[prod(manDim),m,n]); %collapse manifold
                    %collect all jth points, expand manifold again
                    aj = reshape(aj(:,:,j),[manDim,m]);
                    x = this.proxDist(x,aj,lambdait*w(:,j));
                end
                itD = this.dist(x,xold);
                itD = max(itD(:));
            end
        end
    end
end
