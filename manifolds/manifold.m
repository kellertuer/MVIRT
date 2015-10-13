classdef (Abstract) manifold < handle
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
    %    mean(f,q)            : Karcher mean of the values f with inital q
    %    midPoint(x,z)        : Compute the mid point between x and z.
    %
    % ---
    % ManImRes 1.0 ~ R. Bergmann ~ 2014-10-18 | 2015-01-25
    
    properties
        useMex = true; %Whether or not to use mex-files in the manifold-functions
    end
    properties (Abstract)
        type; % Type of manifold
        ItemSize; % Data Item dimension, given as size(q), q from M; might differ from manifold dimensions
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
        % add noise w.r.t. manifold and standard deviation sigma to a
        % (multidimensional) signal X
        Y = addNoise(X,sigma)
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
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19
            if all(g(:) == f(:))
                x=f;
                return
            end
            t = this.dist(g,f);
            if (size(t,1)==1) %column vector
                t = shiftdim(t,1);
            end
            tS = size(t); tS = tS(tS>1);
            v = this.log(g,f);
            % repeat into all manifold dimensions
            %first part produces singleton dimensions - as many as the
            %manifold needs, last reproduces t in permute
            t = repmat(... %permute singleton dimensions (as many as manifold dims) from last to first
                    permute(t,[(length(size(f))-length(this.ItemSize)+1):(length(size(f))),1:(length(tS))]), ...
                    [this.ItemSize(this.ItemSize>1),ones(1,length(tS))]);
            x = this.exp(g, (lambda./(1+lambda)).*t.*v);
        end
        function x = mean(this,varargin)
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
            % ManImRes 1.0, J. Persch 2015-07-24 | R. Bergmann 2015-07-30
            ip = inputParser;
            addRequired(ip,'f');
            addParameter(ip,'Weights',NaN);
            addParameter(ip,'InitVal',NaN);
            addParameter(ip,'MaxIterations',50);
            addParameter(ip,'Epsilon',10^-6);
            parse(ip, varargin{:});
            vars = ip.Results;
            f = vars.f;
            dims = size(f);
            if length(dims) ~= length(this.ItemSize)+2
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
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19 | 2015-01-29
            m = this.exp(x, this.log(x,z)./2);
        end
    end
end