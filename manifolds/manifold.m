classdef (Abstract) manifold < handle & matlab.mixin.Heterogeneous
    % An abstract class representing a Manifold.
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
    %    dot(p,v,w)          : inner product of v,w in TpM
    % _Static
    %    addNoise(X,sigma)   : add Noise on the manifold to the signal/data
    %                          X with standard deviation sigma
    % _normal functions
    %    mean(x)              : Karcher mean of the values x
    %    meadian(x)           : median of the values of x
    %    midPoint(x,y)        : Compute the mid point between x and y.
    %    filter(image,filter) : Filters the image using the Karcher mean
    %    geopoint(x,y,t)      : Compute the geodesic between x and y and
    %                           evaluate it at t
    %    geodesic(x,y)        : functional of the geodesic,
    %                           'pts' for a number of equidistant points
    %                           't' for a vector of sampling points
    %    schildsladder(x,y,xi): Transports the vector xi from x to y with
    %                           Schild's ladder
    %    poleladder(x,y,xi)   : Transports the vector xi from x to y with
    %                           the pole ladder
    %
    % ---
    % Manifold-valued Image Restoration Toolbox 1.0
    % R. Bergmann ~ 2014-10-18 | 2018-02-16
    % see LICENSE.txt

    % Changelog
    % 2018-02-16 ? unifies notation
    properties
        useMex = true; %Whether or not to use mex-files in the manifold-functions
    end
    properties (Abstract)
        type;      % (string) A string representing the manifold
        ItemSize;  % (array) dimension of one manifold point
        Dimension; % (int) dimension of the manifold
        allDims;   % (cell) as many colon operators ({':'}) as the length
                   % ItemSize
    end
    methods (Abstract)
        % exponential map of xi in TxM  at x on the manifold
        q = exp(this,x,xi);
        % inverse exponential map of y at x on the manifold
        v = log(this,x,y);
        % distance between p and q on the manifold
        d = dist(this,p,q);
        % inner product of v,w in TpM
        d = dot(this,p,v,w);
        % addNoise(X,sigma) add noise w.r.t. the manifold itself and a
        % standard deviation sigma to a (multidimensional) signal X of
        % manifold valued pixels.
        Y = addNoise(this,X,sigma,varargin)
    end
    methods
        function y = mean(this,varargin)
            % mean(x) calculates the m means of the input data x [.,m,n]
            % with gradient descent algorithm. The implementation follows
            %
            % B. Afsari, Riemannian Lp center of mass: Existence,
            %    uniqueness, and convexity,
            %    Proc. AMS 139(2), pp.655-673, 2011.
            %
            % INPUT
            %    x :  m x n Data points ([this.Itemsize,m,n]) to compute
            %         m means of n points each, pp.
            % OUTPUT
            %    y :  m data points of the means calculated
            %
            % OPTIONAL
            % 'Weights' : (1/n*ones([m,n]) 1xn or mxn weights for the mean
            %            the first case uses the same weights for all means
            % 'InitVal' : m Initial Data points for the gradient descent
            % 'MaxIterations': maximal number of iterations
            % 'Epsilon'      : change of bewteen to iterates to stop
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.0,
            % J. Persch 2015-07-24 | R. Bergmann 2018-02-16

            % Changelog
            % 2018-02-16 Unifies variables, switches to new gradDesc scheme
            % Changelog
            % 2018-02-16 RB; adapted to the new functional gradient descent
            %            scheme ? clearer code.
            ip = inputParser;
            addRequired(ip,'x');
            addParameter(ip,'Weights',NaN);
            addParameter(ip,'InitVal',NaN);
            addParameter(ip,'MaxIterations',50);
            addParameter(ip,'Epsilon',10^(-5));
            parse(ip, varargin{:});
            vars = ip.Results;
            stoppingCriterion = stopCritMaxIterEpsilonCreator(this,...
                vars.MaxIterations,vars.Epsilon);
            x = vars.x;
            xS = size(x);
            mD = length(this.ItemSize);
            m = size(x,mD+1);
            n = size(x,mD+2);
            if xS(1:mD) ~= this.ItemSize
                error(['input data item size (',num2str(xS(1:mD)),...
                    ') does not fit item size of manifold (',...
                    num2str(this.ItemSize),'.']);
            end
            if n==1 % only one point each - direct return
                y=x;
                return
            end
            %
            if isnan(vars.Weights)
                w = 1/n*ones(m,n);
            elseif isvector(vars.Weights) % repmat to all n means
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
                % normalize
                w = w./repmat(sum(w,2),1,n);
            end
            % functional - dist collapses manifold dims, hence we sum over
            % the second dimension
            % F = @(yk) sum(w.*this.dist(repmat(yk,[ones(1,mD+1),n),x).^2,2);
            gradF = @(yk) sum(shiftdim(w,-mD).*gradDistSquared(this,...
                repmat(yk,[ones(1,mD+1),n]),x),mD+2);
            if isnan(vars.InitVal) %no initial value given
                yI = x(this.allDims{:},:,1); %take all first points
            else
                yI = vars.InitVal;
                if size(yI,mD+2)==1 %only one init point given -> repmat
                    yI = repmat(yI,[ones(1,mD),m]);
                end
                yS = size(yI);
                if any(yS(1:mD)~=this.ItemSize) || size(yI,mD+1) ~= m
                   error(['InitVal has to be of size [',num2str(this.ItemSize),'] or [',...
                        num2str([this.ItemSize,m]),'] but is [',num2str(size(yI)),'].']);
                end
            end
            stepSizeRule = @(x,eta,iter,initial) 1/2;
            y = gradientDescent(this,yI,gradF,stepSizeRule,stoppingCriterion);
        end
        function y = median(this,varargin)
            % median(x) calculates the m medians of x ([.,m,n])
            % of the input data with a gradient descent algorithm.
            % This implementation is based on
            %
            % B. Afsari, Riemannian Lp center of mass: Existence,
            %    uniqueness, and convexity,
            %    Proc. AMS 139(2), pp.655-673, 2011.
            % and adapted to the median defined in
            % P. T. Fletcher, S. Venkatasubramanian, and S. Joshi:
            %    The geometric median on Riemannian manifolds with
            %    application to robust atlas estimation.
            %    NeuroImage. 45 S143?S152
            %
            % INPUT
            %    x :  m x n Data points ([this.Itemsize,m,n]) to compute
            %         m means of n points each, pp.
            % OUTPUT
            %    y :  m data points of the medians calculated
            %
            % OPTIONAL
            % 'Weights' : (1/n*ones([m,n]) 1xn or mxn weights for the mean
            %            the first case uses the same weights for all means
            % 'InitVal' : m Initial Data points for the gradient descent
            % 'MaxIterations': (50) Maximal Number of Iterations
            % 'Epsilon'      : (10^(-5)) Maximal change before stopping
            % ---
            % Manifold-valued Image Restoration Toolbox 1.0
            % J. Persch, R. Bergmann 2015-07-24 | 2018-02-17
            ip = inputParser;
            addRequired(ip,'x');
            addParameter(ip,'Weights',NaN);
            addParameter(ip,'InitVal',NaN);
            addParameter(ip,'MaxIterations',100);
            addParameter(ip,'Epsilon',10^-5);
            addParameter(ip,'StepSize',@(x,descentDir,iter,s) 1/iter);
            parse(ip, varargin{:});
            vars = ip.Results;
            f = vars.f;
            dims = size(f);
            mD = length(this.ItemSize);
            dD = dims((mD+1):end);
            if length(dD) ~= 2
                if size(f,mD+2)==1 %only one median to compute
                    x = f;
                    return
                end
                error('f wrong size');
            end
            % shift manDim in first dimension
            m = dD(1);
            n = dD(2);
            if isempty(vars.Weights)
                w = 1/n*ones(m,n);
            elseif isvector(vars.Weights) % repmat to all n means
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
                %normalize
                w = w./repmat(sum(w,2),1,n);
            end
            if isnan(vars.InitVal)
                x = f(this.allDims{:},:,1);
            else
                x = vars.InitVal;
                if size(x,mD+1) == 1 || size(x,mD+2) == 1
                    x = repmat(x,[ones(mD,1),m,1]);
                elseif size(x,mD+1) ~= m || size(x,mD+2) ~= 1
                    error(['too many initial points (expected ', ...
                        num2str([M.ItemSize,m,1]),' but x is of size ',num2str(size(x)),'.']);
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
                iter = 1000;
            end
            stoppingCriterion = stopCritMaxIterEpsilonCreator(this,iter,epsilon);
            F = @(p) sum( w.*this.dist(repmat(p,[ones(1,mD+1),n]),f),2);
            gradF = @(p) ...
                sum(...
                    -shiftdim(w./(... %the following line avoids division by zero
                        this.dist(repmat(p,[ones(1,mD+1),n]),f) ...
                        + double(this.dist(repmat(p,[ones(1,mD+1),n]),f)==0)),-mD)...
                          .*this.log(repmat(p,[ones(1,mD+1),n]),f),...
                    mD+2);
            y = subGradientDescent(this,x,F,gradF, ...
                vars.StepSize,stoppingCriterion);
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
            % Manifold-valued Image Restoration Toolbox 1.0
            % R. Bergmann ~ 2014-10-19 | 2015-01-29
            m = this.exp(x, this.log(x,z)./2);
        end
        function filteredImage = filter(this,varargin)
            % filteredImage = filter(this,image,filter)
            %    Convoles the this-valued image with the filter using
            %    Karcher mean
            %  INPUT
            %    image   manifold Valued image
            %    filter  filter matrix of size (2*n+1)x(2*n+1)
            %
            %  OPTIONAL
            %    'BoundaryCondition' ['nearest'] specify boundary
            %                        conditions:
            %                        (nearest,symmetric,periodic)
            %
            %  OUTPUT
            %    filteredImage   filtered image with Karcher mean
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.2
            % J. Persch, R. Bergmann ~ 2017-04-06 | 2018-02-17

            % Changelog
            % 2018-02-17 RB; rephrase to avoid reshapes and use allDims.

            ip = inputParser;
            addRequired(ip,'image');
            addRequired(ip,'filter');
            addParameter(ip,'BoundaryCondition','nearest');
            parse(ip, varargin{:});
            vars = ip.Results;
            image = vars.image;
            filter = vars.filter;
            bc = vars.BoundaryCondition;

            dimen = size(image);
            imgDim = dimen(length(this.ItemSize)+1:end);
            assert(length(imgDim)==2,'Works only for manifold valued images.');
            fSize = size(filter);
            assert(all(mod(fSize,2)) || fSize(1)==fSize(2),'Filter needs to be of dimension (2*n+1) x (2*n+1).');
            n = (fSize(1)-1)/2;
            switch bc
                case 'nearest'
                    image = image(this.allDims{:},[ones(1,n),1:imgDim(1),ones(1,n)*imgDim(1)],[ones(1,n),1:imgDim(2),ones(1,n)*imgDim(2)]);
                case 'symmetric'
                    image = image(this.allDims{:},[n:-1:1,1:imgDim(1),imgDim(1):-1:imgDim(1)-n+1],[n:-1:1,1:imgDim(2),imgDim(2):-1:imgDim(2)-n+1]);
                case 'periodic'
                    image = image(this.allDims{:},mod(-n:imgDim(1)+n-1,imgDim(1))+1,mod(-n:imgDim(2)+n-1,imgDim(2))+1);
            end
            range = 0:2*n;
            filteredImage = zeros([this.ItemSize,imgDim]);
            for i = 1:imgDim(1)
                for j = 1:imgDim(2)
                    filteredImage(this.allDims{:},i,j) = ...
                        this.mean(reshape(image(:,i+range,j+range),[this.ItemSize,1,prod(fSize)]),...
                        'Weights',filter(:));
                end
            end
        end
        function W = geopoint(this,x,y,t)
            % geopoint(x,y,t) - Gives the point \gamma_{x,y}(t)
            % placeholder, has many manifolds admit a faster way to compute
            % the combination of exp and log, e.g., SymPosDef
            %
            % INPUT
            %   x,y : a point or set of points on the manifold
            %   t   : a scalar or set of scalars
            %
            %
            % OUTPUT
            %   w : resulting point(s)
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0
            % J. Persch, R. Bergmann | 2017-03-31 | 2018-02-16

            % Changelog
            % 2018-02-16 Unifies notation
            % 2015-04-10 introduced mexfile
            W = this.exp(x,this.log(x,y),t);
        end
        function geo = geodesic(varargin)
            % geo = geodesic(this,x,y)
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
            %   geo : the geodesic between x,y evaluated at several points.
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
        function nu = schildsladder(this,x,y,xi)
            % schildsladder(this,x,y,xi) approximates parallel Transport
            % by mappting xi from TxM to TyM using Schild's ladder
            %
            % INPUT
            % x,y : two point(set)s on the manifold
            % v   : a tangent vector(Set) on TxM
            %
            % OUTPUT
            %   nu : the resulting vectors in TyM that are approximately
            %        parallelTransport(x,y,xi)
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.2
            % J. Persch, R. Bergmann | 2018-01-04 | 2018-02-19
            nu = this.log(y,this.geopoint(x,this.midPoint(this.exp(x,xi),y),2));
        end
        function nu = poleladder(this,x,y,xi)
            % poleladder(this,x,y,xi) approximates parallel Transport
            % by mappting xi from TxM to TyM using the pole ladder
            %
            % INPUT
            % x,y : two point(set)s on the manifold
            % xi   : a tangent vector(Set) on TxM
            %
            % OUTPUT
            %   nu : the resulting vectors in TyM that are approximately
            %        parallelTransport(x,y,xi)
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.2
            % J. Persch, R. Bergmann | 2018-01-04 | 2018-02-19
            nu = -this.log(y,this.geopoint(this.exp(x,xi),this.midPoint(x,y),2));
        end
        function [v, mean_f] = var(this,f)
            % var(f) computes the empirical variance
            %       1/(numel(f)-1) * sum (f-mean(f))^2
            % of f
            % INPUT
            % f      : manifold valued Set
            % OUTPUT
            % v      : variance of the set (scalar)
            % mean_f : mean value of the set
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.2
            % J. Persch, R. bergmann, 2017-01-06
            mD = length(this.ItemSize);
            f = f(this.allDims{:},:);
            num_el = size(f,mD+1);
            mean_f = this.mean(reshape(f,[this.ItemSize,1,num_el]));
            v = 1/(num_el-1)*sum(this.dist(repmat(mean_f,[ones(1,mD),num_el]),f).^2)/this.Dimension;
        end
        function xi = JacobiField(this,varargin)
            % JacobiField(x,y,t,eta) - evaluate a Jacobi field
            %    along the geodesic geo(x,y) at point t, where the
            %    'weight'-function f(k,t,d) determines the boundary
            %    conditions of the field and hence the meaning of eta,
            %
            % INPUT
            %    x   : a point on the manifold Sn (or a vector of points)
            %    y   : a point on the manifold Sn (or a vector of points)
            %    t   : a value from [0,1] indicating the point geo(x,y,t)
            %           (or a vector of values)
            %    eta : an initial condition of the Jacobi field, where the
            %           following weights determine the type of initial
            %           condition.
            %
            % OPTIONAL
            %    'weights' : [@(k,t,d) = (k==0)*t
            %       + (k>0)*sin(sqrt(k)*t*d)/sin(sqrt(k)*d)
            %       + (k<0)*sinh(sqrt(-k)*d*t)/sinh(sqrt(-k)*d)
            %       provides the weight depending on the eigenvalue (k) of
            %       the curvature tensor coresponding to the ONB basis
            %       vector, the position t along the Jacobi field and d the
            %       length of the geodesic.
            %       For the standard value, eta is a tangential vector at 0
            %       and the second boundary condition is J(1)=0, i.e. the
            %       Jacobifield corresponds to D_x\gamma_{xy}(t)[\eta]
            %
            % ---
            % Manifold-valued Image Restoration Toolbox 1.3
            % R. Bergmann | MVIRT | 2017-12-01
            ip = inputParser;
            addRequired(ip,'x');
            addRequired(ip,'y');
            addRequired(ip,'t');
            addRequired(ip,'eta');
            % for the weights we have to include two numerical tricks:
            %   for both k=0 or d=0 the second and third terms yield nan
            %   but since the nominator is also zero we have to avoid to
            %   divide by zero.
            %   Furthermore for d==0 (hence x=y) the result is just the
            %   identity, hence the last term, that for d==0 we have just 1
            addParameter(ip,'weights', @(k,t,d) (k==0).*ones(size(k.*t.*d)).*(1-t) + ...
                (k>0).*sin(sqrt(k).*(1-t).*d)./(sin(sqrt(k).*d) + (k==0)  + (d==0)) + ... %last term avoids division by zero
                (k<0).*sinh(sqrt(-k).*d.*(1-t))./(sinh(sqrt(-k).*d) + (k==0) + (d==0)) ...
            );
            parse(ip, varargin{:});
            vars = ip.Results;
            dDim = (length(this.ItemSize)+1); %dimension where the data lives
            l = size(vars.x,dDim); % number of vectors
            if (size(vars.y,dDim) ~= l)
                error(['The lengths of p (',num2str(l),') and q (',...
                    num2str(size(vars.y,dDim)),') have to be the same']);
            else
                t=vars.t;
                if sum(size(t))==2 % number ? extend to l
                    t = vars.t.*ones(l,1);
                elseif (length(t) ~= l)
                error(['The lengths of p (',num2str(l),') and t (',...
                    num2str(length(vars.t)),') have to be the same']);
                end
            end
            myEps = 10^(-8);
            d = this.dist(vars.x,vars.y);
            xi = zeros(size(vars.eta));
            xi(this.allDims{:},d<myEps) = vars.eta(this.allDims{:},d<myEps);
            if sum(d(:)>myEps)>0
                % only continue with large enough ones
                x = vars.x(this.allDims{:},d>=myEps);
                y = vars.y(this.allDims{:},d>=myEps);
                eta = vars.eta(this.allDims{:},d>=myEps);
                t = t(d>myEps);
                dS = d(d>myEps);
                p = this.geopoint(x,y,t);
                [V,k] = this.TpMONB(x,y);
                W = this.parallelTransport(...%repmat both x and p to num TpONBs
                    repmat(x,[ones(1,dDim),this.Dimension]),...
                    repmat(p,[ones(1,dDim),this.Dimension]),...
                    V);
                weights = vars.weights(k,t,dS);
                % decompose eta into basis of v
                alpha = shiftdim(this.dot(repmat(x,[ones(1,dDim),this.Dimension]),...
                    repmat(eta,[ones(1,dDim),this.Dimension]),V),-length(this.ItemSize));
                % shift weights by mandim dims back such that .* work correctly
                % - sum over all tangential vectors i.e. dataDim+1
                xi(this.allDims{:},d>myEps) = sum(alpha.*W.*permute(weights,[3:(dDim+1),1,2]),(dDim+1));
            end
        end
        function xi = DxGeo(this,x,y,t,eta)
            % DxGeo(x,y,t,eta) - Compute the Derivative D_xGeo(t; x,y)[eta]
            %    i.e. of geo(x,y,t) with respect to the start point x.
            %
            %    For a function f: M \mapsto R and fixed y,t we have for the
            %    gradient of g(x) = f(geo(x,y,t)) that
            %    <grad g, nu>_x = <grad f, DxGeo(.,y,t)(x)[nu]>_g(x,y,t)
            %    hence with the Adjoint we obtain
            %    grad g = AdjDxGeo(.,y,t)(x)[grad f].
            %    This function hence only requires eta=grad f to computed
            %    the chain rule.
            %
            %    INPUT
            %      x   : start point of a geodesic, g(x,y,0)=x
            %      y   : end point of a geodesic, geo(x,y,1) = y
            %      t   : [0,1] a point on the geodesic to be evaluated,
            %            may exceed [0,1] to leave the segment between x and y
            %     eta  : (in Tg(t;x,y)M) direction to take the Adjoint derivative at.
            %
            %    OUTPUT
            %     xi   : ( in TxM ) - the adjoint of DxGeo with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04
            xi = this.JacobiField(x,y,t,eta);
        end
        function xi = DyGeo(this,x,y,t,eta)
            % DxGeo(x,y,t,eta) Derivative of the geodesic(x,y,t) wrt y.
            %
            %
            %    INPUT
            %      x   : start point of a geodesic, g(x,y,0)=x
            %      y   : end point of a geodesic, geo(x,y,1) = y
            %      t   : [0,1] a point on the geodesic to be evaluated,
            %            may exceed [0,1] to leave the segment between x and y
            %     eta  : (in TyM) direction to take the derivative of.
            %
            %    OUTPUT
            %     xi   : ( in Tg(x,y,t)M ) - DyGeo with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04
            xi = this.JacobiField(y,x,1-t,eta);
        end
        function xi = AdjJacobiField(this,varargin)
            % AdjJacobiField(x,y,t,eta) - evaluate a adjoint of a Jacobi field
            %    along the geodesic geo(x,y) at point t, where the
            %    'weight'-function f(k,t,d) determines the boundary
            %    conditions of the field and hence the meaning of eta,
            %
            % INPUT
            %    x   : a point on the manifold Sn (or a vector of points)
            %    y   : a point on the manifold Sn (or a vector of points)
            %    t   : a value from [0,1] indicating the point geo(x,y,t)
            %           (or a vector of values)
            %    eta : an initial condition of the Jacobi field, where the
            %           following weights determine the type of initial
            %           condition.
            %
            % OPTIONAL
            %    'weights' : [@(k,t,d) = (k==0)*t
            %       + (k>0)*sin(sqrt(k)*t*d)/sin(sqrt(k)*d)
            %       + (k<0)*sinh(sqrt(-k)*d*t)/sinh(sqrt(-k)*d)
            %       provides the weight depending on the eigenvalue (k) of
            %       the curvature tensor coresponding to the ONB basis
            %       vector, the position t along the Jacobi field and d the
            %       length of the geodesic.
            %       For the standard value, eta is the Jacobi field at 0
            %       and the second boundary condition is J(1)=0, i.e. the
            %       Jacobifield corresponds to D_x\gamma_{xy}(t)[\eta]
            %
            % ---
            % R. Bergmann | MVIRT | 2017-12-01
            ip = inputParser;
            addRequired(ip,'x');
            addRequired(ip,'y');
            addRequired(ip,'t');
            addRequired(ip,'eta');
            % for the weights we have to include two numerical tricks:
            %   for both k=0 or d=0 the second and third terms yield nan
            %   but since the nominator is also zero we have to avoid to
            %   divide by zero.
            %   Furthermore for d==0 (hence x=y) the result is just the
            %   identity, hence the last term, that for d==0 we have just 1
            addParameter(ip,'weights', @(k,t,d) (k==0).*ones(size(k.*t.*d)).*(1-t) + ...
                (k>0).*sin(sqrt(k).*(1-t).*d)./(sin(sqrt(k).*d) + (k==0)  + (d==0)) + ... %last term avoids division by zero
                (k<0).*sinh(sqrt(-k).*d.*(1-t))./(sinh(sqrt(-k).*d) + (k==0) + (d==0)) ...
            );
            parse(ip, varargin{:});
            vars = ip.Results;
            dDim = (length(this.ItemSize)+1); %dimension where the data lives
            xS = size(vars.x);
            %internally reshape to a vector
            x = reshape(vars.x,[this.ItemSize, prod(xS(dDim:end))]);
            y = reshape(vars.y,[this.ItemSize, prod(xS(dDim:end))]);
            eta = reshape(vars.eta,[this.ItemSize, prod(xS(dDim:end))]);
            l = size(x,dDim); % number of vectors
            if (size(y,dDim) ~= l)
                error(['The lengths of p (',num2str(l),') and q (',...
                    num2str(size(y,dDim)),') have to be the same']);
            else
                t=vars.t;
                if sum(size(t))==2 % number ? extend to l
                    t = vars.t.*ones(l,1);
                elseif (length(t) ~= l)
                error(['The lengths of p (',num2str(l),') and t (',...
                    num2str(length(vars.t)),') have to be the same']);
                end
            end
            myEps = 10^(-8);
            d = this.dist(x,y);
            xi = zeros(size(eta));
            if sum(d(:)>myEps)>0
               % only PT large enough ones
               p = this.geopoint(x,y,t);
               [V,k] = this.TpMONB(x,y);
               W=V;
               W(this.allDims{:},d>=myEps,:) = ...
                  this.parallelTransport(...%repmat both x and p to num TpONBs
                    repmat(x(this.allDims{:},d>=myEps),[ones(1,dDim),this.Dimension]),...
                    repmat(p(this.allDims{:},d>=myEps),[ones(1,dDim),this.Dimension]),...
                    V(this.allDims{:},d>=myEps,:));
                weights = vars.weights(k,t,d);
                % decompose eta into basis of v
                alpha = shiftdim(this.dot(repmat(p,[ones(1,dDim),this.Dimension]),...
                    repmat(eta,[ones(1,dDim),this.Dimension]),W),-length(this.ItemSize));
                % shift weights by mandim dims back such that .* work correctly
                % - sum over all tangential vectors i.e. dataDim+1
                xi = sum(alpha.*V.*permute(weights,[3:(dDim+1),1,2]),(dDim+1));
                xi = reshape(xi,xS);
            end
        end
        function xi = AdjDxGeo(this,x,y,t,eta)
            % AdjDxGeo(x,y,t,eta) Adjoint of the Derivative of geo(x,y,t) wrt x.
            %
            %    For a function f: M \mapsto R and fixed y,t we have for the
            %    gradient of g(x) = f(geo(x,y,t)) that
            %    <grad g, nu>_x = <grad f, DxGeo(.,y,t)(x)[nu]>_g(x,y,t)
            %    hence with the Adjoint we obtain
            %    grad g = AdjDxGeo(.,y,t)(x)[grad f].
            %    This function hence only requires eta=grad f to computed
            %    the chain rule.
            %
            %    INPUT
            %      x   : start point of a geodesic, g(x,y,0)=x
            %      y   : end point of a geodesic, geo(x,y,1) = y
            %      t   : [0,1] a point on the geodesic to be evaluated,
            %            may exceed [0,1] to leave the segment between x and y
            %     eta  : (in Tg(t,x,y)) direction to take the Adjoint derivative at.
            %
            %    OUTPUT
            %     xi   : ( in TxM ) - the adjoint of DxGeo with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04
            xi = this.AdjJacobiField(x,y,t,eta);
        end
        function xi = AdjDyGeo(this,x,y,t,eta)
            % AdjDyGeo(x,y,t,eta) - Adjoint of the Derivative of geo(x,y,t) wrt y.
            %
            %    For a function f: M \mapsto R and fixed x,t we have for the
            %    gradient of g(y) = f(geo(x,y,t)) that
            %    <grad g, nu>_y = <grad f, DyGeo(x,.,t)(y)[nu]>_g(x,y,t)
            %    hence with the Adjoint we obtain
            %    grad g = AdjDxGeo(x,.,t)(y)[grad f].
            %    This function hence only requires eta=grad f to computed
            %    the chain rule.
            %
            %    INPUT
            %      x   : start point of a geodesic, g(x,y,0)=x
            %      y   : end point of a geodesic, geo(x,y,1) = y
            %      t   : [0,1] a point on the geodesic to be evaluated,
            %            may exceed [0,1] to leave the segment between x and y
            %     eta  : (in Tg(x,y,t)) direction to take the Adjoint derivative at.
            %
            %    OUTPUT
            %     xi   : ( in TyM ) - the adjoint of DyGeo with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04
            xi = this.AdjJacobiField(y,x,1-t,eta);
        end
        function nu = AdjDxExp(this,x,xi,eta)
            %   nu = AdjDxExp(x,xi,eta) - Adjoint of the Derivative of Exp with
            %   respect to the basis point
            %    INPUT
            %      x   : base point of the exponential
            %      xi  : direction of the exponential
            %     eta  : (in TExp(x,xi)M) direction to take the Adjoint derivative at.
            %
            %    OUTPUT
            %     nu   : ( in TxM ) - the adjoint of DxExp with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04
            f = @(k,t,d) (k==0).*ones(size(k.*t.*d)) + ...
                (k>0).*cos(sqrt(k).*d.*t) + ...
                (k<0).*cosh(sqrt(-k).*d.*t);
            nu = this.AdjJacobiField(x,this.exp(x,xi),1,eta,'weights',f);
        end
        function nu = AdjDxiExp(this,x,xi,eta)
            %   nu = AdjDxExp(x,xi,eta) - Adjoint of the Derivative of Exp with
            %   respect to the tangential vector xi
            %   INPUT
            %      x   : base point of the exponential
            %      xi  : direction of the exponential
            %     eta  : (in TExp(x,xi)M) direction to take the Adjoint derivative at.
            %
            %    OUTPUT
            %     nu   : ( in TxM, more precisely TTxM ) - the adjoint of DxExp with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04
            f = @(k,t,d) (k==0).*ones(size(k.*t.*d)) + ...
                (k>0).*sin(sqrt(k).*d)./(sqrt(k).*d + (d==0) + (k==0) ) + ... % last terms are again for avoiding division by zero
                (k<0).*sinh(sqrt(-k).*d)./(sqrt(-k).*d + (d==0) + (k==0) );
            nu = this.AdjJacobiField(x,this.exp(x,xi),1,eta,'weights',f);
        end
        function xi = AdjDxLog(this,x,y,eta)
            %   nu = AdjDxLog(x,y,eta) - Adjoint of the Derivative of Log
            %       with respect to the basis point x.
            %   INPUT
            %      x   : base point of the logarithm
            %      y   : argument of the logarithm
            %     eta  : (in TxM) direction to take the Adjoint derivative at.
            %
            %    OUTPUT
            %     nu   : ( in TxM ) - the adjoint of DxLog with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04

            % the following weights of the Jacobi field are derived in
            % [Bredies, Holler, Storath, Weinmann, 2017, Lemma 4.5]
            f = @(k,t,d) -(k==0).*ones(size(k.*t.*d)) + ...
                (k>0).*sqrt(k).*(-d).*cos(sqrt(k).*d)./(sin(sqrt(k).*d) + (d==0) + (k==0) ) + ... % last terms are again for avoiding division by zero
                (k<0).*sqrt(-k).*(-d).*cosh(sqrt(-k).*d)./(sinh(sqrt(-k).*d) + (d==0) + (k==0) );
            xi = this.AdjJacobiField(x,y,0,eta,'weights',f);
        end
        function xi = AdjDyLog(this,x,y,eta)
            %   nu = AdjDyLog(x,y,eta) - Adjoint of the Derivative of Log
            %       with respect to y.
            %   INPUT
            %      x   : base point of the logarithm
            %      y   : argument of the logarithm
            %     eta  : (in TyM) direction to take the Adjoint derivative at.
            %
            %    OUTPUT
            %     nu   : ( in TxM ) - the adjoint of DxLog with respect to eta
            % ---
            % MVIRT R. Bergmann, 2017-12-04

            % the following weights of the Jacobi field are derived in
            % [Bredies, Holler, Storath, Weinmann, 2017, Lemma 4.3]
            f = @(k,t,d) (k==0).*ones(size(k.*t.*d)) + ...
                (k>0).*sqrt(k).*d./(sin(sqrt(k).*d) + (d==0) + (k==0) ) + ... % last terms are again for avoiding division by zero
                (k<0).*sqrt(-k).*d./(sinh(sqrt(-k).*d) + (d==0) + (k==0) );
            xi = this.JacobiField(y,x,1,eta,'weights',f);
        end
    end
end
