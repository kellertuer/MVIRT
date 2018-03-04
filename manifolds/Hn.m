classdef Hn < manifold & handle
    % The n-dimensional hyperbolic manifold H(n) embedded in R^{n+1}
    %    
    % ---
    % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~2015-10-20
    properties
        type = 'Hn';
        tau = 1;
        steps = 50;
        ItemSize;
        Dimension;
        allDims;
    end
    
    methods
        function obj = Hn(n)
            % Create a Manifold for teh n-dimensional hyperbolic manifold
            obj.type = ['hyperbolic manifold of dimension ',num2str(n)];
            obj.ItemSize = n+1;
            obj.Dimension = n;
            obj.allDims = repelem({':'},length(obj.ItemSize));
        end
        function Y = exp(this,X,V,t)
            % exp(P,V) - Exponential map at the point(s) P with respect to
            %    the direction(s) V in TXHn
            %
            % INPUT
            %   X : a point or set of points on the manifold Hn
            %   V : a point or set of point in the tangential spaces TXHn
            %
            % OPTIONAL
            %   t : shorten vectors V by factor t
            %   (given as one value or an array of same length as number of Vs)
            %
            % OUTPUT
            %   Y : resulting point(s) on Hn
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1, R. Bergmann ~ 2015-01-20 | 2015-10-20
            if nargin<4 % define t
                t=1;
            end
            if this.useMex
                Y = HnExp(X,V,t);                
            else
                Y = this.localExp(X,V,t);
            end
        end
        function V = log(this,X,Y)
            % log(X,Y) - Exponential map at the point(s) X of points(s) Y
            %      in TXHn
            %
            % INPUT
            %   X : a point or set of points on the manifold Hn
            %   Y : a point or set of points on the manifold Hn
            %
            % OUTPUT
            %   V : resulting point(s) of X(i) to Y(i) elementwise
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1, R. Bergmann ~ 2015-10-20

            % computes v = (y+x*<x,y>)/sqrt(<x,y>^-1)*d(x,y)
            if this.useMex
                V = HnLog(X,Y);
            else
                V = this.localLog(X,Y);
            end
        end
        function W = parallelTransport(this,X,Y,V)
                % Parallel Transport the vector V from TxM to TyM
                %
                % directional part that changes
                sV = size(V);
                dir = this.log(X,Y);
                norm_dir = repmat(sqrt(sum(dir.^2,1)),[this.ItemSize,ones(1,length(sV(2:end)))]);
                dir = dir./norm_dir;
                scp = sum(dir.*V,1);
                % substract V-part (scp*dir) and add the negative direction
                W = V - repmat(scp,[this.ItemSize,ones(1,sV(2:end))]).*(dir+this.log(Y,X)./norm_dir);
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
            % Manifold-Valued Image Restoration Toolbox 1.1, J. Persch 2015-11-05
           
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
                if all(dims(1:length(this.ItemSize)) == this.ItemSize) && length(dims)<length(this.ItemSize)+2
                    x = f;
                    return
                end
                error('f wrong size');
            end
            % shift manDim in first dimension
            dim = size(f);
            m = dim(end-1);
            n = dim(end);
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
            if ~this.useMex
                x=this.mean@manifold(f,'Weights',w,'Epsilon',epsilon,'MaxIterations',iter);
            else            
            x = HnMean(f,w,epsilon,iter);
            end            
        end
        function d = dist(this,X,Y)
            % d = dist(X,Y) - Compute the distance between points or
            % two sets of points on the manifold H(n).
            %
            % INPUT
            %   X,Y : two points (matrices) or sets of points (matrices)
            %         on P(m)
            %
            % OUTPUT
            %     d : resulting distances of each pair of points of p,q.
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1, R. Bergmann ~ 2014-10-19 | 2015-04-11
            if this.useMex
                d = HnDist(X,Y);
            else
                d = this.localDist(X,Y);
            end
        end
        function ds = dot(this,X,V,W)
            % dot(V,W)
            %     Compute the inner product of two tangent vectors in TXHn
            %
            % INPUT
            %     X  : base point (optional because all TXM are equal)
            %     V  : a first tangent vector( set)
            %     W  : a secod tangent vector( set)
            %
            % OUTPUT
            %     ds : the corresponding value(s) of the inner product of (each triple) V,W at X
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1, R. Bergmann ~ 2015-10-20
            if nargin < 4 
                W = V;
                V = X;
            end
            if this.useMex
                ds = HnDot(V,W);
            else
                ds = this.localDot(V,W);
            end
        end
        function Y = addNoise(~,X,sigma)
            sizes = size(X);
            error('This method is not yet implemented.');
        end
    end
    methods (Access = private)
        function Y = localExp(this,X,V,t)
            dims = size(X);
            X_ = reshape(X,this.ItemSize,[]);
            V_ = reshape(V,this.ItemSize,[]);
            normv = sqrt(max(this.dot(V_,V_),0)).';
            Y = X_;
            % compute Y = cosh(t*norm(v))x + sinh(t*norm(v))*v/norm(v))
            % for all terms, where the norm is large enugh, otherwise, stay
            % at x.
            if any(normv>eps)
                Y(:,normv>eps) = X_(:,normv>eps).*repmat(cosh(t*normv(normv>eps)),[this.ItemSize,1])...
                    + V_(:,normv>eps).*repmat(sinh(t*normv(normv>eps))./normv(normv>eps),[this.ItemSize,1]);
            end
            % be sure to reshape again onto Hn for security reasons
            Y = Y./(sqrt(-repmat(permute(this.dot(Y,Y),[3,1,2]),this.ItemSize,1)));
            Y = reshape(Y,dims);
        end
        function V = localLog(this,X,Y)
            dims = size(X);
                X_ = reshape(X,this.ItemSize,[]);
                Y_ = reshape(Y,this.ItemSize,[]);
            scp = this.dot(X_,Y_).';
            w = Y_ + X_.*repmat(scp,[this.ItemSize,1]);
            normw = sqrt(max(this.dot(X_,Y_).^2-1,0)).';
            V = zeros(size(X_));
            scp = acosh(max(1,-scp)); % = d(x,y)
            % only apply to nonzero distances
            cols = permute((normw>eps),[2:length(size(Y)),1]);
            if any(cols)
                V(:, cols ) = w(:, cols )./repmat(...
                    normw(1,cols)...
                    ,[this.ItemSize,ones(1,length(size(Y))-1)]).*...
                    repmat(...
                    scp(1,cols)...
                    ,[this.ItemSize,ones(1,length(size(Y))-1)]);
            end
            V = reshape(V,dims);
        end
        function d = localDist(this,X,Y)
            dims = size(X);
            X = reshape(X,this.ItemSize,[]);
            Y = reshape(Y,this.ItemSize,[]);
            d = zeros(prod(dims(2:end)),1);
            cols = ~(max(abs(X-Y))<eps);
            if any(cols)
                d(cols) = acosh(max(1,-this.dot(X(:,cols),Y(:,cols))));  
            end
            if ~isscalar(dims(2:end))
                d = reshape(d,dims(2:end));
            end
        end
        function ds = localDot(this,V,W)
            dims = size(V);
            if length(dims)==2
                ds =  ([ones(1,this.ItemSize-1),-1]*(V.*W)).';
            else
                ds = reshape(...
                     [ones(1,this.ItemSize-1),-1]*...
                         (reshape(V,this.ItemSize,[]).*reshape(W,this.ItemSize,[]))...
                     ,dims(2:end));
            end
        end
    end
end