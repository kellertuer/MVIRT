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
            % Create a Manifold for thh n-dimensional hyperbolic manifold
            obj.type = ['hyperbolic manifold of dimension ',num2str(n)];
            obj.ItemSize = n+1;
            obj.Dimension = n;
            obj.allDims = repelem({':'},length(obj.ItemSize));
        end
        function y = exp(this,x,xi,t)
            % exp(x,xi) exponential map at the point(s) x towards xi in T_xHn
            %
            % INPUT
            %   x : a point or set of points on the manifold Hn
            %  xi : a point or set of point in the tangential spaces TxHn
            %
            % OPTIONAL
            %   t : shorten vectors xi by factor t
            %   (given as one value or an array of same length as number of xis)
            %
            % OUTPUT
            %   y : resulting point(s) on Hn
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.
            %  R. Bergmann ~ 2015-01-20 | 2015-10-20
            if nargin<4 % define t
                t=1;
            end
            if this.useMex
                y = HnExp(x,xi,t);
            else
                y = this.localExp(x,xi,t);
            end
        end
        function xi = log(this,x,y)
            % xi = log(x,y) logarithmic map the point(s) x from y
            %
            % INPUT
            %   x : a point or set of points on the manifold Hn
            %   y : a point or set of points on the manifold Hn
            %
            % OUTPUT
            %   xi : resulting point(s) of x(i) to y(i) elementwise
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1
            % R. Bergmann ~ 2015-10-20
            
            % computes v = (y+x*<x,y>)/sqrt(<x,y>^-1)*d(x,y)
            if this.useMex
                xi = HnLog(x,y);
            else
                xi = this.localLog(x,y);
            end
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
        function d = dist(this,x,y)
            % d = dist(x,y) distance between points x and y on Hn
            %
            % INPUT
            %   x,y : two points or sets of points on Hn
            %
            % OUTPUT
            %     d : resulting distances of each pair of points of p,q.
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1
            % R. Bergmann ~ 2014-10-19 | 2015-04-11
            if this.useMex
                d = HnDist(x,y);
            else
                d = this.localDist(x,y);
            end
        end
        function ds = dot(this,x,xi,nu)
            % ds = dot(x,xi,nu) inner product of two tangent vectors in TxHn
            %
            % INPUT
            %     x  : base point (optional because all TXM are equal)
            %    xi  : a first tangent vector( set)
            %    nu  : a secod tangent vector( set)
            %
            % OUTPUT
            %     ds : the corresponding value(s) of the inner product of (each triple) V,W at X
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1
            % R. Bergmann ~ 2015-10-20
            if nargin < 4
                nu = xi;
                xi = x;
            end
            if this.useMex
                ds = HnDot(xi,nu);
            else
                ds = this.localDot(x,xi,nu);
            end
        end
        function eta = parallelTransport(this,x,y,xi)
            % eta = parallelTransport(x,y,xi) parallel transport xi along g(.,x,y)
            %
            % INPUT
            %    x : a point( set) on P(m)
            %    y : a point( set) on P(m)
            %   xi : a tangential vector( set, one) at (each) X
            %
            % OUTPUT
            %  eta : tangential vector( set, one) at (each) Y
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.
            %  R. Bergmann ~ 2018-03-05 | 2018-03-05
            
            % Note that this is very similar -- again 
            sXi = size(xi);
            dir = this.log(x,y);
            norms = sqrt(this.dot(x,dir,dir));
            dir = dir./shiftdim(norms,-1);
            normMask = repmat(norms==0,[this.ItemSize,ones(1,length(sXi(2:end)))]);
            scp = this.dot(x,dir,xi);
            % substract V-part (scp*dir) and add the negative direction
            eta = xi - shiftdim(scp,-1).*(dir+this.log(y,x)./shiftdim(norms,-1));
            eta(normMask) = xi(normMask); %those that need not to be transported
        end
        function Xi = TpMONB(this,x,y)
            % Xi = TpMONB(x,y)
            % Compute an ONB in TpM, where the first vector points to y,
            % whin given.
            %
            % INPUT
            %     x : base point( sets)
            % OPTIONAL:
            %     y : directional indicator( sets) for the first vector(s).
            %
            % OUTPUT
            %    Xi : basiscolumn matrice(s)
            % ---
            % MVIRT 1.0, R. Bergmann ~ 2014-10-19 | 2014-10-23
            if isrow(x)
                x_=x';
            else
                x_=x;
            end
            y_given = 0;
            id = eye(this.ItemSize);
            Xi = permute(repmat(id(:,1:end-1),1,1,size(x_,2)),[1,3,2]);
            if nargin == 3
                y_given =1;
                if isrow(y)
                    y_=y';
                else
                    y_=y;
                end
                Xi(:,:,1) = this.log(x_,y_);
                normsv = sqrt(sum(Xi(:,:,1).^2,1));
                Xi(:,normsv>eps,1) = Xi(:,normsv>eps,1)./repmat(normsv(normsv>eps),[this.ItemSize,1,1]);
            end
            for col=1:size(Xi,2)
                for i = 1+y_given:this.Dimension
                    Xi(:,:,i) = Xi(:,:,i) - x.*repmat(permute(this.dot(Xi(:,:,i),x,Xi(:,:,i))./this.dot(x,x,x),[2,1]),this.ItemSize,1,1);
                    for j = 1:i-1
                        Xi(:,:,i) = Xi(:,:,i) - Xi(:,:,j).*repmat(permute(this.dot(Xi(:,:,i),Xi(:,:,j),Xi(:,:,i)),[2,1]),this.ItemSize,1,1);
                    end
                    Xi(:,:,i) = Xi(:,:,i)./sqrt(repmat(permute(this.dot(Xi(:,:,i),Xi(:,:,i),Xi(:,:,i)),[2,1]),this.ItemSize,1,1));
                end
            end
        end
        function y = addNoise(this,x,sigma)
        % y = addNoise(x,sigma) add noise 
            sizes = size(x);
            x = reshape(x,this.ItemSize,[]);
            ONB = this.TpMONB(x);
            y = this.exp(x,sum(repmat(sigma*randn(1,size(x,2),this.Dimension),this.ItemSize,1,1).*ONB,3));
            y = reshape(y,sizes);
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
        function ds = localDot(this,~,xi,nu)
            dims = size(xi);
            dataDims = dims(2:end);
            p = [ones(this.ItemSize-1,1);-1]; %prefactors in inner product
            ds = reshape(...
               sum(p.*(xi(this.allDims{:},:).*nu(this.allDims{:},:)),1),...
                [dataDims,1]);
        end
    end
end