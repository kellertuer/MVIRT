classdef Sn < manifold & handle
    % The n-dimensional sphere embedded in R(n+1)
    %
    % Additionally provides
    %
    % PROPERTIES
    %    tau   : stepsize in the subgradient descent in the proximal
    %            mappings of the second order differences
    %    steps : number of steps in the just mntioned subgradient descent.
    %
    % FUNCTIONS
    %    midPoint(x,z) : Compute the mid point of x,z on S2
    %    TpMONB(p,q)   : Compute an ONB of the tangential plane at p, where
    %                 the first vector points towards q
    % ---
    % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-01-29
    properties
        type = '';
        ItemSize;
        Dimension;
        allDims;
    end

    methods
        function obj = Sn(n)
            obj.ItemSize = n+1;
            obj.Dimension = n;
            obj.type = ['The ',num2str(n),'-sphere in R',num2str(n+1)];
            obj.allDims = repelem({':'},length(obj.ItemSize));
        end
        function y = exp(this,x,xi,t)
            % exp(x,xi) - Exponential map at x with direction xi in TxM
            %
            % INPUT
            %   x : a point or set of points on the manifold S2
            %   xi : a point or set of point in the tangential spaces TxM
            %
            % OPTIONAL:
            %   t : [] following the geodesics for one time step
            %       a scalar following the geodesics for time t
            %       a set of t following each geodesic_i for time t_i
            %
            % OUTPUT
            %   y : resulting point(s) on S2
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            if isrow(x)
                p_=x';
            else
                p_=x;
            end
            if isrow(xi)
                v_=xi';
            else
                v_=xi;
            end
            if nargin < 4
                t=1;
            elseif isrow(t)
                t = t.';
            end
            if this.useMex
                y = SnExp(p_,v_,t);
            else
                y = this.localExp(p_,v_,t);
            end
        end
        function xi = log(this,x,y)
            % log(x,y) - Inverse Exponential Map at x of y.
            %
            % INPUT
            %    x : point or set of (column) points indicating the
            %    tangential base points
            %    y : point(s) on S2 being put into the tangential plane at
            %    their corresponding p
            %
            % OUTPUT
            %    xi : points on the tangential plane at point(s) p
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            if isrow(x)
                p_=x';
            else
                p_=x;
            end
            if isrow(y)
                q_=y';
            else
                q_=y;
            end
            if this.useMex
                xi = SnLog(p_,q_);
            else
                xi = this.localLog(p_,q_);
            end
        end
        function d = dist(this,x,y)
            % dist(p,q) - Compute the distance between points or a set of
            % points on the manifold S2.
            %
            % INPUT
            %   x,y : a (column) vector from S2 (embd. in R3) or a set of
            %   column vectors
            %
            % OUTPUT
            %     d : resulting distances of each column pair of p,q.
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19 | 2015-03-30

            % Lofgile
            %   2015-03-30 Changed dist to work to collapse first dim
            %   2015-04-11 Extracted to Mex
            if this.useMex
                d = SnDist(x,y);
            else
                d = this.localDist(x,y);
            end
        end
        function fn = addNoise(this,f,sigma)
            fs = size(f);
            if fs(1) ~= this.ItemSize
                error('Wrong data size');
            end
            f = reshape(f,this.ItemSize,[]);
            noise = sigma*randn(fs);
            fn = zeros(this.ItemSize,prod(fs(2:end)));
            for i=1:prod(fs(2:end))
                fn(:,i) = this.exp(f(:,i),noise(:,i) - (permute(noise(:,i),[2,1])*f(:,i))*f(:,i));
            end
            fn = reshape(fn,fs);
        end
        function eta = parallelTransport(this,x,y,xi)
        % eta = parallelTransport(x,y,xi) transport xi from TxM parallel to TyM
        %
        % INPUT
        %   x,y : two (sets of) points on the manifold
        %   xi  : a (set of) vectors from TxM
        %
        % OUTPUT
        %   eta : the parallel transported vectors in TyM
        % ---
        % Manifold-valued Image Restoration Toolbox 1.2
        % R. Bergmann | 2018-03-01
            if this.useMex
                eta = SnParallelTransport(x,y,xi);
            else
                % Parallel Transport the vector V from TxM to TyM
                %
                % directional part that changes
                sXi = size(xi);
                dir = this.log(x,y);
                norms = sqrt(sum(dir.^2,1));
                norm_dir = repmat(norms,[this.ItemSize,ones(1,length(sXi(2:end)))]);
                normMask = repmat(norms==0,[this.ItemSize,ones(1,length(sXi(2:end)))]);
                dir = dir./norm_dir;
                scp = sum(dir.*xi,1);
                % substract V-part (scp*dir) and add the negative direction
                eta = xi - repmat(scp,[this.ItemSize,ones(1,length(sXi(2:end)))]).*(dir+this.log(y,x)./norm_dir);
                eta(normMask) = xi(normMask); %those that need not to be transported
            end
        end
        function [xi,k] = TpMONB(this,x,y)
        % xi = TpMONB(x,y)
        % Compute an ONB in TpM, where the first vector points to q,
        % whin given.
        %
        % INPUT
        %     x : base point( sets)
        %
        % OPTIONAL:
        %     y : directional indicator( sets) for the first vector(s).
        %
        % OUTPUT
        %    xi : basiscolumn matrice(s) for each TpM
        %    k : (optional) curvature coefficients, here
        %           all are 1 except the first which is 0
        % ---
        % MVIRT 1.0, R. Bergmann ~ 2014-10-19 | 2018-03-01
            if isrow(x)
                p_=x';
            else
                p_=x;
            end
            q_given = 0;
            if nargin < 3
                y=x;
            else
                q_given=1;
            end
            pS = size(x);
            p_ = reshape(p_,pS(1),[]);
            q_ = reshape(y, pS(1),[]);
            if q_given && max(this.dist(q_,p_)) > eps
                xi = zeros(this.ItemSize,size(p_,2),this.ItemSize-1);
                xi(:,:,1) = this.log(p_,q_);
                normsv = sqrt(sum(xi(:,:,1).^2,1));
                if ~all(normsv<=eps)
                    xi(:,normsv>eps,1) = xi(:,normsv>eps,1)./repmat(normsv(normsv>eps),[this.ItemSize,1,1]);
                end
            else
                xi = zeros(this.ItemSize,size(p_,2),this.ItemSize-1);
            end
            if this.ItemSize==3 && q_given == 1%S2 -> cross
                xi(:,:,2) = cross( squeeze(xi(:,:,1)), p_);
            else
                for col=1:size(xi,2)
                    if ~q_given || (normsv(col)>eps)
                        % The remaining Tangential vectors are the orthogonal
                        % to p(:,col) and V(:,col,1), i.e. the nullspace of the
                        % matrix p V V ... V
                        xi(:,col,1+q_given:this.Dimension) = null([p_(:,col), xi(:,col,1)].');
                    end
                end
            end
            xi = reshape(xi,[pS,this.ItemSize-1]);
            if nargout > 1
                k = ones(size(p_,2),this.ItemSize-1);
                k(:,1)=0;
                k = reshape(k,[pS(2:end),this.ItemSize-1]);
            end
        end
        function ds = dot(~,x,xi,nu)
            % dot(x,xi,nu)
            %     Compute the inner product of two tangent vectors in T_xM
            %
            % INPUT
            %     x  : a point(Set) in P(n)
            %     xi  : a first tangent vector( set) to (each) x
            %     nu  : a secod tangent vector( set) to (each) x
            %
            % OUTPUT
            %     ds : the corresponding value(s) of the inner product of
            %     (each triple) xi,nu at x
            %
            % ---
            % MVIRT 1.0 ~ J. Persch 2016-06-13
            dimen = size(x);
            if all(size(xi) == dimen & size(nu) == dimen)
                ds = permute(sum(xi.*nu,1),[2:length(dimen),1]);
            else
                error('Dimensions of Input must coincide')
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
            % MVIRT 1.0, J. Persch 2015-11-05
            if ~this.useMex %just pass further on
                x=this.mean@manifold(varargin{:});
            else
            ip = inputParser;
            addRequired(ip,'f');
            addParameter(ip,'Weights',[]);
            addParameter(ip,'InitVal',[]);
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
            if isempty(vars.Weights)
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
            x = SnMean(f,w,epsilon,iter);
            end
        end
    end

    methods (Access = protected)
        function d = localDist(~,p,q)
            if ( size(p,2) > 1) && (size(q,2) > 1) %both nonsingleton
                assert(all(size(p)==size(q)),'p and q have to be of same size');
            end
            d = shiftdim(acos( min( max(sum(p.*q,1), -1), 1) ),1); %remove first (!) singleton dimension
        end
        function q = localExp(this,p_,v_,t)
            if nargin < 4
                t=1;
            elseif ~isscalar(t)
                sP = size(p_);
                sP = sP(2:end);
                sT = size(t);
                sT = sT(sT>1);
                if any(sP(2:end) ~= sT)
                    error('t needs to be of same dimension as p_ and v_ or a scalar')
                end
            end
            if (size(p_,1)~=this.ItemSize) || (size(v_,1)~=this.ItemSize)
                error(['Both sets of vectors have to be vectors from R',num2str(this.ItemSize),'.']);
            end
            if any(size(p_)~=size(v_))
                error(['The number of points p (',num2str(size(p_,2)),...
                    ') and tangential vectors v (',num2str(size(v_,2))...
                    ,')is different.']);
            end
            % computes y = cos(t*norm(v)) + sin(t*norm(v))*v/norm(v);
            %Compute columwise norms of v
            normv = sqrt(sum(v_.^2,1));
            q = p_;
            % only the nonzero vs are directional
            if any(normv>eps)
                % if to save the repmat
                if isscalar(t)
                    q(:,normv>eps) = p_(:,normv>eps).*repmat(cos(t*normv(:,normv>eps)),[this.ItemSize,1])...
                        + v_(:,normv>eps).*repmat(sin(t*normv(:,normv>eps))./normv(:,normv>eps),[this.ItemSize,1]);
                else
                    q(:,normv>eps) = p_(:,normv>eps).*repmat(cos(t(normv>eps).'.*normv(:,normv>eps)),[this.ItemSize,1])...
                        + v_(:,normv>eps).*repmat(sin(t(normv>eps).'.*normv(:,normv>eps))./normv(:,normv>eps),[this.ItemSize,1]);
                end
            end
        end
        function v = localLog(this,p_,q_)
            % Fallback is no mex available
            if any(all(size(p_)~=size(q_)))
                error('p and q have to be of same size');
            end
            dimen = size(p_);
            p_ = reshape(p_,this.ItemSize,[]);
            q_ = reshape(q_,this.ItemSize,[]);
            % Compute all scalar products of points p and q, secure, that
            % they are in -1,1 for numerical reasons
            % computes v = y-<y,x>x/norm(y-<y,x>x) * d(x,y)
            scp = min( max(sum(p_.*q_,1), -1), 1);
            w = q_ - p_.*repmat(scp,[this.ItemSize,1]);
            normw = sqrt(sum(w.^2,1));
            v = zeros(size(p_));
            scp = acos(scp);
            % only apply to nonzero distances
            cols = permute((normw>eps),[2:length(size(q_)),1]);
            if any(cols)
                v(:, cols ) = w(:, cols )./repmat(...
                    normw(1,cols)./scp(1,cols)...
                    ,[this.ItemSize,ones(1,length(size(q_))-1)]);
            end
            v = reshape(v,dimen);
        end
    end
end
