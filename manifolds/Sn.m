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
        tau = 0.01;
        steps = 10;
        one = NaN;
        ItemSize;
        Dimension;
        allDims;
    end
    
    methods
        function obj = Sn(n)
            obj.ItemSize = n+1;
            obj.Dimension = n;
            obj.type = ['The ',num2str(n),'-sphere in R',num2str(n+1)];
            if n == 3
                obj.one = [1;0;0;0];
            end
            obj.allDims = repelem({':'},length(obj.ItemSize));
        end
        function q = exp(this,p,v,t)
            % exp(p,v) - Exponential map at the point p with respect to v in
            % TpM
            %
            % INPUT
            %   p : a point or set of points on the manifold S2
            %   v : a point or set of point in the tangential spaces TpM
            % OPTIONAL:
            %   t : [] following the geodesics for one time step
            %       a scalar following the geodesics for time t
            %       a set of t following each geodesic_i for time t_i
            %
            % OUTPUT
            %   q : resulting point(s) on S2
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            if isrow(p)
                p_=p';
            else
                p_=p;
            end
            if isrow(v)
                v_=v';
            else
                v_=v;
            end
            if nargin < 4
                t=1;
            elseif isrow(t)
                t = t.';
            end
            if this.useMex
                q = SnExp(p_,v_,t);
            else
                q = this.localExp(p_,v_,t);
            end
        end
        function v = log(this,p,q)
            % log(q,p) - Inverse Exponential Map at p of q.
            %
            % INPUT
            %    p : point or set of (column) points indicating the
            %    tangential base points
            %    q : point(s) on S2 being put into the tangential plane at
            %    their corresponding p
            %
            % OUTPUT
            %    v : points on the tangential plane at point(s) p
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            if isrow(p)
                p_=p';
            else
                p_=p;
            end
            if isrow(q)
                q_=q';
            else
                q_=q;
            end
            if this.useMex
                v = SnLog(p_,q_);
            else
                v = this.localLog(p_,q_);
            end
        end
        function d = dist(this,p,q)
            % dist(p,q) - Compute the distance between points or a set of
            % points on the manifold S2.
            %
            % INPUT
            %   p,q : a (column) vector from S2 (embd. in R3) or a set of
            %   column vectors
            %
            % OUTPUT
            %     v : resulting distances of each column pair of p,q.
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19 | 2015-03-30
            
            % Lofgile
            %   2015-03-30 Changed dist to work to collapse first dim
            %   2015-04-11 Extracted to Mex
            if this.useMex
                d = SnDist(p,q);
            else
                d = this.localDist(p,q);
            end
        end
        function G = grad_X_D2(this,X,Y,Z)
            % grad_X_D2_sq(X,Z,Y) Compute the gradient with respect to the first
            % variable of the second order difference term
            % d(c(X,Z),Y). This can also
            % be used for the third, Z, exchanging X and Z
            %
            % INPUT
            %   X : A point( set)
            %   Y : A point( set)
            %   Z : A point( set)
            %
            % OUTPUT
            %   G : A (set of) tangent vector(s, one) at (each) X.
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.1, J. Persch  2016-09-21
            if this.useMex
                G = SnGrad_X_D2(X,Y,Z);
            else
                G = this.localGrad_X_D2(X,Y,Z);
            end
        end
        function G = grad_X_D2_Sq(this,X,Y,Z)
            % grad_X_D2_sq(X,Z,Y) Compute the gradient with respect to the first
            % variable of the squared second order difference term
            % d^2(c(X,Z),Y). This can also
            % be used for the third, Z, exchanging X and Z
            %
            % INPUT
            %   X : A point( set)
            %   Y : A point( set)
            %   Z : A point( set)
            %
            % OUTPUT
            %   G : A (set of) tangent vector(s, one) at (each) X.
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch  2016-09-21
            if this.useMex
                G = SnGrad_X_D2_Sq(X,Y,Z);
            else
                G = this.localGrad_X_D2_Sq(X,Y,Z);
            end
        end
        function [x,y] = prox_midpoint_data_sq(this,varargin)
            % prox_midpoint_data_sq(v,w,f, lambda) - Compute the proximal mapping of the
            % images v,w to the data f with respect to lambda.
            %
            % INPUT
            %      v,w,l : data on S(m) (i.e. m), where [m,l] = size(f) = size(v) = size(w)
            %          and l is the number of parallel proxes.
            % lambda : weight of the proximal map (one value)
            %
            % OUTPUT
            %      x : resulting data of all proximal maps (mxmxl)
            % ---
            % ManImRes 1.0, J. Persch ~ 2016-09-21
            if (length(varargin)==1 && isstruct(varargin{1})) %struct case
                vars = varargin{1};
                assert(isfield(vars,'f')&&isfield(vars,'lambda')&&isfield(vars,'v')&&isfield(vars,'w'),...
                    'Not all required parameters given in the struct');
                if ~isfield(vars,'RegMask')
                    vars.RegMask = [];
                end
            else
                ip = inputParser;
                addRequired(ip,'v');
                addRequired(ip,'w');
                addRequired(ip,'f');
                addRequired(ip,'lambda');
                parse(ip, varargin{:});
                vars = ip.Results;
                [a,~] = size(vars.f);
                if (a~=this.ItemSize)
                    error(['The values of f are not ',num2str(this.ItemSize),' vectors.']);
                end
                if any(size(vars.v)~=size(vars.w)) || any(size(vars.v)~=size(vars.f))
                    error('Inputs v,w,f should be of the same size!');
                end
            end
            v = vars.v;
            x = v;
            w = vars.w;
            y = w;
            if (vars.lambda==0)% no prox
                return
            end
            f = vars.f;
            % Gradient descent
            xopt = x;
            yopt = y;
            func_optvals =  vars.lambda*this.dist(this.midPoint(x,y),f).^2;
            for i=1:this.steps
                % midpoints between firsts and thirds
                % Gradient X
                grad1 = this.log(x,v) + vars.lambda*this.grad_X_D2_Sq(x,f,y);
                % Gradient Y
                grad2 = this.log(y,w) + vars.lambda*this.grad_X_D2_Sq(y,f,x);
                tauit = this.tau/i;
                x = this.exp(x, tauit*grad1);
                y = this.exp(y, tauit*grad2);
                % Testing if new values are better
                xoptt = x;
                yoptt = y;
                func_vals = 1/2*this.dist(v,xoptt).^2+1/2*this.dist(w,yoptt).^2 +...
                    vars.lambda*this.dist(this.midPoint(xoptt,yoptt),f).^2;
                xopt(:,func_vals<func_optvals) = xoptt(:,func_vals<func_optvals);
                yopt(:,func_vals<func_optvals) = yoptt(:,func_vals<func_optvals);
                y = yopt;
                x = xopt;
                func_optvals(func_vals<func_optvals) = func_vals(func_vals<func_optvals);
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
        function W = parallelTransport(this,X,Y,V)
            if this.useMex
                W = SnParallelTransport(X,Y,V);
            else
                % Parallel Transport the vector V from TxM to TyM
                %
                % directional part that changes
                sV = size(V);
                dir = this.log(X,Y);
                norms = sqrt(sum(dir.^2,1));
                norm_dir = repmat(norms,[this.ItemSize,ones(1,length(sV(2:end)))]);
                normMask = repmat(norms==0,[this.ItemSize,ones(1,length(sV(2:end)))]);
                dir = dir./norm_dir;
                scp = sum(dir.*V,1);
                % substract V-part (scp*dir) and add the negative direction
                W = V - repmat(scp,[this.ItemSize,ones(1,length(sV(2:end)))]).*(dir+this.log(Y,X)./norm_dir);
                W(normMask) = V(normMask); %those that need not to be transported
            end
        end
        function [V,k] = TpMONB(this,p,q)
        % V = TpMONB(p,q)
        % Compute an ONB in TpM, where the first vector points to q,
        % whin given.
        %
        % INPUT
        %     p : base point( sets)
        %
        % OPTIONAL:
        %     q : directional indicator( sets) for the first vector(s).
        %
        % OUTPUT
        %    V : basiscolumn matrice(s)
        %    k : (optional) curvature coefficients, here
        %           all are 1 except the first which is 0
        % ---
        % MVIRT 1.0, R. Bergmann ~ 2014-10-19 | 2014-10-23
            if isrow(p)
                p_=p';
            else
                p_=p;
            end
            q_given = 0;
            if nargin < 3
                q=p;
            else
                q_given=1;
            end
            pS = size(p);
            p_ = reshape(p_,pS(1),[]);
            q_ = reshape(q, pS(1),[]);
            if q_given && max(this.dist(q_,p_)) > eps
                V = zeros(this.ItemSize,size(p_,2),this.ItemSize-1);
                V(:,:,1) = this.log(p_,q_);
                normsv = sqrt(sum(V(:,:,1).^2,1));
                if ~all(normsv<=eps)
                    V(:,normsv>eps,1) = V(:,normsv>eps,1)./repmat(normsv(normsv>eps),[this.ItemSize,1,1]);
                end
            else
                V = zeros(this.ItemSize,size(p_,2),this.ItemSize-1);
            end
            if this.ItemSize==3 && q_given == 1%S2 -> cross
                V(:,:,2) = cross( squeeze(V(:,:,1)), p_);
            else
                for col=1:size(V,2)
                    if ~q_given || (normsv(col)>eps)
                        % The remaining Tangential vectors are the orthogonal
                        % to p(:,col) and V(:,col,1), i.e. the nullspace of the
                        % matrix p V V ... V
                        V(:,col,1+q_given:this.Dimension) = null([p_(:,col), V(:,col,1)].');
                    end
                end
            end
            V = reshape(V,[pS,this.ItemSize-1]);
            if nargout > 1
                k = ones(size(p_,2),this.ItemSize-1);
                k(:,1)=0;
                k = reshape(k,[pS(2:end),this.ItemSize-1]);
            end
        end
        function ds = dot(~,P,V,W)
            % Sn.dot(P,V,W)
            %     Compute the inner product of two tangent vectors in T_P M
            %
            % INPUT
            %     X  : a point(Set) in P(n)
            %     V  : a first tangent vector( set) to (each) X
            %     W  : a secod tangent vector( set) to (each) X
            %
            % OUTPUT
            %     ds : the corresponding value(s) of the inner product of (each triple) V,W at X
            %
            % ---
            % MVIRT 1.0 ~ J. Persch 2016-06-13
            dimen = size(P);
            if all(size(V) == dimen & size(W) == dimen)
                ds = permute(sum(V.*W,1),[2:length(dimen),1]);
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
                x = SnMean(f,w,epsilon,iter);
            end
        end
        function W = geopoint(this,X,Y,t)
            % geopoint(X,Y,t) - Give the point \gamma_{XY}(t)
            %
            % INPUT
            %   X,Y : a point or set of points on the manifold Sn(n)
            %   t : a scalar or set of scalars
            %
            %
            % OUTPUT
            %   W : resulting point(s) on Sn(n)
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch ~ 2017-03-31
            
            % Changelog
            %   2015-04-10 introduced mexfile
            W = this.exp(X,this.log(X,Y),t);
        end
        function s = inv(this,s)
            % S1.inv(s) inverts s in the sense of rotations about that angle
            if this.Dimension == 3                
                s(2:4,:) = -s(2:4,:);
            else 
                error('No group structure implemented')
            end
        end
         function s = ga(this,s_1,s_2)
            % S1.ga(s_1,s_2) performs the group action s_1 o s_2
            s = s_1;
            if this.Dimension == 3                
                s(1,:) = s_1(1,:).*s_2(1,:)-sum(s_1(2:4,:).*s_2(2:4,:),1);
                s(2:4,:) = s_1(1,:).*s_2(2:4,:)+s_2(1,:).*s_1(2:4,:)+cross(s_1(2:4,:),s_2(2:4,:));
            else 
                error('No group structure implemented')
            end
         end
         function sv = DL_s(this,s,v)
            % Sn.DL_s(s_1,s_2) performs the derivative of the left group
            % action of s onto v \in T_t S1, which goes into sv\in\T_st S1
            if this.Dimension == 3                
            	sv = this.ga(s,v);
            else 
                error('No group structure implemented')
            end
         end
         function vs = DR_s(this,s,v)
            % S1.DR_s(s_1,s_2) performs the derivative of the right group
            % action of s onto v \in T_t S1, which goes into vs\in\T_ts S1
            if this.Dimension == 3                
            vs = this.ga(v,s);
            else 
                error('No group structure implemented')
            end
        end
    end
    
    methods (Access = protected)
        function G = localGrad_X_D2(this,X,Y,Z)
            l = size(X,2);
            m = this.exp(X, this.log(X,Z)./2);
            % Common directions to middle points y
            r = this.log(m,Y);
            normr = permute(sqrt(sum(r.^2,1)),[2:length(size(r)),1]);
            % Grad X
            V = this.TpMONB(m,Z);
            % Instead of
            % W = this.TpMONB(X,Z);
            % we have would have to use PT
            % W = this.ParallelTransport(m,X,V);
            % but we can also just adapt the first column, because the
            % rest stays as it is
            W = V;
            W(:,:,1) = this.log(X,Z);
            normsw = sqrt(sum(W(:,:,1).^2,1));
            W(:,normsw>eps,1) = W(:,normsw>eps,1)./repmat(normsw(normsw>eps),[this.ItemSize,1,1]);
            alpha = zeros(l,this.ItemSize-1);
            if any(normr>eps)
                alpha(normr>eps,:) = sum(...
                    repmat(r(:,normr>eps),[1,1,this.ItemSize-1]).*...
                    V(:,normr>eps,:)./repmat(permute(normr(normr>eps),[2,1]),[this.ItemSize,1,this.ItemSize-1]),1);
            end
            % and even inplace, i.e. we rename the basis at M (W) also
            % to V
            G = W(:,:,1).*repmat(0.5*permute(alpha(:,1),[2,1]),[this.ItemSize,1]) ...
                + sum(W(:,:,2:this.ItemSize-1).*...
                repmat(...
                repmat(permute( 1./(2*cos(this.dist(X,Z)./2)) ,[3,1,2]),[1,1,this.ItemSize-2]).*...
                permute(alpha(:,2:this.ItemSize-1),[3,1,2]),...
                [this.ItemSize,1,1]),3);
        end
        function G = localGrad_X_D2_Sq(this,X,Y,Z)
            dimen = size(X);
            X = reshape(X,this.ItemSize,[]);
            Y = reshape(Y,this.ItemSize,[]);
            Z = reshape(Z,this.ItemSize,[]);
            l = size(X,2);
            m = this.exp(X, this.log(X,Z)./2);
            % Common directions to middle points y
            r = this.log(m,Y);
            normr = permute(sqrt(sum(r.^2,1)),[2:length(size(r)),1]);
            % Grad X
            V = this.TpMONB(m,Z);
            % Instead of
            % W = this.TpMONB(X,Z);
            % we have would have to use PT
            % W = this.ParallelTransport(m,X,V);
            % but we can also just adapt the first column, because the
            % rest stays as it is
            W = V;
            W(:,:,1) = this.log(X,Z);
            normsw = sqrt(sum(W(:,:,1).^2,1));
            W(:,normsw>eps,1) = W(:,normsw>eps,1)./repmat(normsw(normsw>eps),[this.ItemSize,1,1]);
            alpha = zeros(l,this.ItemSize-1);
            if any(normr>eps)
                alpha(normr>eps,:) = sum(...
                    2*repmat(r(:,normr>eps),[1,1,this.ItemSize-1]).*...
                    V(:,normr>eps,:),1);
            end
            % and even inplace, i.e. we rename the basis at M (W) also
            % to V
            G = W(:,:,1).*repmat(0.5*permute(alpha(:,1),[2,1]),[this.ItemSize,1]) ...
                + sum(W(:,:,2:this.ItemSize-1).*...
                repmat(...
                repmat(permute( 1./(2*cos(this.dist(X,Z)./2)) ,[3,1,2]),[1,1,this.ItemSize-2]).*...
                permute(alpha(:,2:this.ItemSize-1),[3,1,2]),...
                [this.ItemSize,1,1]),3);
            G = reshape(G,dimen);
        end
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
