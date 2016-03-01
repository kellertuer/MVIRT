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
    % ManImRes 1.0, R. Bergmann ~ 2015-01-29
    properties
        type = '';
        tau = 0.01;
        steps = 500;
        ItemSize;
    end
    
    methods
        function obj = Sn(n)
            obj.ItemSize = n+1;
            obj.type = ['The ',num2str(n),'-sphere in R',num2str(n+1)];
        end
        function q = exp(this,p,v,t)
            % exp(p,v) - Exponential map at the point p with respect to v in
            % TpM
            %
            % INPUT
            %   p : a point or set of points on the manifold S2
            %   v : a point or set of point in the tangential spaces TpM
            %
            % OUTPUT
            %   q : resulting point(s) on S2
            % ---
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19
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
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19
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
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19 | 2015-03-30
            
            % Lofgile
            %   2015-03-30 Changed dist to work to collapse first dim
            %   2015-04-11 Extracted to Mex
            if this.useMex
                d = SnDist(p,q);
            else
                d = this.localDist(p,q);
            end
        end
        function m = midPoint(this,x,z)
            % m = midPoint(x,z)
            % Compute the (geodesic) mid point of x and z.
            %
            % INPUT
            %    x,z : two point( sets) as colum(s) in R3
            %
            % OUTPUT
            %      m : resulting mid point( sets)
            % ---
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19
            m = this.exp(x, this.log(x,z)./2);
        end
        function x = proxad(this,varargin)
            % proxad(f, lambda, w) - Compute the proximal mapping of the
            % data f with respect to lambda an the weight w.
            % Any unknown point containing a NaN is inpainted if
            % neighbouring known points or else left unchanged
            %
            % INPUT
            %      f : data on Sn (i.e. R(n+1)), where [n+1,l,d] = size(f) and d
            %          is the number of data points per prox and l is the
            %          number of parallel proxes. A point is unknown if for
            %          fixed d,l any of the three entries is NaN.
            % lambda : weight of the proximal map
            %      w : finite difference weight vector. For S2 up to now
            %          only [-1,1] and [-1,2,1] are supported
            %
            % OUTPUT
            %      x : resulting data of all proximal maps
            % ---
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19
            if (length(varargin)==1 && isstruct(varargin{1})) %struct case
                vars = varargin{1};
                assert(all(isfield(vars,{'f','lambda','w'})),...
                    'Not all required parameters given in the struct');
                if ~isfield(vars,'RegMask');
                    vars.RegMask = [];
                end
            else
                ip = inputParser;
                addRequired(ip,'f');
                addRequired(ip,'lambda');
                addRequired(ip,'w');
                addParameter(ip, 'RegMask', []);
                parse(ip, varargin{:});
                vars = ip.Results;
            end
            [k,l,d] = size(vars.f);
            assert(k==this.ItemSize,['The values of f (',num2str(d),' points in ',...
                num2str(l),' proximal mappings of dimension ',num2str(k),') are not lying in R3']);
            if (isrow(vars.w))
                w = vars.w';
            else
                w = vars.w;
            end
            if (vars.lambda==0)
                x = vars.f;
                return
            end
            %% TODO RegMask
            assert(d==length(w),['The length of the weight (',...
                num2str(length(w)),') does not fit to the data size (',...
                num2str(d),').']);
            x = vars.f;
            x(isnan(vars.f)) = NaN;
            missingpoints = squeeze(any(isnan(vars.f),1)); %sum along manifold dimensions & squeeze
            if l==1
                missingpoints = missingpoints';
            end
            % proxes that can be computed
            proxpts = sum(missingpoints,2)==0; %last dimension d
            % proxes that can be inpainted
            inpaintpts = sum(missingpoints,2)==1;
            if (length(w)==2) && all(w==[-1,1]')
                t = min(vars.lambda, 0.5*this.dist(vars.f(:,proxpts,1),vars.f(:,proxpts,2)));
                % Divide only those that are nonzero.
                dir1 = this.log(vars.f(:,proxpts,1),vars.f(:,proxpts,2));
                l1 = sqrt(sum(dir1.^2,1));
                dir1(:,l1>eps) = dir1(:,l1>eps)./repmat(l1(l1>eps),[this.ItemSize,1]);
                %
                dir2 = this.log(vars.f(:,proxpts,2),vars.f(:,proxpts,1));
                l2 = sqrt(sum(dir2.^2,1));
                dir2(:,l2>eps) = dir2(:,l2>eps)./repmat(l2(l2>eps),[this.ItemSize,1]);
                % permute brings one singleton dimension to the front in t
                
                x(:,proxpts,1) = this.exp(vars.f(:,proxpts,1), repmat(permute(t,[3,1,2]),[this.ItemSize,1]).*dir1);
                x(:,proxpts,2) = this.exp(vars.f(:,proxpts,2), repmat(permute(t,[3,1,2]),[this.ItemSize,1]).*dir2);
                %
                % Inpaint all first missing points, the second one is
                % existing due to limiting the missing number to 1
                x(:,(missingpoints(:,1)>0)&inpaintpts,1) = x(:,(missingpoints(:,1)>0)&inpaintpts,2);
                x(:,(missingpoints(:,2)>0)&inpaintpts,2) = x(:,(missingpoints(:,2)>0)&inpaintpts,1);
            elseif (length(w)==3) && all(w==[1,-2,1]')
                %
                % Iterative subgradient descent
                x = vars.f;
                G = zeros(size(x));
                % Gradient descent
                xopt = x(:,proxpts,:);
                xoptvals = this.dist(this.midPoint(x(:,proxpts,1),x(:,proxpts,3)),x(:,proxpts,2)); %x=f hence first term zero
                for gradcount = 1:this.steps
                    tauit = this.tau*vars.lambda/gradcount;
                    % GradX
                    G(:,proxpts,1) = this.log(x(:,proxpts,1), vars.f(:,proxpts,1)) + ...
                        vars.lambda*this.gradX(x(:,proxpts,1),x(:,proxpts,2),x(:,proxpts,3));
                    % Grad Y
                    V = this.log(x(:,proxpts,2),this.midPoint(x(:,proxpts,1),x(:,proxpts,3)));
                    normV = sqrt(sum(V.^2,1));
                    if any(normV>eps) %norm directions
                        V(:,normV>eps) = V(:,normV>eps)./repmat(normV(normV>eps),[this.ItemSize,1]);
                    end
                    G(:,proxpts,2) = this.log(x(:,proxpts,2),vars.f(:,proxpts,2)) + vars.lambda.*V;
                    % Grad Z
                    G(:,proxpts,3) = this.log(x(:,proxpts,3), vars.f(:,proxpts,3)) + ...
                        vars.lambda*this.gradX(x(:,proxpts,3),x(:,proxpts,2),x(:,proxpts,1));
                    % Gradient step
                    x = reshape(...
                        this.exp(reshape(x,this.ItemSize,[]), tauit*reshape(G,this.ItemSize,[])),...
                        [this.ItemSize,l,d]);
                    x(:,proxpts,:) = x(:,proxpts,:)./repmat(sqrt(sum(x(:,proxpts,:).^2,1)),[this.ItemSize,1,1]); %sec?
                    xoptt = x(:,proxpts,:);
                    % -(l<2) fixes this summation for one point
                    xvals = sum(this.dist(vars.f(:,proxpts,:),xoptt).^2,2) + this.dist(this.midPoint(xoptt(:,:,1),xoptt(:,:,3)),xoptt(:,:,2));
                    xvalsInd = xvals<xoptvals;
                    if any(xvalsInd(:))
                        xopt(:,xvalsInd,:) = xoptt(:,xvalsInd,:);
                        xoptvals(xvalsInd) = xvals(xvalsInd);
                    end
                end
                x(:,proxpts,:) = xopt;
                % No Inpainting for now
            elseif (length(w)==4) && all(w==[-1,1,-1,1]')
                % Iterative subgradient descent
                x = vars.f;
                % Gradient descent
                xopt = x(:,proxpts,:);
                xoptvals = this.dist(this.midPoint(x(:,proxpts,1),x(:,proxpts,3)),this.midPoint(x(:,proxpts,2),x(:,proxpts,4))); %x=f hence first term zero
                for gradcount = 1:this.steps
                    tauit = this.tau/gradcount;
                    %midpoints between pairs of points
                    M = zeros(size(x(:,proxpts,1:2)));
                    M(:,:,1) = this.midPoint(x(:,proxpts,1), x(:,proxpts,3));
                    M(:,:,2) = this.midPoint(x(:,proxpts,2), x(:,proxpts,4));
                    % Common
                    % Grad Xi
                    G = zeros(size(x(:,proxpts,:)));
                    for i=1:4
                        i2 = mod(i+1,4)+1; %other index involved in same mid point, 1-3, 2-4, 3-1, 4-2
                        G(:,:,i) = this.log(x(:,proxpts,i), vars.f(:,proxpts,i))...
                            + vars.lambda*this.gradX(x(:,proxpts,i),...
                            M(:,:,mod(i,2)+1),... %opposite mid point as Y
                            x(:,proxpts,i2));
                    end
                    x(:,proxpts,:) = reshape(...
                        this.exp(reshape(x(:,proxpts,:),this.ItemSize,[]), tauit*reshape(G,this.ItemSize,[])),...
                        [this.ItemSize,sum(proxpts),d]);
                    x(:,proxpts,:) = x(:,proxpts,:)./repmat(sqrt(sum(x(:,proxpts,:).^2,1)),[3,1,1]); %sec?
                    xoptt = x(:,proxpts,:);
                    xvals = sum(this.dist(vars.f(:,proxpts,:),xoptt).^2,2) + this.dist(this.midPoint(xoptt(:,:,1),xoptt(:,:,3)),this.midPoint(xoptt(:,:,2),xoptt(:,:,4)));
                    xopt(:,xvals<xoptvals,:) = xoptt(:,xvals<xoptvals,:);
                    xoptvals(xvals<xoptvals) = xvals(xvals<xoptvals);
                end %end of mixed derivative gradient steps
                x(:,proxpts,:) = xopt;
            else
                warning(['Unknown discrete difference on ',this.type,': ',num2str(w'),'. Returning the input f.']);
                x=vars.f;
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
        function fn = OldaddNoise(this,f,sigma)
            if this.ItemSize==3 %S2
                fs = size(f);
                f = reshape(f,3,[]);
                [AZ,EL,R] = cart2sph(f(1,:),f(2,:),f(3,:));
                angles = sigma*randn([2,prod(fs(2:end))]);
                AZ = AZ+angles(1,:);
                EL = EL+angles(2,:);
                [fx,fy,fz] = sph2cart(AZ,EL,R);
                fn = zeros(size(f));
                fn(1,:) = fx; fn(2,:) = fy; fn(3,:) = fz;
                fn = reshape(fn,fs);
            else
                %Good Idea?
                fn = f + nthroot(sigma,this.ItemSize)*ones(size(f));
                fn = fn./repmat(sqrt(sum(fn.^2,1)),[this.ItemSize,ones(1,length(size(f)))]);
            end
        end
    end
    methods (Access = private)
        function G = gradX(this,X,Y,Z)
%            if this.useMex && (this.ItemSize==3)
%                G = SnGradX(X,Y,Z);
%            else
                l = size(X,2);
                m = this.exp(X, this.log(X,Z)./2);
                % Common directions to middle points y
                r = this.log(m,Y);
                normr = permute(sqrt(sum(r.^2,1)),[2:length(size(r)),1]);
                % Grad X
                V = this.TpMONB(m,Z);
                W = this.TpMONB(X,Z);
                alpha = zeros(l,this.ItemSize-1);
                if any(normr>eps)
                    alpha(normr>eps,:) = sum(...
                        repmat(r(:,normr>eps),[1,1,this.ItemSize-1]).*...
                        V(:,normr>eps,:)./repmat(permute(normr(normr>eps),[2,1]),[this.ItemSize,1,this.ItemSize-1]),1);
                end
                G = W(:,:,1).*repmat(0.5*permute(alpha(:,1),[2,1]),[this.ItemSize,1]) ...
                    + W(:,:,2:this.ItemSize-1).*...
                    repmat(...
                    repmat(1./(2*cos(permute(this.dist(X,Z)./2,[3,1,2]))),[1,this.ItemSize-2]).*...
                    permute(alpha(:,2:this.ItemSize-1),[2,1]),...
                    [this.ItemSize,1,1]);
 %           end
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
                q(:,normv>eps) = p_(:,normv>eps).*repmat(cos(t*normv(normv>eps)),[this.ItemSize,1])...
                    + v_(:,normv>eps).*repmat(sin(t*normv(normv>eps))./normv(normv>eps),[this.ItemSize,1]);
            end
        end
        function v = localLog(this,p_,q_)
            % Fallback is no mex available
            if any(all(size(p_)~=size(q_)))
                error('p and q have to be of same size');
            end
            % Compute all scalar products of points p and q, secure, that
            % they are in -1,1 for numerical reasons
            % scp = shiftdim(min( max(bsxfun(@dot,p,q),-1), 1));
            %
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
        end
        function V = TpMONB(this,p,q)
            % V = TpMONB(p,q)
            % Compute an ONB in TpM, where the first vector points to q.
            %
            % INPUT
            %     p : base point( sets)
            %     q : directional indicator( sets) for the first vector(s).
            %
            % OUTPUT
            %    V : basiscolumn matrice(s)
            % ---
            % ManImRes 1.0, R. Bergmann ~ 2014-10-19, last edit: 2014-10-23
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
            V = zeros(this.ItemSize,size(p_,2),this.ItemSize-1);
            V(:,:,1) = this.log(p_,q_);
            normsv = sqrt(sum(V(:,:,1).^2,1));
            V(:,normsv>eps,1) = V(:,normsv>eps,1)./repmat(normsv(normsv>eps),[this.ItemSize,1,1]);
            if this.ItemSize==3 %S2 -> cross
                V(:,:,2) = cross( squeeze(V(:,:,1)), p_);
            else
                for col=1:size(V,2)
                    if (normsv(col)>eps)
                        % The remaining Tangential vectors are the orthogonal
                        % to p(:,col) and V(:,col,1), i.e. the nullspace of the
                        % matrix p V V ... V
                        V(:,col,2:this.ItemSize-1) = null([p_(:,col), repmat(V(:,col,1),[1,this.ItemSize-1])]);
                    end
                end
            end
        end
    end
end