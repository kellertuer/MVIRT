classdef S1mRn < manifold & handle
    % The product manifold of combined cyclic and real valued data, i.e.
    % (S^1)^n x R^m, so its an k=n+m dimensional space
    % ---
    % ---
    % Manifold-valued Image Restoration Toolbox 1.1
    % R. Bergmann ~ 2014-10-22 | 2016-09-15
    % see LICENSE.txt
    properties
        % Binary vectors indicating S and R components of lenthg of the
        % whole space dimension
        type='S1mRn';
        SComponents;
        RComponents;
        ItemSize;
        Dimension;
    end
    
    methods
        function obj = S1mRn(m,n)
        % SnRm(n,m) Create a product manifold, where
        % - the first n components are cyclic
        % - the last m components are real valued
            if (length(m) == length(n)) && all(islogical(m)) && all(islogical(n))
                obj.SComponents = m;
                obj.RComponents = n;
            else % numbers
                k = n+m;
                obj.SComponents = (1:k)<=m;
                obj.RComponents = (1:k)>m;
                obj.type = ['S',num2str(m),'R',num2str(n)];
            end
            obj.ItemSize = length(obj.SComponents);
            obj.Dimension = obj.ItemSize;
        end
        function q = exp(this,p,v,t)
            % exp(p,v) - Exponential map at the point p with respect to v in
            % TpM on S1.
            %
            % INPUT
            %   p : a point or set (columns) of points on the manifold SnRm
            %   v : a point or set (columns) of point in the tangential spaces TpM
            % OPTIONAL
            %   t : [] following the geodesics for one time step
            %       a scalar following the geodesics for time t
            %       a set of t following each geodesic_i for time t_i
            % OUTPUT
            %   q : resulting point(s) on SnRm
            % ---
            % ManImRes 1.0, R. Bergmann ~ 2014-10-26
            dimen = size(p);
            if nargin < 4
                t = 1;
            end
            if ~isscalar(t)
                assert(all(dimen(2:end)==size(t)),'t has to match data dimeension');
            end
            assert(all(size(p)==size(v)),'p and q have to be of same size');
            assert( (size(p,1)==this.ItemSize) && (size(v,1)==this.ItemSize), ...
                'The manifold dimension of either p or v is wrong');
            p_ = reshape(p,this.ItemSize,[]);
            if isscalar(t)
            v_ = t*reshape(v,this.ItemSize,[]);
            else
                v_ = reshape(repmat(permute(t,[length(dimen),1:length(dimen)-1]),[this.ItemSize,ones(1,length(dimen)-1)]).*v,this.ItemSize,[]);
            end
            q = v_+p_;
            q(this.SComponents,:) = symMod(q(this.SComponents,:),2*pi);
            q = reshape(q,size(p));
        end
        function v = log(this,p,q)
            % log(q,p) - Inverse Exponential Map at p of q.
            %
            % INPUT
            %    p : point or set of (column) points indicating the
            %        tangential base points
            %    q : point(s) on SmRn (i.e. [-pi,pi)^mR^ being put into the
            %        tangential plane at their corresponding p
            %
            % OUTPUT
            %    v : points on the tangential plane at point(s) p
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-22
            assert(all(size(p)==size(q)),'p and q have to be of same size');
            v = q-p;
            v(this.SComponents,:) = symMod(v(this.SComponents,:),2*pi);
        end
        function x = proxDist(this,g,f,lambda)
            % proxDist(g,f,lambda)
            % Proximal step towards f from given data g with parameter
            % lambda on (S1)mRn arbitrary manifold.
            % INPUT
            %  f,g    : data point( sets/columns )
            %  lambda : stepsize towards f
            %
            % OUTPUT
            %       x : result point( sets) of the proximal map
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            t = this.dist(g,f,0); %do not sum
            v = this.log(g,f);
            % t = repmat(t,[this.ItemSize,1]); %repeat t into the dimensions
            % Instead of the line above - for product manifolds we do the
            % following
            v = sign(v); %norm v in each factor of the product manifold
            x = this.exp(g, repmat(lambda./(1+lambda),[this.ItemSize,1]).*t.*v);
        end
        function d = dist(this,p,q,n)
            % dist(a,b) computes the length of the smaller arc of a,b on
            % SnRm and returns the vector of distances. If n is specified,
            % the n-norm of these distances is computed
            %    INPUT
            %        p,q    : 2 point sets (columns) on the SnRm data
            %        n      : (optional on product manifolds) norm to combine the
            %                 distances with n-norm (set to 0 to
            %                 deactivate)
            %
            %    OUTPUT
            %        d      : lengths of the shorter arcs between on S components and abs on R, then the norm
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ created 2013-10-25 ~ last edited 2014-10-22
            if length(this.SComponents)>1 %number of components more than one.
                if isrow(p) % one datum
                    p_ = p';
                else
                    p_ = p;
                end
                if isrow(q)
                    q_ = q';
                else
                    q_ = q;
                end
            else
                p_=p; q_=q;
            end
            assert( all(size(p_)==size(q_)), ...
                'Distances can only be computed for equal length vectors p and q');
            assert( (size(p,1)==this.ItemSize) && (size(q,1)==this.ItemSize), ...
                'The manifold dimension of either p or q is wrong');
            d = reshape(p_-q_,this.ItemSize,[]);
            d(this.SComponents,:) = symMod(d(this.SComponents,:) ,2*pi);
            d = abs(d);
            d = reshape(d,size(p_));
            if (nargin < 4)
                n = 2;
            end
            if n~=0
                d = shiftdim(sum( d.^n, 1).^(1/n)); %eliminate leading zeros
            end
        end
        function x = proxad(this,varargin)
        % data f with respect to lambda an the weight w.
        % Any unknown point containing a NaN is inpainted if
        % neighbouring known points or else left unchanged
        %
        % INPUT
        %      f : data on SnRm, k=n+m, where [k,l,d] = size(f) and d
        %          is the number of data points per prox and l is the
        %          number of parallel proxes. A point is unknown if for
        %          fixed d,l any of the three entries is NaN.
        % lambda : weight of the proximal map
        %      w : finite difference weight vector. For S2 up to now
        %          only [-1,1] and [-1,2,1] are supported
        %
        % OUTPUT
        %      x : resulting data of all proximal maps
        %
        % OPTIONAL PARAMETERS
        %     'RegMask'     : ([]) Specify a binary mask for values affected by
        %                      (1) they are regularized, i.e. moved (or active, say)
        %                      (0) fixed. RegMask is of the same size
        %                      the data point set, i.e. [d,l].
        % ---
        % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-26 | 2014-02-13
        
        % Changelog
        % 2015-02-13 Changed order of data, such that the manifold values
        %            are first
        
        if (length(varargin)==1 && isstruct(varargin{1})) %struct case
            vars = varargin{1};
            assert(all(isfield(vars,{'f','lambda','w'})),...
                'Not all required parameters given in the struct');
            [k,l,d] = size(vars.f);
            if ~isfield(vars,'RegMask');
                 vars.RegMask = ones(l,d); %move all
            end
        else
            ip = inputParser;
            addRequired(ip,'f');
            addRequired(ip,'lambda');
            addRequired(ip,'w');
            addParameter(ip, 'RegMask', []);
            parse(ip, varargin{:});
            vars = ip.Results;
            [k,l,d] = size(vars.f);
            if (numel(vars.RegMask)>0) % Mask given
                if (any(size(vars.RegMask)~=[l,d]))
                    warning(['Size of the regularization mask (',...
                        num2str(size(vars.regMask)),') does not fit the size of f (',...
                        num2str([l,d]),'). Hence it is ignored.']);
                    vars.RegMask = ones(l,d); %move all
                end
            else
                vars.RegMask = ones(l,d); %move all
            end
            assert(k==length(this.SComponents),['The values of f (',...
                num2str(d),' points in ',num2str(l),...
                ' proximal mappings of dimension ',num2str(k),...
                ') are not of dimension ',num2str(length(this.SComponents)),'.']);
            if (isrow(vars.w))
                vars.w = vars.w';
            else
                vars.w = vars.w;
            end
            if (vars.lambda==0)
                x = vars.f;
                return
            end
            assert(d==length(vars.w),['The length of the weight (',...
                num2str(length(vars.w)),') does not fit to the data size (',...
                num2str(d),').']);
        end
            x = vars.f;
            x(isnan(vars.f)) = NaN;
            missingpoints = any(isnan(vars.f),1);
            % proxes that can be computed
            proxsets = (sum(missingpoints,3)==0)';
            % proxes that can be inpainted
            inpaintpts = (sum(missingpoints,3)==1)'; %count all, where just one is missing
            missingpoints = squeeze(missingpoints);
            if sum(inpaintpts)>0 && length(vars.w)==2 && all(vars.w==[-1;1]) %data exists, that can be inpainted with first order - do so.
                % first point missing - set to second
                x(:, (missingpoints(:,1)>0)&inpaintpts,1) = x(:, (missingpoints(:,1)>0)&inpaintpts,2);
                % second point missing - set to first
                x(:, (missingpoints(:,2)>0)&inpaintpts,2) = x(:, (missingpoints(:,2)>0)&inpaintpts,1);
            end
            if sum(proxsets)>0 %any points to prox?
                % scalar products <f,w> for each row/column in 1&3 dim ->
                % lxk result and also take S-components into account
                sp =  dot(repmat(permute(vars.w,[3,2,1]),[k,sum(proxsets),1]),vars.f(:,proxsets,:),3);
                % sum over last dimension -> kxlxd result
                sp(this.SComponents,:,:) = symMod(sp(this.SComponents,:,:),2*pi);
                spn = sqrt(sum(sp.^2,1)); %sp norm -> 1xl
                s = sp./repmat(spn,[k,1,1]); %norm -> s - kxl
                s(isnan(s)) = 0; %else case
                % ||w^a||_2^2
                nwsq = sum(permute(vars.RegMask(proxsets,:),[2,1]).*repmat(vars.w,[1,sum(proxsets)]).^2,1); % 1xl
                % m -> 1xlx1 (one for each prox
                mins = min( ones(size(spn)).*vars.lambda,spn./nwsq); % 1xl
                x(:,proxsets,:) = vars.f(:,proxsets,:) - ...
                    repmat(permute(...%into all dimensions
                        permute(vars.RegMask(proxsets,:),[2,1]).*repmat(vars.w,[1,sum(proxsets)]),...
                        [3,2,1]),[k,1,1]).*...
                    repmat(s.*repmat(mins,[k,1,1]),[1,1,d]);
            end
            % modulo all computed on S components
            x(this.SComponents,inpaintpts|proxsets,:) = ...
                symMod(x(this.SComponents,inpaintpts|proxsets,:),2*pi);
        end
        function fn = addNoise(this,f,sigma)
            fn = f + sigma*randn(size(f));
            fn = reshape(fn,this.ItemSize,[]);
            fn(this.SComponents,:) = symMod(fn(this.SComponents,:),2*pi);
            fn = reshape(fn,size(f));
        end
        function ds = dot(this,P,V,W)
            % Sm1Rn.dot(P,V,W)
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
            % MVIRT 1.1, J. Persch 2016-06-13
            %
            dimen = size(P);
            if all(size(V) == dimen & size(W) == dimen)
                ds = permute(sum(V.*W,1),[2:length(dimen),1]);
            else
                error('Dimensions of Input must coincide')
            end
        end
        function V = TpMONB(this,p)
            % V = TpMONB(p,q)
            % Compute an ONB in TpM, where the first vector points to q,
            % whin given.
            %
            % INPUT
            %     p : base point( sets)
            % OPTIONAL:
            %     q : directional indicator( sets) for the first vector(s).
            %
            % OUTPUT
            %    V : basiscolumn matrice(s)
            % ---
            % MVIRT 1.1, J. Persch 2016-06-24
            if isrow(p)
                p_=p';
            else
                p_=p;
            end
            dimen = size(p_);
            V = permute(eye(this.Dimension),[1,3:length(dimen)+1,2]);
            V = repmat(V,[1,dimen(2:end),1]);
        end
    end
end
