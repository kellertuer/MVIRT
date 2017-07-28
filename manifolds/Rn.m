classdef Rn < manifold & handle
    % The manifold of the usual n-dimensional euclidean space.
    % ---
    % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~2014-10-22
    properties
        type='Rn';
        ItemSize;
        Dimension;
    end
    
    methods
        function obj = Rn(n)
            obj.type = ['R',num2str(n)];
            obj.ItemSize = n;
            obj.Dimension = n;
        end
        function q = exp(this,p,v,t)
            % exp(p,v) - Exponential map at the point p with respect to v in
            % TpM.
            %
            % INPUT
            %   p : a point or set (columns) of points on the manifold Rn
            %   v : a point or set (columns) of point in the tangential spaces TpM
            %
            % OUTPUT
            %   q : resulting point(s) on Rn
            % ---
            %  % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-26
            %               J. Persch ~ 2016-10-26
            if nargin < 4
                t=1;
            end
            if isscalar(t)
                q = p+t*v;
            else
                assert([this.ItemSize,size(t)] == size(v),'t should be scalar or match dimension of v');
                q = p+repmat(permute(t,[length(size(t))+1,1:length(size(t))]),[this.ItemSize,ones(size(size(t)))])*v;
            end
        end
        function v = log(~,p,q)
            % log(q,p) - Inverse Exponential Map at p of q.
            %
            % INPUT
            %    p : point or set of (column) points indicating the
            %        tangential base points (in Rn)
            %    q : point(s) on Rn being put into the
            %        tangential plane at their corresponding p
            %
            % OUTPUT
            %    v : points on the tangential plane at point(s) p
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-22
            assert(all(size(p)==size(q)),'p and q have to be of same size');
            v = q-p;
        end
        function x = proxDist(this,g,f,lambda)
            % proxDist(g,f,lambda)
            % Proximal step towards f from given data g with parameter
            % lambda on an arbitrary manifold.
            % INPUT
            %  f,g    : data point( sets/columns )
            %  lambda : stepsize towards f
            %
            % OUTPUT
            %       x : result point( sets) of the proximal map
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            t = this.dist(g,f,0);
            v = this.log(g,f);
            v = sign(v); %norm v in each factor of the product manifold
            x = this.exp(g, lambda/(1+lambda).*t.*v);
        end
        function d = dist(this,p,q,norm)
            % dist(a,b) computes the length of the smaller arc of a,b on
            % SnRm and returns the vector of distances. If n is specified,
            % the n-norm of these distances is computed
            %    INPUT
            %        p,q    : 2 point sets (columns) on the SnRm data
            %    
            %    OPTIONAL
            %       norm     : (2) norm to use (0 to disable, i.e. stay elementwise)
            %
            %    OUTPUT
            %        d      : lengths of the shorter arcs between on S components and abs on R, then the norm
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ created 2013-10-25 ~ last edited 2014-10-22
            assert(all(size(p)==size(q)), ...
                'Distances can only be computed for equal length vectors p and q');
            assert((size(p,1)==this.ItemSize) && (size(q,1)==this.ItemSize), ...
                'The manifold dimension of either p or q is wrong');
            d = abs(p-q);
            if (nargin < 4)
                norm = 2;
            end
            if norm~=0
                d = shiftdim(sum( d.^norm, 1).^(1/norm)); %eliminate leading zeros
            end
        end
        function W = parallelTransport(this,X,Y,V)
            % W = parallelTransport(X,Y,V) parallel transport a tangential
            % vector V at X along a geodesic from X to Y to Y
            %
            % INPUT
            %    X : a point( set) on P(m)
            %    Y : a point( set) on P(m)
            %    V : a tangential vector( set, one) at (each) X
            %
            % OPTIONAL
            %    t : value or vector to only take a fraction along the
            %    geodesic(s)
            %
            % OUTPUT
            %    W : tangential vector( set, one) at (each) Y
            %
            % ---
            % ManImRes 1.0, R. Bergmann ~ 2015-01-29 | 2015-04-10
            
            % Changelog
            %   2015-04-10 Introduced Mex-Files
            if nargin > 4
                error('Too many input arguments for parallelTransport');
            elseif nargin< 4
                error('Not enough input arguments for parallelTransport');
            end
            W = V;
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
            assert(k==length(this.ItemSize),['The values of f (',...
                num2str(d),' points in ',num2str(l),...
                ' proximal mappings of dimension ',num2str(k),...
                ') are not of dimension ',num2str(length(this.ItemSize)),'.']);
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
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch, 2016-10-21
            G = -1/2*(X+Z-2*Y);
        end
        function ds = dot(this,P,V,W)
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
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch 2016-10-23
            %
            dimen = size(P);
            if all(size(V) == dimen & size(W) == dimen)
                ds = permute(sum(V.*W,1),[2:length(dimen),1]);
            else
                error('Dimensions of Input must coincide')
            end
        end
        function fn = addNoise(~,f,sigma)
            fn = f + sigma*randn(size(f));
        end
    end
end
