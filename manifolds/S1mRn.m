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
        allDims;
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
            obj.allDims = repelem({':'},length(obj.ItemSize));
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
        function fn = addNoise(this,f,sigma)
            fn = f + sigma*randn(size(f));
            fn = reshape(fn,this.ItemSize,[]);
            fn(this.SComponents,:) = symMod(fn(this.SComponents,:),2*pi);
            fn = reshape(fn,size(f));
        end
        function ds = dot(~,P,V,W)
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
