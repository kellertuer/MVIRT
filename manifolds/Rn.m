classdef Rn < manifold & handle
    % The manifold of the usual n-dimensional euclidean space.
    % ---
    % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~2014-10-22
    properties
        type='Rn';
        ItemSize;
        Dimension;        
        allDims;
    end
    
    methods
        function obj = Rn(n)
            obj.type = ['R',num2str(n)];
            obj.ItemSize = n;
            obj.Dimension = n;
            obj.allDims = repelem({':'},length(obj.ItemSize));
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
                if isvector(t)
                    s=[this.ItemSize,length(t)];
                else
                    s = [this.ItemSize,size(t)];
                end
                assert(all(s == size(v)),'t should be scalar or match dimension of v');
                q = p+repmat(permute(t,[length(size(t))+1,1:length(size(t))]),[this.ItemSize,ones(size(size(t)))]).*v;
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
                d = shiftdim(sum( d.^norm, 1).^(1/norm),1); %eliminate leading zeros
            end
        end
        function W = parallelTransport(~,~,~,V)
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
        function [V,k] = TpMONB(this,p,~)
           dimen = size(p);
           V = permute(repmat(eye(this.ItemSize),[1,1,dimen(2:end)]),[1,3:length(dimen)+1,2]);
           if nargout > 1
                k = zeros([dimen(2:end),this.Dimension]);
            end 
        end
        function G = grad_X_D2_Sq(~,X,Y,Z)
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
