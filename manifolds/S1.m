classdef S1 < manifold & handle
    % The manifold of the circle or S1
    %
    % ---
    % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~2014-10-22
    properties
        type = 'S1';
        ItemSize = 1;
        Dimension =1;        
        allDims;
    end
    
    methods
        function obj = S1()
            obj.allDims = repelem({':'},length(obj.ItemSize));
        end
        function q = exp(~,p,v,t)
            % exp(p,v) - Exponential map at the point p with respect to v in
            % TpM on S1.
            %
            % INPUT
            %   p : a point or set of points on the manifold S1
            %   v : a point or set of point in the tangential spaces TpM
            %
            % OUTPUT
            %   q : resulting point(s) on S1
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            if nargin<4
                t=1;
            end
            if ~isscalar(t)
                size_p = size(p);
                assert(all(size_p(2:end)==size(t)),'t should either be a scalar or have the same size as p'); 
            end
            assert(all(size(p)==size(v)),'p and q have to be of same size');
            q = symMod(t*v+p,2*pi);
        end
        function v = log(~,p,q)
        % log(q,p) - Inverse Exponential Map at p of q.
        %
        % INPUT
        %    p : point or set of (column) points indicating the
        %        tangential base points
        %    q : point(s) on S1 (i.e. [-pi,pi) being put into the
        %        tangential plane at their corresponding p
        %
        % OUTPUT
        %    v : points on the tangential plane at point(s) p
        % ---
        % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-22
            assert(all(size(p)==size(q)),'p and q have to be of same size');
            v = symMod(q-p,2*pi);
        end
        function d = dist(~,p,q)
        % dist(a,b) computes the length of the smaller arc of a,b on S1
        %   also works for arbitrary sized array a and b componentwise
        %    INPUT
        %        p,q    : 2 point sets on the S1 = [-pi,pi)
        %
        %    OUTPUT
        %        d      : lengths of the shorter arcs between a and b
        % ---
        % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2013-10-25 | 2014-10-22
            if iscolumn(p)
                p_ = p';
            else
                p_ = p;
            end
            if iscolumn(q)
                q_ = q';
            else
                q_ = q;
            end
            if (any(size(p_)~=size(q_)))
                error('Distances can only be computed for equal length vectors p and q');
            end
            d = shiftdim(abs(symMod(p_-q_,2*pi)),1); %squeeze first dimension if both are not vectors
        end
        function [V,k] = TpMONB(this,p,~)
           dimen = size(p);
           V = permute(repmat(eye(this.ItemSize),[1,1,dimen(2:end)]),[1,3:length(dimen)+1,2]);
           if nargout > 1
                k = zeros([dimen(2:end),this.Dimension]);
           end
        end
        function fn = addNoise(~,f,sigma)
            fn = symMod(f + sigma*randn(size(f)),2*pi);
        end
        function ds = dot(~,P,V,W)
            % S1.dot(P,V,W)
            %     Compute the inner product of two tangent vectors in T_P M
            %
            % INPUT
            %     X  : a point(Set) in S1
            %     V  : a first tangent vector( set) to (each) X
            %     W  : a secod tangent vector( set) to (each) X
            %
            % OUTPUT
            %     ds : the corresponding value(s) of the inner product of (each triple) V,W at X
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch 2016-12-06
            %
            dimen = size(P);
            if all(size(V) == dimen & size(W) == dimen)
                ds = permute(V.*W,[2:length(dimen),1]);
            else
                error('Dimensions of Input must coincide')
            end
        end
    end
end
