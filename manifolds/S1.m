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
        function y = exp(~,x,xi,t)
            % exp(x,xi) - Exponential map at x with respect to xi in TxS1
            %
            % INPUT
            %   x : a point or set of points on the manifold S1
            %   xi : a point or set of point in the tangential spaces TpM
            %
            % OUTPUT
            %   y : resulting point(s) on S1
            %
            % OPTIONAL
            %   t : shorten xi by a factor.
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19
            if nargin<4
                t=1;
            end
            if ~isscalar(t)
                size_p = size(x);
                assert(all(size_p(2:end)==size(t)),'t should either be a scalar or have the same size as p'); 
            end
            assert(all(size(x)==size(xi)),'p and q have to be of same size');
            y = symMod(t*xi+x,2*pi);
        end
        function xi = log(~,x,y)
        % log(x,y) - Inverse Exponential Map at x to y.
        %
        % INPUT
        %    x : point or set of (column) points indicating the
        %        tangential base points
        %    y : point(s) on S1 (i.e. [-pi,pi) being put into the
        %        tangential plane at their corresponding p
        %
        % OUTPUT
        %    xi : points on the tangential plane at point(s) p
        % ---
        % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-22
            assert(all(size(x)==size(y)),'p and q have to be of same size');
            xi = symMod(y-x,2*pi);
        end
        function d = dist(~,x,y)
        % dist(a,b) computes the length of the smaller arc of a,b on S1
        %   also works for arbitrary sized array a and b componentwise
        %    INPUT
        %        x,y    : 2 point sets on the S1 = [-pi,pi)
        %
        %    OUTPUT
        %        d      : lengths of the shorter arcs between a and b
        % ---
        % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2013-10-25 | 2014-10-22
            if iscolumn(x)
                x_ = x';
            else
                x_ = x;
            end
            if iscolumn(y)
                y_ = y';
            else
                y_ = y;
            end
            if (any(size(x_)~=size(y_)))
                error('Distances can only be computed for equal length vectors p and q');
            end
            d = shiftdim(abs(symMod(x_-y_,2*pi)),1); %squeeze first dimension if both are not vectors
        end
        function [Xi,k] = TpMONB(this,p,~)
        % [V,k] = TpMONB(p,q) an onb of the tangent space including log_pq
        %   since on S1 there is only one direction, this is always
        %   included and 1; the curvature coefficients are all zero, since
        %   this is always along the geodesic
        %
        % INPUT
        %   p,q : two point(sets) of cyclic data
        %
        % OUTPUT
        %    Xi : basiscolumn matrice(s).
        %    k : (optional) curvature coefficients, here all are 0.
        % ---
        % MVIRT 1.0, R. Bergmann ~ 2018-03-03
           dimen = size(p);
           Xi = permute(repmat(eye(this.ItemSize),[1,1,dimen(2:end)]),[1,3:length(dimen)+1,2]);
           if nargout > 1
                k = zeros([dimen(2:end),this.Dimension]);
           end
        end
        function fn = addNoise(~,f,sigma)
            % fn = addNoise(f,sigma) add wrapped Gaussian noise.
            %
            % INPUT
            %   f   : cyclic data
            % sigma : standard deviation
            %
            % OUTPUT
            %  fn   :  noisy data (f+sigma*randn(size(f))_{2pi}
            % ---
            % Manifold-valued Image Restoration Toolbox
            % R. Bergmann ~ 2018-03-03
            fn = symMod(f + sigma*randn(size(f)),2*pi);
        end
        function eta = parallelTransport(~,~,~,xi)
            % W = parallelTransport(x,y,xi) parallel transport a tangential
            % vector xi at x along a geodesic from x to y
            %
            % INPUT
            %    x : a point( set) on P(m)
            %    y : a point( set) on P(m)
            %    xi : a tangential vector( set, one) at (each) X
            %
            %
            % OUTPUT
            %    eta : tangential vector( set, one) at (each) Y
            %
            % ---
            % Manifold-valued Image Restoration Toolbox
            % R. Bergmann ~ 2015-01-29 | 2015-04-10
            
            % Changelog
            %   2015-04-10 Introduced Mex-Files
            if nargin > 4
                error('Too many input arguments for parallelTransport');
            elseif nargin< 4
                error('Not enough input arguments for parallelTransport');
            end
            eta = xi;
        end
        function ds = dot(~,x,xi,nu)
            % S1.dot(x,xi,nu)
            %     Compute the inner product of two tangent vectors in T_P M
            %
            % INPUT
            %     x  : a point(Set) in S1
            %     xi  : a first tangent vector( set) to (each) x
            %     nu  : a secod tangent vector( set) to (each) x
            %
            % OUTPUT
            %     ds : the corresponding value(s) of the inner product of (each triple) V,W at X
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch 2016-12-06
            %
            dimen = size(x);
            if all(size(xi) == dimen & size(nu) == dimen)
                ds = permute(xi.*nu,[2:length(dimen),1]);
            else
                error('Dimensions of Input must coincide')
            end
        end
    end
end
