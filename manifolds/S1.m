classdef S1 < manifold & handle
    % The manifold of the circle or S1
    %
    % ---
    % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~2014-10-22
    properties
        type = 'S1';
        ItemSize = 1;
        Dimension =1;        
    end
    
    methods
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
        function x = proxad(this,varargin)
            % proxad(f,lambda,w) compute the [PROX]imal mapping for
            %       [A]bsolute ([C]yclic) [D]ifferences for each column of f with
            %       respect to weight w, where the length of w has to be the same as
            %       the column length (rows) of f. Unknown (NaN) entries of f are
            %       computed if possible, i.e. if they are the only unknown value in
            %       one column.
            % proxad(problem)
            %       Call the same function but with a struct containing both
            %       the mandatory (f,lambda,w) and optional parameters (see
            %       below)
            %
            %
            % INPUT
            %    f      : each column of f is threatened as one data set
            %             if f is a row vector, it isjust threated as one data set
            %    lambda : usual value of the proximum
            %    w      : weights corresponding to each column of f and with sum 0
            %
            %   OUTPUT
            %    x      : proximal mappings
            %
            % OPTIONAL PARAMETERS
            % 'RegMask'     : ([]) Specify a binary mask for values affected by
            %                  (1) they are regularized, i.e. moved (or active, say)
            %                  (0) fixed. RegMask is of the same size
            %                  the data point set, i.e. the number of rows.
            %
            %
            % See also:
            %    cppa_ad_1D, cppa_ad_2D, cppa_ad_nD
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0 ~ R. Bergmann ~ 2014-03-23 | 2014-11-29
            if length(varargin)==1 && isstruct(varargin{1}) %struct case - less checks
                vars = varargin{1};
                assert(all(isfield(vars,{'f','lambda','w'})),...
                    'Not all required parameters given in the struct');
                f = vars.f;
                [~,l,d] = size(f); %number of sets (l) and length of each data set (d)
                if ~isfield(vars,'RegMask');
                  regMask = ones(l,d); %move all
                else
                  regMask = vars.RegMask;
                end
                w = vars.w;
            else
                ip = inputParser;
                addRequired(ip,'f');
                addRequired(ip,'lambda');
                addRequired(ip,'w');
                addParameter(ip, 'RegMask', []);
                parse(ip, varargin{:});
                vars = ip.Results;
                if isrow(vars.w) %if w is a single row -> use as column.
                    vars.w = vars.w';
                end
                if sum(vars.w)~=0
                    error('The sum of the weights w has to be zero.');
                end
                w = vars.w;
                f = vars.f;
                if (size(f,2)==1) && (length(size(vars.f))==2) %one column data
                    f = f';
                end
                if (isrow(f)) % 1xd - ein Datensatz
                    f = permute(f,[3,1,2]); %1x1xd
                elseif (length(size(f))==2) %lxd
                    f = permute(f,[3,1,2]); % to 1xlxd
                end
                [k,l,d] = size(f); %number of sets (l) and length of each data set (d)
                if k~=this.ItemSize
                    error('No S1 data!');
                end
                if length(w)~=d
                    error(['Wrong amount of weights w (length: ',...
                    num2str(length(w)),') for data f (length: ',num2str(d),').']);
                end
                regMask = vars.RegMask;
                if numel(regMask)>0 % Mask given
                    if any(size(regMask)~=[size(f,2),size(f,3)])
                        warning(['Size of the regularization mask (',...
                            num2str(size(regMask)),') does not fit the size of f (',...
                            num2str([l,d]),'). Hence it is ignored.']);
                        regMask = ones(l,d); %move all
                    end
                else
                    regMask = ones(l,d); %move all
                end
            end
            %% Divide f into 3 sets:
            % those, where it is possible to inpaint
            inpaintcols = permute(sum(isnan(f(1,:,:)),3)==1,[2,1]);
            l1 = sum(inpaintcols);
            % those that do not contain NaNs
            proxcols = permute(sum( isnan(f(1,:,:)),3)==0,[2,1]);
            l2 = sum(proxcols);
            % and the rest is not handled at all
            x = zeros(1,l,d);
            % third: Those with more than one NaN -> keep them unchanged
            if (l1+l2)<l
                x(:,~(inpaintcols|proxcols),:) = f(:,~(inpaintcols|proxcols),:);
            end
            %
            %% Inpaint the first case
            if (l1>0) % there are possible values to pe intainted - not yet debugged
                    W = permute(repmat(w,[1,l1]),[2,1]);
                    prex = permute(f(1,inpaintcols,:),[2,3,1]);
                    %missing column entries, i.e. rows of '
                    [missingW,~] = find(isnan(prex'));
                    divW = diag(W(1:l1,missingW));
                    predot = prex.*W;
                    predot(isnan(predot))=0;
                    prex = prex.'; %work on rows
                    prex(isnan(prex)) = symMod(-sum(predot,2)./(divW),2*pi);
                    prex = prex.'; %after filling, get back
                    x(1,inpaintcols,:) = prex;
            end
            %% prox the second case
            % scalar product <f,w> for each row
            if (l2>0) % there are proxes...
                sp = symMod( dot(repmat(w,[1,l2]),permute(f(:,proxcols,:),[3,2,1])), 2*pi);
                % By Theorem 3.5 of BSW14 - ommitting the + case for
                s = sign(sp);
                % elementwise minima, denoted m in Thm 3.5
                %
                % If lambda is a number - it is applied to all elements, if it is a
                % vector, it has to mach the dimensions (1,m)
                nwsq = sum(permute(regMask(proxcols,:),[2,1]).*repmat(w,[1,l2]).^2,1);
                mins = min(ones(1,l2).*vars.lambda, abs(sp)./nwsq);
                % Multiply elementwise the sign (s) - the rest works the same as
                % usual mult
                x(:,proxcols,:) = symMod(f(:,proxcols,:) - permute(regMask(proxcols,:),[3,1,2]).*permute(w*(s.*mins),[3,2,1]),2*pi); % f + s*m*w for each row of f
            end
            % Rearrange if input was squeezed somehow
            if (size(vars.f,2)==1) && (length(size(vars.f))==2) %one column data
                x = permute(x,[3,1,2]);
            end
            if (isrow(vars.f)) || (length(size(vars.f))==2) %was mxd % 1xd - ein Datensatz
                x = permute(x,[2,3,1]); % 1xmxd -> mxdx1
            end
        end
        function fn = addNoise(~,f,sigma)
            fn = symMod(f + sigma*randn(size(f)),2*pi);
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
            G  = -1/2*symMod((X+Z)-2*Y,2*pi);
        end
        function ds = dot(this,P,V,W)
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
