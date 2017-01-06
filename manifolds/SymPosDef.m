classdef SymPosDef < manifold & handle
    % The manifold of m-by-m symmetric positive definite matrices
    %   which is denoted P(m) for short
    %
    % Additionally provides
    %
    % PROPERTIES
    %    tau   : stepsize in the subgradient descent inside each of the
    %            proximal mappings of the second order differences
    %    steps : number of steps in the just mentioned subgradient descent.
    %
    % FUNCTIONS
    %    parallelTransport(X,Z) : Compute the mid point on the geodesic
    %    dot(X,V,W)             : Riemannian dot product of V,W in TXS(m)
    %    JacobianEigenFrame(X,V): A frame for TXP(m) based on (but not
    %                             including) V.
    %
    % ---
    % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-01-28
    properties
        type = 'SymPosDef';
        tau = 1;
        steps = 10;
        ItemSize;
        Dimension;
    end
    
    methods
        function obj = SymPosDef(m)
            % Create a Manifold for mxm symmetric positive definite matrices
            obj.type = ['symmetric positive definite ',num2str(m),'x',num2str(m),' matrices'];
            obj.ItemSize = [m,m];
            obj.Dimension = m*(m+1)/2;
        end
        function x = mean(this,varargin)
            % mean(f) calculates the mean of the input data with a gradient
            % descent algorithm. This implementation is based on
            %
            % B. Afsari, Riemannian Lp center of mass: Existence,
            %    uniqueness, and convexity,
            %    Proc. AMS 139(2), pp.655-673, 2011.
            %
            % and this method overwrites the general case, because it is
            % implemented in C++.
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
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch 2015-11-05
            if ~this.useMex || getDebugLevel('MatlabMean')==1 % fallback nonMex
                x=this.mean@manifold(varargin{:});
            else
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
                    epsilon = 10^-12;
                end
                if vars.MaxIterations > 0
                    iter = vars.MaxIterations;
                else
                    warning('Iterations should be larger than zero, set Iterations to 100')
                    iter = 100;
                end
                x = SPDMean(f,w,epsilon,iter);
            end
        end
        function Y = exp(this,X,V,t)
            % exp(X,V) - Exponential map at the point(s) P with respect to
            %    the direction(s) V in TXP(m)
            %
            % INPUT
            %   X : a point or set of points on the manifold P(m)
            %   V : a point or set of point in the tangential spaces TpM
            %
            % OPTIONAL
            %   t : shorten vectors V by factor t
            %   (given as one value or an array of same length as number of Vs)
            %
            % OUTPUT
            %   Y : resulting point(s) on P(m)
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-01-20 | 2015-04-10
            
            % Changelog
            %   2015-04-10 introduced mexfile
            if nargin<4 % define t
                t=1;
            end
            if this.useMex
                Y = SPDExp(X,V,t);
            else
                Y = this.localExp(X,V,t);
            end
        end
        function V = log(this,X,Y)%,verifyFlag)
            % log(X,Y) - Exponential map at the point(s) X of points(s) Y
            %      in TXP(m)
            %
            % INPUT
            %   X : a point or set of points on the manifold P(m)
            %   Y : a point or set of points on the manifold P(m)
            %
            % OUTPUT
            %   V : resulting point(s) of X(i) to Y(i) elementwise
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-01-20 | 2015-04-10
            
            % Changelog
            % 2015-04-10 Changed to Mex.
            if this.useMex
                V =  SPDLog(X,Y);
            else
                V = this.localLog(X,Y);
            end
        end
        function d = dist(this,X,Y)
            % d = dist(X,Y) - Compute the distance between points or
            % two sets of points on the manifold P(m).
            %
            % INPUT
            %   X,Y : two points (matrices) or sets of points (matrices)
            %         on P(m)
            %
            % OUTPUT
            %     d : resulting distances of each pair of points of p,q.
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2014-10-19 | 2015-04-11
            
            % Changelog
            %   2015-04-11 Added mex-file version
            % TODO: DEbug dist
            if this.useMex
                d = SPDDist(X,Y);
            else
                d = this.localDist(X,Y);
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
            if this.useMex
                W = SPDParallelTransport(X,Y,V);
            else
                W = this.localParallelTransport(X,Y,V);
            end
        end
        function ds = dot(this,X,V,W)
            % PnDot(X,V,W)
            %     Compute the inner product of two tangent vectors in TXP(m)
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
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-01-29 | 2015-04-11
            
            % Changelog
            %   2015-04-11 Created Mex-File
            if this.useMex
                ds = SPDDot(X,V,W);
            else
                ds = this.localDot(X,V,W);
            end
        end
        function [lambda,W] = JacobianEigenFrame(this,X,V)
            % [lambda,W] = JacobianEigenFrame(X,V) Compute the eigenframe and
            % eigenvalues for the (each) tangential plane of the point( set) X
            % w.r.t (one of each) V.
            %
            % INPUT
            %    X : A point( set)
            %    V : A tangential vector( set)
            %
            % OUTPUT
            %     lambda : eigenvalue( set)s of the riemannian tensor at (each) X
            %     W      : correspoinding eigenmatrices
            %
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-01-20 | 2015-04-11
            
            % Changelog
            %   2015-04-11 Extracted to local; optimized GradX instead of
            %              writing a Mex here o spare memory
            [lambda,W] = this.localJacobianEigenFrame(X,V);
        end
        function G = grad_X_D2(this,X,Y,Z)
            % grad_X_D2(X,Y,Z) Compute the gradient with respect to the first
            % variable of the second order difference term. This can also
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
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann  2015-01-29
            % Compute mid point(s)
            if this.useMex
                G = SPDGrad_X_D2(X,Y,Z);
            else
                G = this.localGrad_X_D2(X,Y,Z);
            end
        end
        function G = grad_X_D2_Sq(this,X,Z,Y)
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
            % ManImRes 1.0, R. Bergmann  2015-01-29
            % Compute mid point(s)
            if this.useMex
                G = SPDGrad_X_D2_Sq(X,Y,Z);
            else
                G = this.localGrad_X_D2_Sq(X,Y,Z);
            end
        end
        function x = proxad(this,varargin)
            % proxad(f, lambda, w) - Compute the proximal mapping of the
            % data f with respect to lambda an the weight w.
            % Any unknown point containing a NaN is inpainted if
            % neighbouring known points or else left unchanged
            %
            % INPUT
            %      f : data on P(m) (i.e. mxm), where [m,m,l,d] = size(f) and d
            %          is the number of data points per prox and l is the
            %          number of parallel proxes. A point is unknown if for
            %          fixed d,l any of the three entries is NaN.
            % lambda : weight of the proximal map (one value)
            %      w : finite difference weight vector. For P(m) up to now
            %          only [-1,1] and [1,-2,1] are supported
            % OPTIONAL PARAMETERS
            %     'UnknownMask'  : ([]) Specify a mask, which values of f are unknown
            %                      (1) and initialized in the cycles otherwise, when
            %                         possible (depending on the neighbours).
            %     'RegMask'      : ([]) Specify a binary mask for values affected by
            %                     (1) they are regularized, i.e. moved (or active, say)
            %                     (0) fixed. RegMask is of the same size
            %                         the data point set, i.e. [n,m].
            %      If specified all Masks have to be the same size as f.
            %
            % OUTPUT
            %      x : resulting data of all proximal maps (mxmxdxl)
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-01-20
            if (length(varargin)==1 && isstruct(varargin{1})) %struct case
                vars = varargin{1};
                assert(isfield(vars,'f')&&isfield(vars,'lambda')&&isfield(vars,'w'),...
                    'Not all required parameters given in the struct');
                if ~isfield(vars,'RegMask')
                    vars.RegMask = [];
                end
                d = size(vars.f,4);
                l = size(vars.f,3);
                w = vars.w; %no check!
            else
                ip = inputParser;
                addRequired(ip,'f');
                addRequired(ip,'lambda');
                addRequired(ip,'w');
                addParameter(ip, 'RegMask', []);
                parse(ip, varargin{:});
                vars = ip.Results;
                [a,b,l,d] = size(vars.f);
                if any([a,b]~=this.ItemSize)
                    error(['The values of f (',num2str(d),' points in ',...
                        num2str(l),' proximal mappings of dimension ',num2str(a),'x',num2str(b),...
                        ') are not ',num2str(this.ItemSize(1)),'x',num2str(this.ItemSize(2)),' matrices.']);
                end
                if (isrow(vars.w))
                    w = vars.w';
                else
                    w = vars.w;
                end
            end
            if (vars.lambda==0)
                x = vars.f;
                return
            end
            if sum(size(vars.RegMask))==0
                vars.RegMask = ones(l,d);
            end
            assert(d==length(w),['The length of the weight (',...
                num2str(length(w)),') does not fit to the data size (',...
                num2str(d),').']);
            x = vars.f;
            x(isnan(vars.f)) = NaN;
            missingpoints = squeeze(any(any(isnan(vars.f),1),2)); %sum along manifold dimensions & squeeze
            if l==1
                missingpoints = missingpoints';
            end
            % proxes that can be computed
            proxpts = sum(missingpoints,2)==0; %last dimension d
            proxpts = proxpts & (sum(vars.RegMask,2)==d); %only affected proxes
            
            % proxes that can be inpainted
            inpaintpts = sum(missingpoints,2)==1;
            if (length(w)==2) && all(w==[-1,1]')
                vecls = this.dist(vars.f(:,:,proxpts,1), vars.f(:,:,proxpts,2));
                t = zeros(size(vecls));
                t(vecls>eps) = min(vars.lambda./(vecls(vecls>eps)), 0.5); %Not a unit speed geodesic -> divide by length
                % Compute new values
                x(:,:,proxpts,1) = this.exp(vars.f(:,:,proxpts,1), repmat(reshape(t,[ones(1,length(this.ItemSize)),length(t)]),[this.ItemSize,1]).*this.log(vars.f(:,:,proxpts,1),vars.f(:,:,proxpts,2)));
                x(:,:,proxpts,2) = this.exp(vars.f(:,:,proxpts,2), repmat(reshape(t,[ones(1,length(this.ItemSize)),length(t)]),[this.ItemSize,1]).*this.log(vars.f(:,:,proxpts,2),vars.f(:,:,proxpts,1)));
                %
                % Inpaint all first missing points, the second one is
                % existing due to limiting the missing number to 1
                x(:,:,(missingpoints(:,1)>0)&inpaintpts,1) = x(:,:,(missingpoints(:,1)>0)&inpaintpts,2);
                x(:,:,(missingpoints(:,2)>0)&inpaintpts,2) = x(:,:,(missingpoints(:,2)>0)&inpaintpts,1);
            elseif (length(w)==3) && all(w==[1,-2,1]')
                %
                % Iterative subgradient descent
                x = vars.f;
                % Gradient descent
                xopt = x(:,:,proxpts,:);
                xoptvals = vars.lambda*this.dist(this.midPoint(x(:,:,proxpts,1),x(:,:,proxpts,3)),x(:,:,proxpts,2)); %x=f hence first term zero
                for i=1:this.steps
                    % midpoints between firsts and thirds
                    % Gradient X
                    grad1 = this.log(x(:,:,proxpts,1),vars.f(:,:,proxpts,1)) + vars.lambda*this.grad_X_D2(x(:,:,proxpts,1),x(:,:,proxpts,2),x(:,:,proxpts,3));
                    % Gradient Y (easy)
                    M = this.exp(x(:,:,proxpts,1), this.log(x(:,:,proxpts,1),x(:,:,proxpts,3))/2);
                    V = this.log(x(:,:,proxpts,2),M);
                    Vl = sqrt(this.dot(x(:,:,proxpts,2),V,V));
                    % rearrange dimensions, let ./ do the singleton
                    % expansion
                    V(:,:,Vl>eps) = V(:,:,Vl>eps)./...
                        repmat(reshape(Vl(Vl>eps),[ones(1,length(this.ItemSize)),length(Vl(Vl>eps))]),[this.ItemSize,1]);
                    V(:,:,~(Vl>eps)) = 0;
                    grad2 = this.log(x(:,:,proxpts,2),vars.f(:,:,proxpts,2)) + vars.lambda*V;
                    % Gradient Z
                    grad3 = this.log(x(:,:,proxpts,3),vars.f(:,:,proxpts,3)) + vars.lambda*this.grad_X_D2(x(:,:,proxpts,3),x(:,:,proxpts,2),x(:,:,proxpts,1));
                    tauit = this.tau/i;
                    x(:,:,proxpts,1) = this.exp(x(:,:,proxpts,1), tauit*grad1);
                    x(:,:,proxpts,2) = this.exp(x(:,:,proxpts,2), tauit*grad2);
                    x(:,:,proxpts,3) = this.exp(x(:,:,proxpts,3), tauit*grad3);
                    xoptt = x(:,:,proxpts,:);
                    xvals = 1/2*sum(this.dist(vars.f(:,:,proxpts,:),xoptt).^2,2) + vars.lambda*this.dist(this.midPoint(xoptt(:,:,:,1),xoptt(:,:,:,3)),xoptt(:,:,:,2));
                    xopt(:,:,xvals<xoptvals,:) = xoptt(:,:,xvals<xoptvals,:);
                    xoptvals(xvals<xoptvals) = xvals(xvals<xoptvals);
                    x(:,:,proxpts,:) = xopt;
                end
                %we do not interpolate up to now!
            elseif (length(w)==4) && all(w==[-1,1,-1,1]')
                x=vars.f;
                xopt = x(:,:,proxpts,:);
                xoptvals = vars.lambda*this.dist(this.midPoint(x(:,:,proxpts,1),x(:,:,proxpts,3)),this.midPoint(x(:,:,proxpts,2),x(:,:,proxpts,4))); %x=f hence first term zero
                for i=1:this.steps
                    % Two arrays of midpoints
                    M = zeros(this.ItemSize(1),this.ItemSize(2),size(x(1,1,proxpts,1),3),2);
                    M(:,:,:,1) = this.exp(x(:,:,proxpts,1), this.log(x(:,:,proxpts,1),x(:,:,proxpts,3))/2);
                    M(:,:,:,2) = this.exp(x(:,:,proxpts,2), this.log(x(:,:,proxpts,2),x(:,:,proxpts,4))/2);
                    grad = zeros(this.ItemSize(1),this.ItemSize(2),size(M,3),4); % all four gradients
                    for j=1:4
                        j2 = mod(j+1,4)+1; %other index involved in same mid point, 1->3, 2->4, 3->1, 4->2
                        mi = mod(j,2)+1; %opposite midpoint index (1,3)->2, (2,4)->1
                        % similar to the approach for 1,-2,1 just that y is
                        % the other middle point
                        grad(:,:,:,j) = this.grad_X_D2(x(:,:,proxpts,j),M(:,:,:,mi),x(:,:,proxpts,j2));
                    end
                    %perform a gradient step
                    tauit = this.tau/i;
                    for j=1:4
                        x(:,:,proxpts,j) = this.exp(x(:,:,proxpts,j), tauit*grad(:,:,:,j));
                    end
                    xoptt = x(:,:,proxpts,:);
                    xvals = sum(this.dist(vars.f(:,:,proxpts,:),xoptt).^2,2)/2 + vars.lambda*this.dist(this.midPoint(xoptt(:,:,:,1),xoptt(:,:,:,3)),this.midPoint(xoptt(:,:,:,2),xoptt(:,:,:,4)));
                    xopt(:,:,xvals<xoptvals,:) = xoptt(:,:,xvals<xoptvals,:);
                    xoptvals(xvals<xoptvals) = xvals(xvals<xoptvals);
                    x(:,:,proxpts,:) = xopt;
                end
            else
                warning(['Unknown discrete difference on ',this.type,': ',num2str(w'),'. Returning the input f.']);
                x=vars.f;
            end
        end
        function Y = addNoise(this,X,sigma,varargin)
            ip = inputParser;
            addParameter(ip,'Distribution','Rician');
            parse(ip,varargin{:});
            vars = ip.Results;
            sizes = size(X);
            manDim = sizes(1:length(this.ItemSize));
            assert(all(manDim==this.ItemSize),...
                ['Manifold Dimnsions of signal X (',num2str(manDim)...
                ,') do not fit to this manifold (',num2str(this.ItemSize),').']);
            dataDim = sizes((length(this.ItemSize)+1):end);
            if strcmp(vars.Distribution,'Rician')
                X = reshape(X,[manDim,prod(dataDim)]);
                Y = zeros(size(X));                
                for l=1:prod(dataDim)
                    T = chol(X(:,:,l)) + triu( sqrt(sigma)*randn(this.ItemSize));
                    Y(:,:,l) = T*permute(T,[2,1]);
                end
                Y = reshape(Y,sizes);
            elseif strcmp(vars.Distribution,'Gaussian')
                % Should be done in C
                ONB_Id = zeros([this.ItemSize,this.Dimension]);
                n = this.ItemSize(1);
                n_sq = n^2;
                ONB = reshape(eye(n_sq),n,n,[]);
                ONB_Id(:,:,1:n) = ONB(:,:,(1:n:(n_sq-n+1))+(0:(n-1)));
                ONB = ONB(:,:,find(tril(ones(n),-1)));
                ONB_Id(:,:,n+1:end) = 1/sqrt(2)*(ONB+permute(ONB,[2,1,3]));
                X = reshape(X,[manDim,prod(dataDim)]);
                V = sum(repmat(sigma*randn(1,1,prod(dataDim),this.Dimension),n,n,1,1).*repmat(permute(ONB_Id,[1,2,4,3]),1,1,prod(dataDim),1),4);
                Y = this.exp(X,this.parallelTransport(repmat(eye(n),1,1,prod(dataDim)),X,V));
                Y = reshape(Y,sizes);
            else
                warning('Distribution unknown, output equals input');
                Y = X;
            end
        end
        function [ONB] = TpMONB(this,p)
            % V = TpMONB(p)
            % Compute an ONB in TpM
            %
            % INPUT
            %     p : base point( sets)
            % OUTPUT
            %    ONB : orthonormal bases ( n x n x SetDim x n(n+1)/2 )
            % ---
            % ManImRes 1.0, J. Persch ~ 2016-06-13
            dimen = size(p);
            [~,ONB] = this.localJacobianEigenFrame(p,this.log(p,repmat(eye(dimen(1)),[1,1,dimen(3:end)])));
            ONB = permute(ONB,[1:2, 3+(1:max(length(dimen(3:end)),1)), 3]);
        end
    end
    methods (Access = private)
        %
        % local functions as fallback from Mex
        %
        function d = localDist(~,X,Y)
            % Instead of sqrt(trace(log(sqrt(X)^-1 Y sqrt(X)^-1)^2)) we can
            % compute the generalized eigenvalues of X and Y, i.e., 
            % l_i X x_i = Y x_i and take sqrt(sum(log^2(l_i)) idea from
            % Torben Fetzer
            dims = size(X);
            X = reshape(X,dims(1),dims(2),[]);
            Y = reshape(Y,dims(1),dims(2),[]);
            d = zeros(1,size(X,3));
            for i=1:size(X,3)
                S =real(eig(X(:,:,i),Y(:,:,i)));
                if any(S(:)<=eps)
                    d(i) = 0; % Either X or Y not SPD
                else
                    d(i) = sqrt(sum(abs(log(S)).^2));
                end
            end
            if length(dims)==3
                d = reshape(d,dims(3),1);
            elseif length(dims)>3
                d = reshape(d,dims(3:end));
            end
        end
        function ds = localDot(~,X,V,W)
            dims = size(X);
            X = reshape(X,dims(1),dims(2),[]);
            V = reshape(V,dims(1),dims(2),[]);
            W = reshape(W,dims(1),dims(2),[]);
            ds = zeros(size(X,3),1);
            for i=1:size(X,3)
                ds(i) = trace(V(:,:,i)/X(:,:,i) * W(:,:,i)/X(:,:,i));
            end
            if length(dims)==3
                ds = reshape(ds,dims(3),1);
            elseif length(dims)>3
                ds = reshape(ds,dims(3:end));
            end
        end
        function Y = localExp(~,X,V,t)
            if isrow(t)
                pt = t';
            end
            dims = size(X);
            X = reshape(X,dims(1),dims(2),[]);
            V = reshape(V,dims(1),dims(2),[]);
            pt = pt(:);
            %put vector dimension of t to the array dimension of V (i.e. 3)
            Vt = bsxfun(@times,V,permute(pt,[3,2,1]));
            Y = zeros(size(X));
            for i=1:size(X,3)
                [U,S,~]=svd(X(:,:,i));
                if all(all(Vt(:,:,i)==0))
                    Y(:,:,i)=X(:,:,i);
                elseif all(all(X(:,:,i)==0))
                    Y(:,:,i)=X(:,:,i);
                else
                    S= diag(sqrt(max(diag(S),0))); %Avoid rounding errors, elementwise sqrt
                    sX=U*S*U';
                    A = sX\(Vt(:,:,i))/sX;
                    [U,S]=eig(0.5*(A+A'));
                    S= diag(exp(diag(S))); %Avoid rounding errors, elementwise exp
                    Y(:,:,i) = sX*(U*S*U')*sX;
                end
            end
            Y = reshape(Y,dims);
        end
        function V = localLog(~,X,Y)
            dims = size(X);
            X = reshape(X,dims(1), dims(2), []);
            Y = reshape(Y,dims(1), dims(2), []);
            V = zeros(size(X));
            for i=1:size(X,3)
                if all(all(X(:,:,i)==0))
                    V(:,:,i)=X(:,:,i);
                else
                    [U,S,~]=svd(X(:,:,i));
                    S= diag(sqrt(max(diag(S),0))); %Avoid rounding errors, elementwise sqrt
                    sX=U*S*U';
                    S = diag(1./diag(S)); %sX^-1
                    A = U*S*U'*(Y(:,:,i))*U*S*U';
                    [U,S,~]=svd(0.5*(A+A'));
                    S= diag(log(max(diag(S),eps))); %Avoid rounding errors, elementwise log
                    V(:,:,i) = sX*U*S*U'*sX;
                end
            end
            V = reshape(V,dims);
        end
 function W = localParallelTransport(~,X,Y,V)
            dims = size(X);
            X = reshape(X,dims(1),dims(2),[]);
            Y = reshape(Y,dims(1),dims(2),[]);
            V = reshape(V,dims(1),dims(2),[]);
            W = zeros(size(X));
            for i=1:size(X,3)
                [U,S,~]=svd(X(:,:,i));
                S= diag(sqrt(max(diag(S),0)));
                sX=U*S*U';
                if all(all(V(:,:,i)==0))
                    A = V(:,:,i);
                else
                    A = sX\(V(:,:,i))/sX; % most inner term
                end
                B = sX\Y(:,:,i)/sX;
                [U,S,~]=svd(0.5*(B+B'));
                S= diag(log(max(diag(S),eps))); %log(X,Y) without outer sqrtm(x)e, they vanish in the following exp
                [U,S,] = eig((0.5*0.5*( (U*S*U') + (U*S*U')')));%first 0.5 from formula second for numerical stability
                S = diag(exp(diag(S))); % i.e. Z = U*S*U'
                
                W(:,:,i) = sX*(U*S*U')*(0.5*(A+A'))*(U*S*U')*sX;
            end
            W = reshape(W,dims);
        end
        function [lambda,W] = localJacobianEigenFrame(this,X,V)
            n = (this.ItemSize(1)+1)*this.ItemSize(2)/2;
            dims = size(X);
            X = reshape(X,dims(1),dims(2),[]);
            V = reshape(V,dims(1),dims(2),[]);
            W = zeros(dims(1),dims(2),n,size(X,3));
            lambda = zeros(n,size(X,3));
            for i=1:size(X,3)
                [lambda(:,i),W(:,:,:,i)] = this.localEigenFrame(X(:,:,i),V(:,:,i));
            end
            if length(dims)>2 % back to original shape
                W = reshape(W,[dims(1:2),n,dims(3:end)]);
                lambda = reshape(lambda,[n,dims(3:end)]);
            end
        end
        function [lambda,W] = localEigenFrame(this,X,V)
            % W = localEigenframe(X,Y) Compute the Eigenframe w.r.t. the
            %     eigenvectors of V (from TXP(m)) at X
            %
            % INPUT
            %     X : a point on P(m)
            %     V : a corresponding tangential vector
            %
            % OUTPUT
            %      W : a set of vectors forming the eigenframe
            % lambda : values characterizing the curvature (-\lambda =
            %          kappa)
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~2015-01-29
            n = this.ItemSize(1)*(this.ItemSize(1)+1)/2;
            [U,S,~]=svd(X);
            S= diag(sqrt(max(diag(S),0))); %Avoid rounding errors, elementwise sqrt
            lambda = zeros(1,n);
            W = zeros(this.ItemSize(1),this.ItemSize(2),n);
            if ~(all(S(:)==0)) %|| all(V(:)==0))
                sX=U*S*U';
                VE = 0.5*(sX\V/sX + (sX\V/sX)');
                [cpointV,D] = eig(VE);
                cpointlambda = diag(D);
                % lazy double for loop - optimize later
                s = 1;
                for i=1:this.ItemSize(1)
                    for j=i:this.ItemSize(2)
                        lambda(s) = abs(cpointlambda(i)-cpointlambda(j))*norm(V,'fro');
                        t = 0.5;
                        if (i~=j)
                            t = sqrt(t);
                        end
                        thisV = t*( cpointV(:,i)*cpointV(:,j)' + cpointV(:,j)*cpointV(:,i)');
                        W(:,:,s) =  sX*(thisV)*sX;
                        s = s+1; %also lazy
                    end
                end
            end
        end
        function G = localGrad_X_D2(this,X,Y,Z)
            M = this.exp(X, this.log(X,Z)/2);
            n = this.ItemSize(1)*(this.ItemSize(1)+1)/2;
            % Solve Eigenvalue problem in Z (in Z with direction towards X)
            [lambda, Vx] = this.JacobianEigenFrame(X,this.log(X,Z));
            Vm = zeros(size(Vx));
            for i=1:n %transport frame(s) into tangential planes of X and M
                %                Vx(:,:,i,:) = this.parallelTransport(Z,X,squeeze(Vz(:,:,i,:)));
                Vm(:,:,i,:) = this.parallelTransport(X,M,squeeze(Vx(:,:,i,:)));
            end
            R = this.log(M,Y); % log_c(x) y
            Rl = sqrt(this.dot(M,R,R)); % norm of log_c(x) y
            alpha = zeros(n,size(M,3));
            for i=1:n
                alpha(i,Rl>eps) = ...
                    this.dot( M(:,:,Rl>eps),R(:,:,Rl>eps),squeeze(Vm(:,:,i,Rl>eps)))...
                    ./(Rl(Rl>eps)); % <normed_log..., B_{ij}(T/2)>
            end
            G = zeros(this.ItemSize(1),this.ItemSize(2),size(X,3));
            for i=1:n %the gradient is a weighted sum of the tangential vectors in X
                % w.,r.t the just computed weights and the eigenvalues lambda
                nonZeroL = squeeze(lambda(i,:)>eps); % times the effect of D_vc[b_{ij}]
                if sum(nonZeroL)>0
                    G(:,:,nonZeroL) = G(:,:,nonZeroL) ...
                        + repmat(...
                        permute(...
                        (sinh(lambda(i,nonZeroL)*0.5))./sinh(lambda(i,nonZeroL)) .* alpha(i,nonZeroL), ...
                        [4,3,2,1]), [this.ItemSize,1,1]).*permute(Vx(:,:,i,nonZeroL),[1,2,4,3]); %nonzero eigenvalue
                end
                if sum(~nonZeroL)>0
                    G(:,:,~nonZeroL) = G(:,:,~nonZeroL) ...
                        + repmat(permute(0.5*alpha(i,~nonZeroL),[4,3,2,1]),[this.ItemSize,1,1]).*permute(Vx(:,:,i,~nonZeroL),[1,2,4,3]);
                    % zero eigenvalue
                end
            end
        end
        function G = localGrad_X_D2_Sq(this,X,Y,Z)
            M = this.exp(X, this.log(X,Z)/2);
            n = this.ItemSize(1)*(this.ItemSize(1)+1)/2;
            % Solve Eigenvalue problem in Z (in Z with direction towards X)
            [lambda, Vx] = this.JacobianEigenFrame(X,this.log(X,Z));
            Vm = zeros(size(Vx));
            for i=1:n %transport frame(s) into tangential planes of X and M
                %                Vx(:,:,i,:) = this.parallelTransport(Z,X,squeeze(Vz(:,:,i,:)));
                Vm(:,:,i,:) = this.parallelTransport(X,M,squeeze(Vx(:,:,i,:)));
            end
            R = this.log(M,Y); % log_c(x) y
            alpha = zeros(n,size(M,3));
            for i=1:n
                % <log..., B_{ij}(T/2)>
                alpha(i,:) = 2* this.dot( M,R,squeeze(Vm(:,:,i,:)));
            end
            G = zeros(this.ItemSize(1),this.ItemSize(2),size(X,3));
            for i=1:n %the gradient is a weighted sum of the tangential vectors in X
                % w.,r.t the just computed weights and the eigenvalues lambda
                nonZeroL = squeeze(lambda(i,:)>eps); % times the effect of D_vc[b_{ij}]
                if sum(nonZeroL)>0
                    G(:,:,nonZeroL) = G(:,:,nonZeroL) ...
                        + repmat(...
                        permute(...
                        (sinh(lambda(i,nonZeroL)*0.5))./sinh(lambda(i,nonZeroL)) .* alpha(i,nonZeroL), ...
                        [4,3,2,1]), [this.ItemSize,1,1]).*permute(Vx(:,:,i,nonZeroL),[1,2,4,3]); %nonzero eigenvalue
                end
                if sum(~nonZeroL)>0
                    G(:,:,~nonZeroL) = G(:,:,~nonZeroL) ...
                        + repmat(permute(0.5*alpha(i,~nonZeroL),[4,3,2,1]),[this.ItemSize,1,1]).*permute(Vx(:,:,i,~nonZeroL),[1,2,4,3]);
                    % zero eigenvalue
                end
            end
        end
    end
end
