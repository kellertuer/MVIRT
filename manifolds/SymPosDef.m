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
        allDims;
    end
    
    methods
        function obj = SymPosDef(m)
            % Create a Manifold for mxm symmetric positive definite matrices
            obj.type = ['symmetric positive definite ',num2str(m),'x',num2str(m),' matrices'];
            obj.ItemSize = [m,m];
            obj.Dimension = m*(m+1)/2;
            obj.allDims = repelem({':'},2);
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
            if ~this.useMex
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
        function W = geopoint(this,X,Y,t)
            % geopoint(X,V,t) - Give the point \gamma_{XY}(t) of point(s) P 
            %
            % INPUT
            %   X,Y : a point or set of points on the manifold P(m)
            %   t : a scalar or set of scalars 
            %
           %
            % OUTPUT
            %   W : resulting point(s) on P(m)
            % ---
            % Manifold-Valued Image Restoration Toolbox 1.0, J. Persch ~ 2017-03-31
            
            % Changelog
            %   2015-04-10 introduced mexfile            
            if this.useMex
                W = SPDGeo(X,Y,t);
            else
                W = this.exp(X,this.log(X,Y),t);
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

        function [x,y] = prox_midpoint_data_sq(this,varargin)
            % prox_midpoint_data_sq(v,w,f, lambda) - Compute the proximal mapping of the
            % images v,w to the data f with respect to lambda.            %
            % INPUT
            %      v,w,f : data on P(m) (i.e. mxm), where [m,m,l] = size(f) = size(v) = size(w)
            %          and l is the number of parallel proxes.
            % lambda : weight of the proximal map (one value)
            % OUTPUT
            %      x : resulting data of all proximal maps (mxmxl)
            % ---
            % ManImRes 1.0, R. Bergmann ~ 2015-01-20
            if (length(varargin)==1 && isstruct(varargin{1})) %struct case
                vars = varargin{1};
                assert(isfield(vars,'f')&&isfield(vars,'lambda')&&isfield(vars,'v')&&isfield(vars,'w'),...
                    'Not all required parameters given in the struct');
                if ~isfield(vars,'RegMask');
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
                [a,b,~] = size(vars.f);
                if any([a,b]~=this.ItemSize)
                    error(['The values of f are not ',num2str(this.ItemSize(1)),'x',num2str(this.ItemSize(2)),' matrices.']);
                end
                if any(size(vars.v)~=size(vars.w)) || any(size(vars.v)~=size(vars.f))
                    error('Inputs v,w,f should be of the same size!');
                end
            end
            if (vars.lambda==0)
                x = vars.v;
                y = vars.w;
                return
            end
            v = vars.v;
            x = v;
            w = vars.w;
            y = w;
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
                xopt(:,:,func_vals<func_optvals) = xoptt(:,:,func_vals<func_optvals);
                yopt(:,:,func_vals<func_optvals) = yoptt(:,:,func_vals<func_optvals);
                y = yopt;
                x = xopt;
                func_optvals(func_vals<func_optvals) = func_vals(func_vals<func_optvals);
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
        function [W,k] = TpMONB(this,Xp,Yp)
            % V = TpMONB(p)
            % Compute an ONB in TpM
            %
            % INPUT
            %     p : base point( sets)
            % OPTIONAL:
            %     q : directional indicator( sets) for the first vector(s).
            % OUTPUT
            %    W : orthonormal bases ( n x n x SetDim x Dimension )
            %    k : (optional) curvature coefficients (Dimension x SetDim)
            % ---
            % MVIRT R. Bergmann, 2017-12-03
            n = this.Dimension;  
            dataDim  =size(Xp);
            dataDim = dataDim(length(this.ItemSize)+1:end);  
            if isempty(dataDim)
                dataDim = 1;
            end
            X = reshape(Xp,[this.ItemSize,prod(dataDim)]);
            m = size(X,3);
            dimen = size(X);
            if nargin < 3
               V = repmat(eye(dimen(1)),[1,1,dimen(3:end)]);
            else
               Y = reshape(Yp,[this.ItemSize,prod(dataDim)]);
               V = this.log(X,Y);
            end
            W = zeros(this.ItemSize(1),this.ItemSize(2),m,n);
            k = zeros(m,n);
            for l=1:m
                [U,S,~]=svd(X(this.allDims{:},l));
                S = diag(sqrt(max(diag(S),0))); %Avoid rounding errors, elementwise sqrt
                if ~(all(S(:)==0)) 
                    sX=U*S*U';
                    VE = 0.5*(sX\V(:,:,l)/sX + (sX\V(:,:,l)/sX)');
                    [cpointV,D] = eig(VE);
                    cpointlambda = diag(D);
                    % lazy double for loop - optimize later
                    s = 1;
                    for i=1:this.ItemSize(1)
                        for j=i:this.ItemSize(2)
                            k(l,s) = -1/4*abs(cpointlambda(i)-cpointlambda(j)).^2;%*norm(V,'fro');
                            t = 0.5;
                            if (i~=j)
                                t = sqrt(t);
                            end
                            thisV = t*( cpointV(:,i)*cpointV(:,j)' + cpointV(:,j)*cpointV(:,i)');
                            W(:,:,l,s) =  sX*(thisV)*sX;
                            s = s+1; %also lazy
                        end
                    end
                end
            end
            W = reshape(W,[this.ItemSize,dataDim,this.Dimension]);
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
            else 
                pt = t;
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
                    Y(:,:,i) = 1/2*(Y(:,:,i)+ Y(:,:,i).');% Symmetrize to avoid numerical errors
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
                    V(:,:,i) = 1/2*(V(:,:,i)+V(:,:,i).');% Symmetrize to avoid numerical errors
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
                S= 0.5*diag(log(max(diag(S),eps))); %log(X,Y) without outer sqrtm(x)e, they vanish in the following exp
                S = diag(exp(diag(S))); % i.e. Z = U*S*U'
                
                W(:,:,i) = sX*(U*S*U')*(0.5*(A+A'))*(U*S*U')*sX;
            end
            W = reshape(W,dims);
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
                Vm(:,:,i,:) = this.parallelTransport(X(:,:,:),M(:,:,:),squeeze(Vx(:,:,i,:)));
            end
            R = this.log(M,Y); % log_c(x) y
            alpha = zeros(n,size(M(:,:,:),3));
            for i=1:n
                % <log..., B_{ij}(T/2)>
                alpha(i,:) = 2* this.dot( M(:,:,:),R(:,:,:),squeeze(Vm(:,:,i,:)));
            end
            G = zeros(this.ItemSize(1),this.ItemSize(2),size(X(:,:,:),3));
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
            G = reshape(G,size(X));
        end
    end
end
