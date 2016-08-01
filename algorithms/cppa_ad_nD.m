function [x,fctValSeq] = cppa_ad_nD(varargin)
% cppa_ad_nD(M,f,alpha,beta,lambda)
%      Compute the cyclic proximal point algorithm for first and
%      second order absolute differences, weighted by alpha and
%      beta on a manifold for n-dimensional discrete data
% cppa_ad_nD(problem) performs as the one above but all parameters are
%      provided in a struct problem for performance reasons.
%
% INPUT
%   M      ; a manifold.
%   f      : an image from M^{m_1,m_2,...,m_n} of phase valued data
%   alpha  : weight for first order absolute difference terms
%            either a number, 2 numbers (straight and
%            diagonals) or one for each
%            (i.e. n + 2(n-1 + ... + 1))
%   beta   : a weight for each weight for second order absolute difference terms
%            either two nubers (straight & mixed differences)
%            or one for each (n + n-1+...+1)
%   lambda : weight for the proximal mappings
%
%   OUTPUT:
%      x   : result of the cyclic proximal point algorithm
%
%   OPTIONAL PARAMETERS
%      'MaxIterations' : (400) Maximal number of Iterations
%      'Epsilon'       : (10^{-6}) Lower bound for the max change in one cycle
%
%           one of the parameters MaxIterations/Epsilon can be dactivated
%           by specifying them as Inf but not both.
%
% 'UnknownMask' : ([]) Specify a mask, which values of f are unknown
%                 (1) and initialized in the cycles, when possible.
%
%           If specified all Masks have to be the same size as f.
% 'RegMask'     : ([]) Specify a binary mask for values affected by
%                  (1) they are regularized, i.e. moved (or active, say)
%                  (0) fixed. RegMask is of the same size
%                      the data point set, i.e. [m_1,...,n_m].
%  'EvalFct'    : ([]) specify a function to be evaluated in every
%                          iteration and accepts two parameter, the x_k
%                          and the distance data array from x_k to x_{k-1} measured on
%                          the manifold
%  See also:
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann,
%     A second order non-smooth variational model for restoring
%     manifold-valued images,
%     SIAM Journal on Scientific Computing, 38, (1), A567?A597, 2016.
%
%    manifold.proxad, cppa_ad_1D, cppa_ad_2D
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-10-22 | 2016-02-22
% see LICENSE.txt

%
% Changelog
% ---
% 2015-05-09 Adapted to new function interfaces having manifolds first.
if (length(varargin)==1 && isstruct(varargin{1}))
    % short hand call to surpass parser
    vars = varargin{1};
    assert(all(isfield(vars,{'M','f','alpha','beta','lambda'})),'There are required input fields missing'); % at least required
    assert(isa(vars.M,'manifold'),'Input M is not a manifold');
    dimen = size(vars.f);
    manDims = length(vars.M.ItemSize);
    imgDim = dimen((manDims+1):end);
    imgDims = length(imgDim);
    manDim = dimen(1:manDims);
    % Fill optionals
    useRegMask = isfield(vars,'RegMask');
    if ~useRegMask
        rM = true(imgDim); % regularize all
    else
        rM = vars.RegMask;
    end
    useUnknownMask = isfield(vars,'UnknownMask');
    if ~useUnknownMask
        uM = false(imgDim);
    else
        uM = vars.UnknownMask;
    end
    if ~isfield(vars,'MaxIterations')
        maxiter = 400;
    else
        maxiter = vars.MaxIterations;
    end
    if ~isfield(vars,'Epsilon')
        epsilon = 10^(-6);
    else
        epsilon = vars.Epsilon;
    end    % Fill optionals
    if ~isfield(vars,'EvalFct')
        evalFct = [];
    else
        evalFct = vars.EvalFct;
    end
    alphanum = imgDims + (imgDims-1)*imgDims;
    alpha = vars.alpha;
    betanum = imgDims + (imgDims-1)*imgDims/2;
    beta = vars.beta;
else % Parse Variables
    ip = inputParser;
    addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
    addRequired(ip,'f');
    addRequired(ip,'alpha');
    addRequired(ip,'beta');
    addRequired(ip,'lambda');
    addParameter(ip,'RegMask', []);
    addParameter(ip,'UnknownMask',[]);
    addParameter(ip,'MaxIterations',400);
    addParameter(ip,'Epsilon',10^(-6));
    addParameter(ip,'EvalFct',[]);
    parse(ip, varargin{:});
    vars = ip.Results;
    %
    % Validate optional stuff
    % a) Iteratiosn & Epsilon
    if sum(size(vars.MaxIterations))>2
        error('The maximal number of iterations has to be a number.');
    end
    if floor(vars.MaxIterations)~=vars.MaxIterations
        error('The maximal number of iterations has to be an integer.');
    end
    maxiter = vars.MaxIterations;
    if sum(size(vars.Epsilon))>2
        error('The upper bound for cyclic step change has to be a number');
    end
    epsilon = vars.Epsilon;
    % Error checks
    if ( isinf(maxiter) && (isinf(epsilon)) )
        error('You can not set both the maximal Iterations and Epsilon to Inf');
    end
    % (b) Regularization parameters
    dimen = size(vars.f);
    manDims = length(vars.M.ItemSize);
    imgDim = dimen(manDims+1:end);
    imgDims = length(imgDim);
    manDim = dimen(1:manDims);
    % Number of alphas
    alphanum = imgDims + (imgDims-1)*imgDims;
    assert(isvector(vars.alpha),...
        ['alpha has to be a number or a vector of length 2 or of length ',num2str(alphanum),'.']);
    assert(all(vars.alpha>=0) && all(vars.beta>=0),...
        'All entries of alpha and beta have to be nonnegative.');
    betanum = imgDim + (imgDim-1)*imgDim/2;
    assert(isvector(vars.beta),...
        ['beta has to be a number or a vector of length 2 or of length ',num2str(betanum),'.']);
    %
    % Handle optional input RegMask
    useRegMask = false;
    rM = true(imgDim); %Standard: All are regularized, none is fixed
    if (numel(vars.RegMask)>0) % Mask given
        if ( any(size(vars.UnknownMask) ~= m) )
            warning(['Size of the regularization mask (',...
                num2str(size(vars.RegMask)),') does not fit the size of image f (',...
                num2str(m),'). Hence it is ignored.']);
        else
            rM = vars.UnknownMask;
            useRegMask =true;
        end
    end
    useUnknownMask = false;
    uM = false(imgDim); %Standard: All are known
    if (numel(vars.UnknownMask)>0) % Mask given
        if ( any(size(vars.UnknownMask) ~= m) )
            warning(['Size of the unknown data mask (',...
                num2str(size(vars.UnknownMask)),') does not fit the size of image f (',...
                num2str(imgDim),'). Hence it is ignored.']);
        else
            uM = vars.UnknownMask;
            useUnknownMask =true;
        end
    end
end
% Expand alpha/beta
    if (length(alpha)==1) % one number fits all
        talpha = alpha;
        alpha = zeros(1,alphanum);
        alpha(1:imgDims) = talpha;
    elseif length(alpha)==2 % set the first ones
        talpha = alpha;
        alpha = zeros(1,alphanum);
        alpha(1:imgDims) = talpha(1);
        alpha((imgDims+1):end) = talpha(2);
    elseif length(vars.alpha)~=alphanum
        error(['Length of alpha (',num2str(length(alpha)),' does not fit data size, i.e. is none of the values 1,2,',num2str(alphanum),'.']);
    end
    if (length(beta)==1) % one number fits all
        beta = ones(1,betanum)*beta;
    elseif length(vars.beta)==2 % set the first ones
        tbeta = beta;
        beta = zeros(1,betanum);
        beta(1:imgDims) = tbeta(1);
        beta((imgDims+1):end) = tbeta(2);
    elseif length(vars.beta)~=betanum
        error(['Length of beta (',num2str(length(beta)),' does not fit data size, i.e. is none of the values 1,2,',num2str(betanum),'.']);
    end           % N is the number of data points, k the dimension of the manifold, and they

% --- INIT x
x = vars.f;
%
% Handle optional input UnknownMask
if (useUnknownMask) %kill unknown pixels completely
    x(repmat(uM==1,[manDim,ones(size(imgDim))])) = NaN;
end
if useRegMask&&useUnknownMask
    fixedAndUnknown = uM&(~rM);
    assert( ~any(fixedAndUnknown(:)), 'There exists at least one pixel that is unknown and fixed by the two masks, which is not allowed.');
end
% init values for the loop
recordFct = (nargout==2) && ( sum(size(evalFct)) > 0) && isa(evalFct,'function_handle');
if recordFct
    fctValSeq = [];
end
i=0;
lambdait=vars.lambda;
itD = [Inf, NaN]; %iteration Distances
debug('time',3,'StartTimer','Cyclic proximal point algorithm on manifolds');
% Help shifters
smfc = [3,1,2];
smlc = [2,3,1];
while ( (any(isnan(itD(:))) || (max(itD(~isnan(itD)))>=epsilon)) && (i<maxiter))
    % work on manifold first arrays
    xold = x;
    % First step : proximal mappings of the distances - only on
    % known points
    Inds = ones(size(x));
    if useRegMask
        Inds = Inds&repmat(permute(rM==1,[imgDims+(1:manDims),1:imgDims]),...
            [manDim,ones(1,imgDims)]);
    end
    if useUnknownMask
        Inds = Inds&repmat(permute(~uM,...
            [imgDims+(1:manDims),1:imgDims]),...
        [manDim,ones(1,imgDims)]);
    end
    if useRegMask || useUnknownMask
        % Prox only existing/affected
        x(Inds) = vars.M.proxDist(reshape(x(Inds),[vars.M.ItemSize,numel(x(Inds))/prod(vars.M.ItemSize)]),...
                reshape(vars.f(Inds),[vars.M.ItemSize,numel(vars.f(Inds))/prod(vars.M.ItemSize)]),lambdait);
    else
        x = vars.M.proxDist(x,vars.f,lambdait);
    end
    % (A) internally for reshaping use collapsed mainfold dims
    x = reshape(x,[prod(manDim),imgDim]); %linearize manifold dimensions into first dim
    for l4=0:1 % first and second split term for first differences
        for l = 1:imgDims % each dimension
            mlt = floor( (imgDim(l) - l4)/2 );
            remnum = prod(imgDim(1:(l-1)))*prod(imgDim(l+1:end));
            % finite difference along x_l
            if alpha(l) > 0
                % (a) permute the subsampled dimension to second (most inner
                % permute, not necessary above)
                p1t = 1:imgDims;
                p1 = [1,l+1, 1+p1t((p1t~=l))];
                invp1 = zeros(size(p1));
                invp1(p1) = 1:(imgDims+1);
                xshift = permute(x,p1); % shift manifold to first, l to second               
                % (b) reshape to manifoldx2xRestsize(
                % (c) permute 2 to third dim for easier proxes
                problem.f = reshape(xshift(:,(l4+1):(2*mlt+l4),:),[manDim,remnum*mlt,2]);
                problem.lambda = alpha(l)*lambdait;
                if useRegMask
                    rMshift = permute(rM,p1(2:end)-1); %same shift in rM but ignoring manifold dimension
                    problem.regMask = reshape(rMshift((l4+1):(2*mlt+l4),:),[2,remnum*mlt]);
                end
                problem.w = [-1,1]';
                xshift(:,(l4+1):(2*mlt+l4),:) ... %1: manifold - keep, 2: difference dim - sub, 3-end keep in one
                    = reshape(vars.M.proxad(problem),[prod(manDim),mlt*2,remnum]); %reshape to l-dims, rest, manifold and shift manifold first again
                x = permute(xshift,invp1); % undo shift from beginning of for loop
            end
        end
        % All Diagonals
        alphaind=imgDims;
        for l1 = 1:imgDims
            for l2 = (l1+1):imgDims
                ml1 = floor((imgDim(l1)-l4)/2); %Mt
                ml2 = floor((imgDim(l2)-l4)/2); %Nt
                for l3=0:1 % both diagonals
                    alphaind = alphaind+1;
                    if (alpha(alphaind)>0)
                        % shift manifold 1, l1 to 2, l2 to 3, rest ordered
                        p1t = 1:imgDims;
                        p1 = [1,l1+1,l2+1, 1+p1t( (p1t~=l1)&(p1t~=l2) )];
                        invp1 = zeros(size(p1));
                        invp1(p1) = 1:(imgDims+1); %inverse
                        xshift = permute(x,p1);
                        problem.f = cat(2,...
                        reshape(xshift(:,(1+l4):2:2*ml1+l4,(1+l3):ml2-1+l3,:),prod(manDim),1,[]),... %all even rows,
                        reshape(xshift(:,(2+l4):2:2*ml1+l4,(2-l3):ml2-l3,:),prod(manDim),1,[])); %all odd rows
                        problem.lambda = alpha(alphaind)*lambdait;
                        problem.w = [-1,1];
                        if useRegMask
                            rMshift = permute(rM,p1(2:end));
                            problem.RegMask = [reshape(rMshift((1+l4):2:2*ml1+l4,(1+l3):ml2-1+l3,:),1,[]);... %all even rows,
                        reshape(rMshift((2+l4):2:2*ml1+l4,(2-l3):ml2-l3,:),1,[]) ]; %all odd rows
                        end
                        res = vars.M.proxad(problem);
                        % rearrange manifold first again
                        xshift(:,(1+l4):2:2*ml1+l4,(1+l3):ml2-1+l3,:) = reshape(res(:,1,:), size(xshift(:,(1+l4):2:2*ml1+l4,(1+l3):ml2-1+l3,:)) );
                        xshift(:,(2+l4):2:2*ml1+l4,(2-l3):ml2-l3,:) = reshape(res(:,2,:), size(xshift(:,(2+l4):2:2*ml1+l4,(2-l3):ml2-l3,:)) );
                        % undo the first shift
                        x = permute(xshift,invp1);
                    end
                end
            end
        end
    end % end j, both splitting blocks for first differnces
    %
    % Third: Second order differences
    for l4=0:2 % all splits
        % Xl
        for l = 1:imgDims % each dimension
            mlt = floor( (imgDim(l) - l4)/3 );
            remnum = prod(imgDim(1:(l-1)))*prod(imgDim(l+1:end));
            % finite difference along x_l
            if beta(l) > 0
                p1t = 1:imgDims;
                p1 = [1,l+1, 1+p1t((p1t~=l))];
                invp1 = zeros(size(p1));
                invp1(p1) = 1:(imgDims+1);
                xshift = permute(x,p1); % shift l to second
                problem.f = reshape(xshift(:,(l4+1):(3*mlt+l4),:),[prod(manDim),3,remnum*mlt]);
                problem.lambda = beta(l)*lambdait;
                problem.w = [1,-2,1];
                if useRegMask
                    rMshift = permute(rM,p1(2:end));
                    problem.RegMask = reshape(rMshift((l4+1):(3*mlt+l4),:),[3,remnum*mlt]);
                end
                xshift(:,(l4+1):(3*mlt+l4),:) ... %1: manifold - keep, 2: difference dim - sub, 3-end keep in one
                  = reshape(vars.M.proxad(problem),[prod(manDim),mlt*3,remnum]);
                x = permute(xshift,invp1); % undo shift from beginning of for loop
            end
        end
    end
    betaind = imgDims;
    for l1 = 1:imgDims
        for l2 = (l1+1):imgDims
            betaind = betaind + 1;
            if (beta(betaind) > 0)
                ml1 = [2*floor(m(l1)/2),2*floor((m(l1)-1)/2)]; % NT
                ml2 = [2*floor(m(l2)/2),2*floor((m(l2)-1)/2)]; % MT
                % 3 shifts as above
                p1t = 1:n;
                p1 = [1,l1+1,l2+1, 1+p1t( (p1t~=l1)&(p1t~=l2) )];
                invp1 = zeros(size(p1));
                invp1(p1) = 1:(n+1); %inverse
                %for adressing: manifold first, before
                xshift = permute(x,p1);
                rMshift = permute(rM,p1(2:end));
                for l3=0:1
                    for l4=0:1
                        %reshape with manifold last
                        problem.f = cat(2,...
                            reshape(xshift(:,1+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:),prod(manDim),1,[]),...
                            reshape(xshift(:,2+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:),prod(manDim),1,[]),...
                            reshape(xshift(:,2+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:),prod(manDim),1,[]),...
                            reshape(xshift(:,1+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:),prod(manDim),1,[]));
                        problem.lambda = beta(betaind)*lambdait;
                        problem.w = [-1,1,-1,1];
                        if useRegMask
                        problem.RegMask = [reshape(rMshift(1+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:),1,[]);...
                            reshape(rMshift(1+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:),1,[]);...
                            reshape(rMshift(2+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:),1,[]);...
                            reshape(rMshift(2+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:),1,[])];
                        end
                        tx = vars.M.proxad(problem);
                        %reset to xshift, i.e. permute right to
                        %manifold first
                        xshift(:,1+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:) = reshape(tx(:,:,1),size(xshift(:,1+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:)));
                        xshift(:,2+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:) = reshape(tx(:,:,2),size(xshift(:,2+l3:2:ml1(l3+1)+l3, 1+l4:2:ml2(l4+1)+l4,:)));
                        xshift(:,2+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:) = reshape(tx(:,:,3),size(xshift(:,2+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:)));
                        xshift(:,1+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:) = reshape(tx(:,:,4),size(xshift(:,1+l3:2:ml1(l3+1)+l3, 2+l4:2:ml2(l4+1)+l4,:)));
                    end
                end
                % shift back 3
                x = permute(xshift,invp1);
            end
        end
    end
    i = i + 1;
    lambdait = vars.lambda/(i);
    x = reshape(x,[manDim,imgDim]); %expand manifold dimensions from first dim
    itD = vars.M.dist(x,xold);
    if recordFct
        fctValSeq = [fctValSeq,evalFct(x,itD)]; %#ok<AGROW>
    end
    %
    % shift back
    debug('text',3,'text',...
        ['i=',num2str(i),' where lastdiff is ',num2str(max(itD(~isnan(itD)))),' and lambdait=',num2str(lambdait)]);
end
debug('time',3,'StopTimer','Cyclic proximal point algorithm on manifolds');
debug('text',2,'Text',...
    [num2str(i),' iterations, last difference:',...
    num2str(max(itD(~isnan(itD))))]);
end
