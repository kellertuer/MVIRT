function [x,fctValSeq] = cppa_ad_2D(varargin)
% cppa_ad_2D(M,f,alpha,beta,lambda)
%     Compute the cyclic proximal point algorithm for first and
%     second order absolute differences, weighted by alpha and
%     beta on a manifold
% cppa_ad_2D(problem)
%     performs the same, but the parameters and optional parameters are
%     given in the struct 'problem' havinf the fields of the parameter
%     names below. These are just minimally checked for performance reasons.
%
% INPUT
%   M      : a manifold
%   f      : an image from M^{n,m} manifold valued data [mfD,n,m]
%   alpha  : weight for first order absolute difference terms
%   beta   : weight for second order absolute difference terms
%   lambda : weight for the proximal mappings
%
%   OUTPUT:
%      x   : result of the cyclic proximal point algorithm
%
%   OPTIONAL PARAMETERS
%     'MaxIterations' : (400) Maximal number of Iterations
%     'Epsilon'       : (10^{-6}) Lower bound for the max change in one cycle
%
%           one of the parameters MaxIterations/Epsilon can be dactivated
%           by specifying them as Inf but not both.
%
%     'UnknownMask'  : ([]) Specify a mask, which values of f are unknown
%                        (1) and initialized in the cycles otherwise, when
%                         possible (depending on the neighbours).
%
%           If specified all Masks have to be the same size as f.
%     'RegMask'      : ([]) Specify a binary mask for values affected by
%                     (1) they are regularized, i.e. moved (or active, say)
%                     (0) fixed. RegMask is of the same size
%                         the data point set, i.e. [n,m].
%     'EvalFct'      : ([]) specify a function to be evaluated in every
%                          iteration and accepts two parameter, the x_k
%                          and the distance image from x_k to x_{k-1}
%                          measured on the manifold
%
% See also:
%    manifold.proxad, cppa_ad_1D, cppa_ad_nD
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-10-22 | 2016-02-22
% see LICENSE.txt

%
% Changelog
% 2015-01-30 Optimized to new array structure
% 2014-12-01 Changed internal calls of proxAD from inputParser to Struct
% 2014-11-29 Extracted Algorithms from the manifold class to be more
%            independent and added structs

if (length(varargin)==1 && isstruct(varargin{1}))
    % short hand call to surpass parser
    vars = varargin{1};
    assert(all(isfield(vars,{'M','f','alpha','beta','lambda'})),'There are required input fields missing'); % at least required
    assert(isa(vars.M,'manifold'),'Input M is not a manifold');
    dimen = size(vars.f);
    manDim = dimen(1:length(vars.M.ItemSize));
    imgDim = dimen(length(manDim)+1:end);
    if all(size(imgDim) == [1 1])
       imgDim = [imgDim,1]; 
    end
    % Fill optionals
    useRegMask = isfield(vars,'RegMask');
    if ~useRegMask
        rM = ones(imgDim);
    else
        rM = vars.RegMask;
    end
    useUnknownMask = isfield(vars,'UnknownMask');
    if ~useUnknownMask
        uM = zeros(imgDim);
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
    end
    if ~isfield(vars,'EvalFct')
        evalFct = [];
    else
        evalFct = vars.EvalFct;
    end
    if length(vars.alpha)>1
        alpha = zeros(1,4);
        alpha(1:length(vars.alpha)) = vars.alpha;
    else
        alpha = vars.alpha.*[1,1,0,0];
    end
    beta = vars.beta.*[1,1,1]; %just the simplest extension
else % Parse Variables and check input
    ip = inputParser;
    addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
    addRequired(ip,'f');
    addRequired(ip,'alpha');
    addRequired(ip,'beta');
    addRequired(ip,'lambda');
    addParameter(ip,'DataMask', []);
    addParameter(ip,'RegMask', []);
    addParameter(ip,'UnknownMask',[]);
    addParameter(ip,'EvalFct',[]);
    addParameter(ip,'Border',0);
    addParameter(ip,'MaxIterations',400);
    addParameter(ip,'Epsilon',10^(-6));
    parse(ip, varargin{:});
    vars = ip.Results;
    %
    % Validate Optional Stuff
    epsilon = vars.Epsilon;
    maxiter = vars.MaxIterations;
    if sum(size(vars.MaxIterations))>2
        error('The maximal number of iterations has to be a number.');
    end
    if floor(vars.MaxIterations)~=vars.MaxIterations
        error('The maximal number of iterations has to be an integer.');
    end
    if sum(size(vars.Epsilon))>2
        error('The upper bound for cyclic step change has to be a number');
    end
    % Error checks
    if ( isinf(maxiter) && (isinf(epsilon)) )
        error('You can not set both the maximal Iterations and Epsilon to Inf');
    end
    % alpha (Check and form to a vector for future simplicity)
    assert(isvector(vars.alpha),'alpha has to be a number or a vector of length 4 at most.');
    assert(all(vars.alpha>=0) && all(vars.beta>=0),...
        'All entries of alpha and beta have to be nonnegative.');
    assert(isvector(vars.beta),...
        'beta has to be a number or a vector of length 3 at most.');
    % handle simplifications
    if (length(vars.alpha)==1) % use Axis TV
        alpha = [1,1,0,0]*vars.alpha;
    elseif (length(vars.alpha)<5) % set the first ones
        alpha = zeros(1,4);
        alpha(1:length(vars.alpha)) = vars.alpha;
    else
        error('the vector alpha is too long');
    end
        %
    % beta (Check and form to a vector for future simplicity)
    if (length(vars.beta)==1) % number
        beta = [1,1,1]*vars.beta;
    elseif (length(vars.beta)==2) %
        beta = [vars.beta(1),vars.beta(1),vars.beta(2)];
    elseif (length(vars.beta)==3)
        beta = vars.beta;
    else
        error('beta is a vector longer than 3 and hence too long.');
    end
    % N is the number of data points, k the dimension of the manifold, and they
    % are given as columns
    dimen = size(vars.f);
    manDim = dimen(1:length(vars.M.ItemSize));
    imgDim = dimen(length(manDim)+1:end);
    if all(size(imgDim) == [1 1])
       imgDim = [imgDim,1]; 
    end
    assert(all(manDim == vars.M.ItemSize),...
    ['Data Ite dimensions (',num2str(dimen(1:end-1))...
        ,' not consistent with manifold specification ',...
        num2str(vars.M.ItemSize),'.']);
    %
    % Handle optional input RegMask
    useRegMask = false;
    rM = ones(imgDim); %Standard: All are regularized, none is fixed
    if (numel(vars.RegMask)>0) % Mask given
        if any(size(vars.RegMask) ~= imgDim)
            warning(['Size of the regularization mask (',...
                num2str(size(vars.RegMask)),') does not fit the size of image f (',...
                num2str(imgDim),'). Hence it is ignored.']);
        else
            rM = vars.RegMask;
            useRegMask =true;
        end
    end
    %
    % Handle optional input UnknownMask
    useUnknownMask = false;
    uM = zeros(imgDim(1),imgDim(2)); %Standard: All are known
    if (numel(vars.UnknownMask)>0) % Mask given
        if any(size(vars.UnknownMask) ~= imgDim)
            warning(['Size of the unknown data mask (',...
                num2str(size(vars.UnknownMask)),') does not fit the size of image f (',...
                num2str(imgDim),'). Hence it is ignored.']);
        else
            uM = vars.UnknownMask;
            useUnknownMask =true;
        end
    end
    evalFct = vars.EvalFct;
end %end parser & checks
% record l2TV & eps?
recordFct = (nargout==2) && ( sum(size(evalFct)) > 0) && isa(evalFct,'function_handle');
if recordFct
    fctValSeq = [];
end
i=0;
lambdait=vars.lambda;
x = vars.f;
if (useUnknownMask) %kill unknown pixels completely
    x(repmat(permute(uM==1,[(length(imgDim)+(1:length(vars.M.ItemSize))),1:length(imgDim)]),[manDim,ones(1,length(imgDim))])) = NaN;
end
itD = [Inf, NaN];
debug('time',3,'StartTimer','Cyclic proximal point algorithm on manifolds');
while ( (any(isnan(itD(:))) || (max(itD(~isnan(itD)))>=epsilon)) && (i<maxiter))
    xold = x;
    % First step : proximal mappings of the distances - only on
    % known points
    Inds = ones(size(x));
    if useRegMask
        Inds = Inds&repmat(permute(rM==1,[length(imgDim)+(1:length(vars.M.ItemSize)),1:length(imgDim)]),...
            [manDim,ones(1,length(imgDim))]);
    end
    if useUnknownMask
        Inds = Inds&repmat(permute(~uM,...
            [(length(imgDim)+(1:length(vars.M.ItemSize))),1:length(imgDim)]),...
        [manDim,ones(1,length(imgDim))]);
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
    for j=0:1
        Nt = floor((imgDim(1)-j)/2);
        Mt = floor((imgDim(2)-j)/2);
        % reshape transposed into two columns
        % X
        if (alpha(1) > 0)
            % (a) permute the subsampled dimension to second (most inner
            % permute, not necessary above)
            % (b) reshape to manifoldx2xRestsize(
            % (c) permute 2 to third dim for easier proxes
            problem.f = permute(reshape(permute(x(:,:,j+1:2*Mt+j),[1,3,2]),[prod(manDim),2,imgDim(1)*Mt]),[1,3,2]);
            problem.f = reshape(problem.f,[manDim,imgDim(1)*Mt,2]); %expand manifold dims again
            problem.lambda = alpha(1)*lambdait;
            problem.w = [-1,1]';
            problem.RegMask = permute(reshape(permute(rM(:,j+1:2*Mt+j),[2,1]),[2,imgDim(1)*Mt]),[2,1]);
            x(:,:,j+1:2*Mt+j) = permute(reshape(permute(...
                reshape(vars.M.proxad(problem),[prod(manDim),imgDim(1)*Mt,2]),... %contract manifold dims again
                [1,3,2]),[prod(manDim),2*Mt,imgDim(1)]),[1,3,2]); %rearrange for x (invert a-c)
        end
        % Y
        if (alpha(2) > 0)
            % (a) reshape to manifoldx2xRest
            % (b) permute 2 to third dim for easier proxes
            problem.f = permute(reshape(x(:,j+1:2*Nt+j,:),[prod(manDim),2,Nt*imgDim(2)]),[1,3,2]);
            problem.f = reshape(problem.f,[manDim,Nt*imgDim(2),2]); %expand manifold dims again
            problem.lambda = alpha(2)*lambdait;
            problem.w = [-1,1]';
            problem.RegMask = permute(reshape(rM(j+1:2*Nt+j,:),[2,Nt*imgDim(2)]),[2,1]);
            x(:,j+1:2*Nt+j,:) = reshape(permute(...
                reshape(vars.M.proxad(problem),[prod(manDim),Nt*imgDim(2),2]),... %contract manifold dims again
                [1,3,2]),[prod(manDim),Nt*2,imgDim(2)]); %invert b-a
        end
        % both diagonals starting again from j even or odd 
        for l=0:1
            if (alpha(3+l)>0)
                Nt = floor((imgDim(1)-j)/2);
                problem.f = permute(reshape(... %internally take manifold as onedimensional for easier selection
                    ...% then this reshape shapes them back to multdim
                    cat(2,... %
                        reshape(x(:,(1+j):2:2*Nt+j,(1+l):imgDim(2)-1+l),prod(manDim),1,[]),... %all even rows,
                        reshape(x(:,(2+j):2:2*Nt+j,(2-l):imgDim(2)-l),prod(manDim),1,[])... %all odd rows
                    ),[manDim,2,Nt*(imgDim(2)-1)]),[(1:length(vars.M.ItemSize)),length(vars.M.ItemSize)+[2,1]]); %and permute diff-size (2) last
                problem.w = [-1,1]';
                problem.lambda = alpha(3+l)*lambdait;
                problem.RegMask = permute([ reshape(rM((1+j):2:2*Nt+j,(1+l):imgDim(2)-1+l),1,[]);... %all even rows,
                                   reshape(rM((2+j):2:2*Nt+j,(2-l):imgDim(2)-l),1,[]) ],[2,1]); %all odd rows
                res = reshape(vars.M.proxad(problem),[prod(manDim),Nt*(imgDim(2)-1),2]); %contract manifold dims again, keep diffs last
%                reshape result to old dimensions (i.e. esp. one manifold dim) and save it.
                x(:,(1+j):2:2*Nt+j,(1+l):imgDim(2)-1+l) = reshape(res(:,:,1),size(x(:,(1+j):2:2*Nt+j,(1+l):imgDim(2)-1+l,:)));
                x(:,(2+j):2:2*Nt+j,(2-l):imgDim(2)-l) = reshape(res(:,:,2),size(x(:,(2+j):2:2*Nt+j,(2-l):imgDim(2)-l,:)));
            end
        end
    end % end j, both splitting blocks
    %
    % Third: Second order differences
    for j=0:2
        Nt = floor(abs(imgDim(1)-j)/3);
        Mt = floor(abs(imgDim(2)-j)/3);
        % X
        if (beta(1) > 0)
            % (a) permute the subsampled dimension to second (most inner
            % permute, not necessary above)
            % (b) reshape to manifoldx2xRest
            % (c) permute 2 to third dim for easier proxes
            problem.f = permute(...
                reshape(...
                    permute(x(:,:,(j+1):(3*Mt+j)),[1,3,2]),...
                [prod(manDim),3,Mt*imgDim(1)]),...
            [1,3,2]);
            problem.f = reshape(problem.f,[manDim,Mt*imgDim(1),3]); %expand manifold Dims
            problem.lambda = beta(1)*lambdait;
            problem.w = [1,-2,1]';
            problem.RegMask = permute(reshape(permute(rM(:,j+1:3*Mt+j),[2,1]),[3,imgDim(1)*Mt]),[2,1]);
            x(:,:,j+1:3*Mt+j) = permute(reshape(permute(...
                    reshape(vars.M.proxad(problem),[prod(manDim),Mt*imgDim(1),3]),... %contract manifld dims
                [1,3,2]),[prod(manDim),3*Mt,imgDim(1)]),[1,3,2]);
        end
        % Y
        if (beta(2) > 0)
            % (a) reshape to manifoldx3xRest
            % (b) permute 3 last
            problem.f = permute(reshape(x(:,j+1:3*Nt+j,:),[prod(manDim),3,Nt*imgDim(2)]),[1,3,2]);
            problem.f = reshape(problem.f,[manDim,Nt*imgDim(2),3]); %expand manifold Dims
            problem.lambda = beta(2)*lambdait;
            problem.w = [1,-2,1]';
            problem.RegMask = permute(reshape(rM(j+1:3*Nt+j,:),[3,Nt*imgDim(2)]),[2,1]);
            x(:,j+1:3*Nt+j,:) = reshape(permute(...
                reshape(vars.M.proxad(problem),[prod(manDim),Nt*imgDim(2),3]),...%collapse manifold dims
                [1,3,2]),[prod(manDim),Nt*3,imgDim(2)]);
        end
    end
    %
    % Fourth: b_11
    if (beta(3) > 0)
        Nt = [2*floor(imgDim(1)/2),2*floor((imgDim(1)-1)/2)];
        Mt = [2*floor(imgDim(2)/2),2*floor((imgDim(2)-1)/2)];
        for l=0:1
            for j=0:1
                problem.f = cat(2,...
                    reshape(x(:,1+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j,:),prod(manDim),1,[]),...
                    reshape(x(:,2+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j,:),prod(manDim),1,[]),...
                    reshape(x(:,2+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j,:),prod(manDim),1,[]),...
                    reshape(x(:,1+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j,:),prod(manDim),1,[]));
                problem.f= permute(reshape(problem.f,[manDim,4,(Mt(j+1)/2)*(Nt(l+1)/2)]),[(1:length(vars.M.ItemSize)),length(vars.M.ItemSize)+[2,1]]); %expand manifold dims and put diffs last
                problem.w = [-1,1,-1,1]';
                problem.lambda = beta(3)*lambdait;
                problem.RegMask = permute([reshape(rM(1+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j),1,[]);...
                    reshape(rM(2+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j),1,[]);...
                    reshape(rM(2+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j),1,[]);...
                    reshape(rM(1+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j),1,[])],[2,1]); %put again number of terms per diff lasts
                tx = reshape(vars.M.proxad(problem),[prod(manDim),(Mt(j+1)/2)*(Nt(l+1)/2),4]); % collapse for rearranging
                x(:,1+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j) = reshape(tx(:,:,1),size( x(:,1+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j)));
                x(:,2+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j) = reshape(tx(:,:,2),size( x(:,2+l:2:Nt(l+1)+l, 1+j:2:Mt(j+1)+j)));
                x(:,2+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j) = reshape(tx(:,:,3),size( x(:,2+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j)));
                x(:,1+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j) = reshape(tx(:,:,4),size( x(:,1+l:2:Nt(l+1)+l, 2+j:2:Mt(j+1)+j)));
            end
        end
    end
    i = i + 1;
    lambdait = vars.lambda/(i);
    % revert (A) (roughly line 199)
    x = reshape(x,[manDim,imgDim]); %expand manifold dimensions from first dim
    itD = NaN*zeros(imgDim);
    Inds = ~permute(any(isnan(reshape(xold,[prod(manDim),imgDim])),1),[2:(length(imgDim)+1),1]);
    % only DIst nonNan-ones
    itD(Inds) = vars.M.dist(...
            reshape(x(repmat(permute(Inds,[length(imgDim)+1:length(imgDim)+length(manDim),1:length(imgDim)]),...
            [manDim,ones(1,length(imgDim))])),...
            [manDim,sum(Inds(:))]),...
            reshape(xold(repmat(permute(Inds,[length(imgDim)+1:length(imgDim)+length(manDim),1:length(imgDim)]),...
            [manDim,ones(1,length(imgDim))])),...
            [manDim,sum(Inds(:))])...
        );
    % Record ?
    if recordFct
           fctValSeq = [fctValSeq,evalFct(x,itD)]; %#ok<AGROW>
    end
    if mod(i,getDebugLevel('IterationStep'))==0
        debug('text',3,'text',...
        ['i=',num2str(i),' where lastdiff is ',num2str(max(itD(~isnan(itD)))),' and lambdait=',num2str(lambdait)]);    
    end
end
debug('time',3,'StopTimer','Cyclic proximal point algorithm on manifolds');
debug('text',2,'Text',...
    [num2str(i),' iterations, last difference:',...
    num2str(max(itD(~isnan(itD))))]);
end
