function [x,fctValSeq] = cppa_ad_1D(varargin)
% cppa_ad_1D(M,f,alpha,beta,lambda)
%   Compute the cyclic proximal point algorithm for first and
%   second order absolute differences, weighted by alpha and
%   beta on a manifold
% cppa_ad_1D(problem)
%     performs the same, but the parameters and optional
%     parameters are given in the struct problem.
%
% INPUT
%   M      : a manifold.
%   f      : a vector from M^d of M-valued data, where size is [M.ItemSize,d]
%   alpha  : weight for first order absolute difference terms
%   beta   : weight for second order absolute difference terms
%   lambda : weight for the proximal mappings
%
%   OUTPUT:
%      x   : result of the cyclic proximal point algorithm (again
%      [M.ItemSize,d])
%
%   OPTIONAL PARAMETERS
%      'MaxIterations' : (400) Maximal number of Iterations
%      'Epsilon'       : (10^{-6}) Lower bound for the max change in one cycle
%
%           one of the parameters MaxIterations/Epsilon can be dactivated
%           by specifying them as Inf but not both.
%
% 'RegMask'     : ([]) Specify a binary mask for values affected by
%                      (1) and not affected by the difference proximal
%                      mappings (0) (i.e. fixed input values).
% 'UnknownMask' : ([]) Specify a mask, which values of f are unknown
%                      (1) and initialized in the cycles, when possible.
%
%   'EvalFct'   : ([]) specify a function to be evaluated in every
%                          iteration and accepts two parameter, the x_k
%                          and the distance array from x_k to x_{k-1} measured on
%                          the manifold
%
%  See also:
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann,
%     A second order non-smooth variational model for restoring
%     manifold-valued images,
%     SIAM Journal on Scientific Computing, 38, (1), A567?A597, 2016.
%
%    manifold, manifold.proxad, cppa_ad_2D, cppa_ad_nD
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-10-22 | 2016-02-22
% see LICENSE.txt

%
% Changelog
% 2015-01-20 Revised data arrangement to fit matlab conventions
% 2014-12-01 Changed the proxAD calls to use the struct
% 2014-11-29 Extracted Algorithms from the manifold class to be more
%            independent and added structs
if (length(varargin)==1 && isstruct(varargin{1}))
    % short hand call to surpass parser
    vars = varargin{1};
    assert(all(isfield(vars,{'M','f','alpha','beta','lambda'})),'There are required input fields missing'); % at least required
    assert(isa(vars.M,'manifold'),'Input M is not a manifold');
    % Fill optionals
    dimen = size(vars.f);
    N = dimen(end); %length of data
    if (N==1) && (length(dimen)==2) % zero-dimensional-manifold i.e. values given as column -> prime
        vars.f = vars.f';
        N = dimen(1);
    end
    useRegMask = isfield(vars,'RegMask');
    if ~useRegMask
        rM = ones(1,N); %Standard: All are regularized, none is fixed
    else
        rM = vars.RegMask;
    end
    useUnknownMask = isfield(vars,'UnknownMask');
    if ~useUnknownMask
        uM = zeros(1,N);
    else
        uM = vars.UnknownMask;
    end
    if ~isfield(vars,'MaxIterations')
        vars.MaxIterations = 400;
    end
    maxiter = vars.MaxIterations;
    if ~isfield(vars,'Epsilon')
        vars.Epsilon = 10^(-6);
    end
    if ~isfield(vars,'EvalFct')
        evalFct = [];
    else
        evalFct = vars.EvalFct;
    end
    epsilon = vars.Epsilon;
else % Parse Variables & check
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
    % Validate
    %
    % Validate optional stuff
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
    assert( ~ (isinf(maxiter) && (epsilon==0)),'You can not set both the maximal Iterations and Epsilon to Inf');
    assert(epsilon>= 0, 'Epsilon has to be positive');
    assert( isinf(maxiter) || (round(maxiter)==maxiter) && (maxiter > 0), 'MaxIterations has to be a positive integer.');
    assert(vars.lambda>0,'Lambda has to be a positive real number');
    
    % N is the number of data points, k the dimension of the manifold, and they
    % are given as columns
    dimen = size(vars.f);
    N = dimen(end); %length of data
    if (N==1) && (length(dimen)==2) % zero-dimensional-manifold i.e. values given as column -> prime
        vars.f = vars.f';
        dimen  = size(vars.f);
    end
    assert(all(dimen(1:end-1) == vars.M.ItemSize),...
        ['Data Ite dimensions (',num2str(dimen(1:end-1))...
        ,' not consistent with manifold specification ',...
        num2str(vars.M.ItemSize),'.']);
    %
    % Optional parameter: RegMask
    useRegMask = false;
    rM = ones(1,N); %Standard: All are regularized, none is fixed
    if (numel(vars.RegMask)>0) % Mask given
        if ( ~isvector(vars.RegMask) || (length(vars.RegMask) ~= N) )
            warning(['Length of the regularization mask (',...
                num2str(length(vars.RegMask)),') does not fit the size of signal f (',...
                num2str(N),'). Hence it is ignored.']);
        else
            if (isrow(vars.RegMask))
                rM = vars.RegMask;
            else
                rM = vars.RegMask';
            end
            useRegMask =true;
        end
    end
    %
    % Optional parameter: Unknownmask
    useUnknownMask=false;
    uM = zeros(1,N); %Standard: All are known
    if (numel(vars.UnknownMask)>0) % Mask given
        if ( ~isvector(vars.UnknownMask) || (length(vars.UnknownMask) ~= N) )
            warning(['Length of the unknown data mask (',...
                num2str(length(vars.UnknownMask)),') does not fit the size of signal f (',...
                num2str(N),'). Hence it is ignored.']);
        else
            if (isrow(vars.UnknownMask))
                uM = vars.UnknownMask;
            else
                uM = vars.UnknownMask';
            end
            useUnknownMask =true;
        end
    end
    
end
xold = zeros(size(vars.f));
x = vars.f;
if useRegMask&&useUnknownMask
    fixedAndUnknown = uM&(~rM);
    assert( ~any(fixedAndUnknown(:)), 'There exists at least one pixel that is unknown and fixed by the two masks, which is not allowed.');
end
if (useUnknownMask) %kill unknown pixels completely
    x(repmat(uM==1,[k,1])) = NaN;
end
recordFct = (nargout==2) && ( sum(size(evalFct)) > 0) && isa(evalFct,'function_handle');
if recordFct
    fctValSeq = [];
end
i=0;
lambdait=vars.lambda;
itD = [Inf,NaN];
debug('time',3,'StartTimer','Cyclic proximal point algorithm on manifolds');
while ( (any(isnan(itD(:))) || (max(itD(~isnan(itD)))>epsilon)) && (i<maxiter))
    xold = x;
    % First step : proximal mappings of the distances
    proxD = vars.M.proxDist(x,vars.f,lambdait);
    if (useUnknownMask)
        x(repmat(~uM,[3,1])) = proxD(repmat(~uM,[3,1]));
    else
        x = proxD;
    end
    %
    x = reshape(x,[prod(vars.M.ItemSize),N]); %linearize manifold dimensions into first dim
    if vars.alpha > 0
        for j=0:1 %even and odd splitting
            Nt = floor((N-j)/2);
            % X
            % reshape transposed into two columns
            % reshape permute to get ordering right
            problem.f = permute(reshape(x(:,j+1:2*Nt+j),[prod(vars.M.ItemSize),2,Nt]),[1,3,2]); %only for input back to M
            problem.f = reshape(problem.f, [vars.M.ItemSize,Nt,2]); %expand after permute
            problem.lambda = vars.alpha*lambdait;
            problem.w = [-1,1]';
            problem.RegMask = permute(reshape(rM(j+1:2*Nt+j),[2,Nt]),[2,1]);
            % collapse manifold dims, permute back, expand.
            x(:,j+1:2*Nt+j) = reshape(permute(vars.M.proxad(problem),[(1:length(vars.M.ItemSize)),length(vars.M.ItemSize)+[2,1]]),[prod(vars.M.ItemSize),2*Nt]);
        end %end split even odd
    end
    if vars.beta > 0
        for j=0:2
            Nt = floor((N-j)/3);
            problem.f = permute(reshape(x(:,j+1:3*Nt+j),[prod(vars.M.ItemSize),3,Nt]),[1,3,2]);
            problem.f = reshape(problem.f, [vars.M.ItemSize,Nt,3]); %expand after permute
            problem.lambda = vars.beta*lambdait;
            problem.w = [1,-2,1]';
            problem.RegMask = permute(reshape(rM(j+1:3*Nt+j),[3,Nt]),[2,1]);
            x(:,j+1:3*Nt+j) = reshape(permute(vars.M.proxad(problem),[(1:length(vars.M.ItemSize)),length(vars.M.ItemSize)+[2,1]]), [prod(vars.M.ItemSize),3*Nt]);
        end
    end
    x = reshape(x,[vars.M.ItemSize,N]); %back to multiple manifold dimensions
    i = i + 1;
    lambdait = vars.lambda/(i);
    itD = vars.M.dist(x,xold);
    if recordFct
        fctValSeq = [fctValSeq,evalFct(x,itD)]; %#ok<AGROW>
    end
    if mod(i,getDebugLevel('IterationStep'))==0
        debug('text',3,'text',...
            ['i=',num2str(i),' where lastdiff is ',num2str(max(itD(~isnan(itD)))...
            ),' and lambdait=',num2str(lambdait)]);
    end
end
debug('time',3,'StopTimer','Cyclic proximal point algorithm on manifolds');
debug('text',2,'Text',...
    [num2str(i),' iterations, last difference:',...
    num2str(max(vars.M.dist(x,xold)))]);
end
