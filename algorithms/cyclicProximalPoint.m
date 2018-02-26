function [y,recData] = cyclicProximalPoint(varargin)
% CyclicProximalPoint(M, f, proximalMaps, stoppingCriterion)
%    Compute the cyclic proximal point algorithm starting with data f,
%    cyclic through the proximal points
%  Input
%   M                 : a manifold
%   x                 : given (initial) data
%   proximalMaps      : cell arrays of functionals computing the proximal maps,
%                         i.e. can be called as prox = @(x,lambda)
%   stoppingCriterion : a functional @(x,xold,lambda,iter) as a stopping
%                       criterion based on x, the last iterate xold,
%                       lambda,
%                       that indicates, then to stop, i.e., is 1 for stops.
%
% OPTIONAL
%   lambdaIterate     : (@ (iter)) a functional @(lambda,iter) computing the iterates value of lambda
%   Order       : ('') while the standard order is ascending in proximalMaps,
%                'random' takes a random proximal map in each iteration
%                   (note that then one iteration is only _one_ proximal map
%                   evaluated)
%                 'permute' : keeps a cyclic manner but permutes the order
%                             for each cycle anew.
%   Debug             : ([]) a function @ (x,xold,iter) producing debug
%                          output
%   Record            : (@(x,iter)) a function returning a column vector,
%                       if there's a second return value and this function
%                       is given, data is recorded in an array and returned
%
% OUTPUT
%   y       : result of the CPP algorithm
%   recData : recorded data array of the by the record-function
% ---
% R. Bergmann | MVIRT | 2018-01-22

if (length(varargin)==1 && isstruct(varargin{1})) % as struct
    % short hand call to surpass parser
    vars = varargin{1};
    assert(isfield(vars,{'M'}),'A manifold is not given');
    assert(isa(vars.M,'manifold'),'Input M is not a manifold');
    assert(isfield(vars,{'x'}),'no data x given');
    assert(isfield(vars,{'proximalMaps'}),'The cell array of proxial map functionals is missing');
    assert(isfield(vars,{'stoppingCriterion'}),'a stopping criterion function is missing');
    if ~isfield(vars,{'Debug'})
        vars.Debug = [];
    end
    if ~isfield(vars,{'Order'})
        vars.RandomOrder = '';
    end
    if ~isfield(vars,{'lambdaIterate'})
        if isfield(vars,{'lambda'}) % fallback for old cases
            vars.lambdaIterate = @(iter) vars.lambda/iter;
        else
            error('Please specify either a lambdaIterate function or an initial (vars.)lambda');
        end
    end
    if ~isfield(vars,{'Record'})
        vars.Record = [];
    end
else % parse Input
    ip = inputParser;
    addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
    addRequired(ip,'x');
    addRequired(ip,'proximalMaps');
    addRequired(ip,'stoppingCriterion');
    addOptional(ip,'Order','acsending');
    addOptional(ip,'Debug',[]);
    addOptional(ip,'Record',[]);
    addOptional(ip,'lambda',[]);
    addOptional(ip,'lambdaIterate',[]);
    parse(ip, varargin{:});
    vars = ip.Results;
end
record = ~isempty(vars.Record) && nargout > 1;
if record
    recData = [];
end
if ~isempty(vars.lambda) % priority over functional
    vars.lambdaIterate = @(iter) vars.lambda/iter;
end
if isempty(vars.lambdaIterate)
            error('Please specify either a lambdaIterate function or an initial (vars.)lambda');
end
x = vars.x;
xold=x;
iter = 0;
lambdak = vars.lambdaIterate(1);
numProx = length(vars.proximalMaps);
if strcmp(vars.Order,'random')
    numProx = 1; %choose randomly every prox call
elseif strcmp(vars.Order,'permute')
    order = randperm(numProx);
else
    order = 1:numProx;
end
 while ~vars.stoppingCriterion(x,xold,lambdak,iter)
    if strcmp(vars.Order,'random')
        order = randi(length(vars.proximalMaps)); %choose randomly
    elseif strcmp(vars.Order,'permute')
        order = randperm(numProx);
    end
    iter=iter+1;
    xold = x;
    lambdak = vars.lambdaIterate(iter);
    for i=1:numProx
        x = vars.proximalMaps{order(i)}(x,lambdak);
    end
    if ~isempty(vars.Debug)
        vars.Debug(x,xold,iter);
    end
    if record
        recData = cat(2,recData,vars.Record(x,xold,iter));
    end
end
y = x;
end

