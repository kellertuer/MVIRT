function [x,recData] = subGradientDescent(varargin)
% subGradientDescent(M,x,F,gradF,stepSizeRule,stoppingCriterion)
% compute a sub gradient descent on F using x as inital data and subgradF
% as the (sub)gradient. If F returns a number, the classical subgradient
% descent is performed; if it returs a vector or matrix; a parallel
% subgradient descent is computed, where the (last) data dimensions of x
% are handled in parallel.
%
%   Ferreira, O.P.; Oliverira, P.R.: Subgradient Algorithm on Riemannian
%       manifolds
%       Journal of Optimization Theory and Applications 97.1, 93-104, 1998.
%
% INPUT
%    M                : a manifold M
%    x                : initial data
%    F                : the function F for the different gradients in
%                       parallel or just a value for standard gradient
%                       descent
%    subgradF         : a function @(x) returning the gradient of F at x
%    stepSize         : a function @(x,eta,iter,initial) returning the
%                       stepsize at x with direction eta corresponding to
%                       a certain rule (e.g. Amijo) based on iterate point
%                       x, descent direction eta, iteration iter and
%                       initial stepsize initial.
%   stoppingCriterion : a function @(x,xold,iter) determining whether to
%                       stop or not
% OPTIONAL
%   Debug             : ([]) a function @ (x,xold,iter) producing debug
%                          output
%   Record            : (@(x,iter)) a function returning a column vector,
%                       if there's a second return value and this function
%                       is given, data is recorded in an array and returned
% ---
% MVIRT | R. Bergmann | 2018-03-15

ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addRequired(ip,'F');
addRequired(ip,'subgradF');
addRequired(ip,'stepSize');
addRequired(ip,'stoppingCriterion');
addOptional(ip,'Debug',[]);
addOptional(ip,'Record',[]);
parse(ip, varargin{:});
vars = ip.Results;
if isempty(vars.Debug)
    debugF = @(x,xold,iter) [num2str(iter),' last change: ',...
        num2str(sum(vars.M.dist(x(vars.M.allDims{:},:),xold(vars.M.allDims{:},:)))/numel(x)*prod(vars.M.ItemSize)),'.'];
else
    debugF = vars.Debug;
end
record = ~isempty(vars.Record) && nargout > 1;
if record
    recData = [];
end
sX = size(vars.x);
mD = length(vars.M.ItemSize); %manifold item size
if length(sX)>mD
    dD = sX(mD+1:end); %data dimensions
else
    dD=1;
end
dL = length(dD); %numer of data dimensions
InitValue = vars.F(vars.x);
parallel = ~isscalar(InitValue);
if parallel
    pD = size(InitValue);
    pD = pD(pD>1); % omit last dimension if we are parallel on a vector
    pL = length(pD); %number of parallel dimensions
    % check that these are the last dimensions of dD; init corresponding
    % adressers
    if any(pD ~= dD( (dL-pL+1):end) )
        error(['The dimensions F provides (',num2str(pD),...
            ') does not fit the last data dimensions (',...
            num2str(dD( (dL-pL):end)),')']);
    end
    % remaining data adressor
    dataF = repelem({':'},dL-pL);
end
x = vars.x;
xold=x;
xOpt = x;
iter = 0;
s = 1;
while ~vars.stoppingCriterion(x,xold,s,iter)
    iter=iter+1;
    xold = x;
    descentDir = -vars.subgradF(x);
    s = vars.stepSize(x,descentDir,iter,s);
    x = vars.M.exp(x,s.*descentDir);
    if parallel
        update = vars.F(x)<vars.F(xOpt);
        xOpt(vars.M.allDims{:},dataF{:},update) = x(vars.M.allDims{:},dataF{:},update);
    elseif vars.F(x)<vars.F(xOpt)
        xOpt = x;
    end
    if record
        recData = cat(2,recData,vars.Record(x,xold,iter));
    end
    if ~isempty(vars.Debug)
        debugF(x,xold,iter);
    end
end
end