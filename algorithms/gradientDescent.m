function [x,recData] = gradientDescent(varargin)
% gradientDescent(M,x,gradF,stepSizeRule,stoppingCriterion)
%    compute a gradient descent based on
%
% INPUT
%    M                : a manifold M
%    x                : initial data
%    gradF            : a function @(x) returning the gradient of F at x
%    stepSize         : a function @(x,eta,iter,initial) returning the
%                       stepsize at x with direction eta corresponding to a
%                       certain rule (e.g. Amijo) based on iterate point x,
%                       descent direction eta, iteration iter and initial
%                       stepsize initial.
%   stoppingCriterion : a function @(x,xold,iter) determining whether to
%                       stop or not
%
% OPTIONAL
%   Debug             : ([]) a function @ (x,xold,iter) producing debug
%                          output
%   Record            : (@(x,iter)) a function returning a column vector,
%                       if there's a second return value and this function
%                       is given, data is recorded in an array and returned
% ---
% MVIRT | R. Bergmann | 2018-01-22

ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addRequired(ip,'gradF');
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
x = vars.x;
xold=x;
iter = 0;
s = 1;
while ~vars.stoppingCriterion(x,xold,s,iter)
    iter=iter+1;
    xold = x;
    descentDir = -vars.gradF(x);
    s = vars.stepSize(x,descentDir,iter,s);
    x = vars.M.exp(x,s.*descentDir);
    if record
        recData = cat(2,recData,vars.Record(x,xold,iter));
    end
    if ~isempty(vars.Debug)
        debugF(x,xold,iter);
    end
end
end