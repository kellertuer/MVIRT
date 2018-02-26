function t = SecondOrderMidpointModel(varargin)
%  gradientSecondOrderLog(M,f,alpha,beta) --- compute the gradient of
% d(x,f)^2 + alpha*TV(x) + beta*TV2(x), where TV(x) is d(x_i,x_i+1) and
% the second order terms are given by ||log_x_i x_i-1 + log_x_i x_i+1 ||
%
% INPUT
%    M : a manifold
%    f : given (original) data
%    x : value to take the gradient at
%    alpha : weight of TV
%    beta : weight of TV2
%
% OPTIONAL
% Epsilon: Relaxation
% ---
% MVIRT - R. Bergmann, 2017-12-07
ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'f');
addRequired(ip,'x');
addRequired(ip,'alpha');
addRequired(ip,'beta');
addOptional(ip,'p',1);
addOptional(ip,'Epsilon',0);
parse(ip, varargin{:});
vars = ip.Results;

data = 1/2*vars.M.dist(vars.f,vars.x).^2;
t = sum(data(:)) ...
    + TV(vars.M,vars.x,'Epsilon',vars.Epsilon,'Weights',vars.alpha) ...
    + TV2Midpoint(vars.M,vars.x,'Epsilon',vars.Epsilon,'Weights',vars.beta);
end

