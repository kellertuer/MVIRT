function t = SecondOrderMidpointModel(varargin)
%  gradientSecondOrderLog(M,f,alpha,beta) the l^2-TV-TV2 mid point model
%
% INPUT
%    M : a manifold
%    f : given (original) data
%    x : value to take the gradient at
%    alpha : weight of TV
%    beta : weight of TV2
%
% OPTIONAL
%   'p'      : (p=1) compute TV with p-norm coupling in the dimensions of
%              the data, i.e. anisotropic TV for p=1 and isotropic for p=2
%  'epsilon' : compute the gradient of the epsilon-relaxed TV and TV2
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

