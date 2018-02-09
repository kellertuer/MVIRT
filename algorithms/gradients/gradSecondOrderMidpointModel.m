function eta = gradSecondOrderMidpointModel(varargin)
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
%   'p       : (p=1) compute TV with p-norm coupling in the dimensions of the
%             data, i.e. anisotropic TV for p=1 and isotropic for p=2
%  epsilon   : compute the gradient of the epsilon-relaxed TV
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

eta = gradData(vars.M,vars.f,vars.x,2) ...
    + vars.alpha.*gradTV(vars.M,vars.x,'p',vars.p,'Epsilon',vars.Epsilon) ...
    + vars.beta*gradTV2Midpoint(vars.M,vars.x,'p',vars.p,'Epsilon',vars.Epsilon);
end

