function V = TV2_Midpoint(varargin)
% TV(M,x) - compute all TV2 terms (and sum them
%
% INPUT
%   M      :  a manifold
%   x      : data (size [manDims,dataDims])
%   
% OPTIONAL
%   'p       : (p=1) compute TV with p-norm coupling in the dimensions of the
%             data, i.e. anisotropic TV for p=1 and isotropic for p=2
%  epsilon   : compute the gradient of the epsilon-relaxed TV
%  weights   : (ones(dataDims) exclude certain data points from all gradient terms
%  Sum : (true) return a value (true) or a matrix of TV terms (false)
% ---
% MVIRT, R. Bergmann, 2017-12-08
ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addOptional(ip,'p',1);
addOptional(ip,'Epsilon',0);
addOptional(ip,'Weights',[]);
addOptional(ip,'Sum',true);
parse(ip, varargin{:});
vars = ip.Results;
sX = size(vars.x);
dataDims = sX( (length(vars.M.ItemSize)+1):end );
n = length(dataDims);
if isempty(vars.Weights)
    weights = ones([dataDims,1]);
else
    weights = vars.Weights;
end
tv2 = zeros([dataDims,1]);
for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    center = vars.x;
    forward = vars.x(vars.M.allDims{:},preFill{:},[1 3:(dataDims(i)) (dataDims(i))],postFill{:});
    backward = vars.x(vars.M.allDims{:},preFill{:},[1 1:(dataDims(i)-2) (dataDims(i))],postFill{:});
    d = vars.M.dist(center,vars.M.geopoint(forward,backward,0.5));
    if vars.Epsilon>0 && vars.p==1
        d = sqrt(d.^2+vars.Epsilon^2);
    end
    tv2 = tv2 + d.^vars.p;
end
if vars.Epsilon>0 && vars.p>1
    tv2 = tv2+epsilon^2;
end
    tv2 = weights.*(tv2.^(1/vars.p));
if vars.Sum
    V = sum(tv2(:));
else
    V = tv2;
end