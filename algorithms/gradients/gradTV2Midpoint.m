function eta = gradTV2Midpoint(varargin)
% gradTV2Log ? compute the second order Log model ||Log_y x + Log_y z ||
% for all 3 successive terms in x.
%
% INPUT
% M : Manifold
% x : data to evaluate the gradient at
%
% OPTIONAL
%   'p       : (p=1) compute TV with p-norm coupling in the dimensions of the
%             data, i.e. anisotropic TV for p=1 and isotropic for p=2
%  epsilon   : compute the gradient of the epsilon-relaxed TV
%  weights   : (ones(dataDims) exclude certain data points from all gradient terms
%
% OUTPUT
%   eta : the gradient
ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addOptional(ip,'p',1);
addOptional(ip,'Epsilon',0);
addOptional(ip,'Weights',[]);
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
eta = zeros(size(vars.x));
if vars.p>1
    prefactors = TV2Midpoint(vars.M,vars.x,'Sum',false,'p',vars.p,'Epsilon',vars.Epsilon);
else
    prefactors = ones([dataDims,1]);
end
for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    center = vars.x;
    forward = vars.x(vars.M.allDims{:},preFill{:},[1 3:(dataDims(i)) (dataDims(i))],postFill{:});
    backward = vars.x(vars.M.allDims{:},preFill{:},[1 1:(dataDims(i)-2) (dataDims(i))],postFill{:});
    forwardweights = weights(preFill{:},[1 3:(dataDims(i)) (dataDims(i))],postFill{:});
    backwardweights = weights(preFill{:},[1 1:(dataDims(i)-2) (dataDims(i))],postFill{:});
    mid = vars.M.midPoint(backward,forward);
    d = permute(sqrt(vars.M.dist(mid,center).^2+vars.Epsilon^2),...
        [(n+1):(n+length(vars.M.ItemSize)),1:n]);
    inner = -vars.M.log(mid,center).*(d>0)./(d+(d==0));
    outer = vars.M.AdjDxGeo(backward,forward,0.5,inner) ... % center appears as x in a term above
        - vars.M.log(center,mid).*(d>0)./(d+(d==0)) ... % center appears as y (first)
        + vars.M.AdjDxGeo(forward,backward,0.5,inner); % center appears as z
    w = permute(...
            weights.*backwardweights.*forwardweights.*(prefactors~=0)./(vars.p*prefactors+(prefactors==0)),...
        [(n+1):(n+length(vars.M.ItemSize)),1:n]);
    eta = eta + w.*outer; %divide by outer prefactor and avoid division by zero
end
end