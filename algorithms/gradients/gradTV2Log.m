function eta = gradTV2Log(varargin)
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
    prefactors = TV2(vars.M,vars.x,'Sum',false,'p',vars.p,'Epsilon',vars.Epsilon);
end
for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    center = vars.x;
    forward = vars.x(vars.M.allDims{:},preFill{:},[2:(dataDims(i)) (dataDims(i))],postFill{:});
    backward = vars.x(vars.M.allDims{:},preFill{:},[1 1:(dataDims(i)-1)],postFill{:});
    forwardweights = weights(preFill{:},[2:(dataDims(i)) (dataDims(i))],postFill{:});
    backwardweights = weights(preFill{:},[1 1:(dataDims(i)-1)],postFill{:});
    if vars.p==1
        logsum = vars.M.log(center,forward) + vars.M.log(center,backward);
        prefactors = sqrt(vars.M.dot(center,logsum,logsum)+vars.Epsilon^2);
    end
    inner = -2*vars.M.log(center,forward) -2*vars.M.log(center,backward);
    denom = permute( sqrt(vars.M.dot(center,inner,inner)+vars.Epsilon.^2), [(n+1):(n+length(vars.M.ItemSize)),1:n]);
    outer = vars.M.AdjDyLog(center,forward,inner) ... % appears as x in a term above
        + vars.M.AdjDxLog(center,forward,inner) ... % as y (first)
        + vars.M.AdjDxLog(center,backward,inner) ... %as y
        + vars.M.AdjDyLog(center,backward,inner);
    w = permute(...
            weights.*backwardweights.*forwardweights.*(prefactors~=0)./(vars.p*prefactors+(prefactors==0)),...
        [(n+1):(n+length(vars.M.ItemSize)),1:n]);
    eta = eta + w.*outer.*(denom~=0)./(denom+(denom==0)); %divide by outer prefactor and avoid division by zero
end
end