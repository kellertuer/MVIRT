function eta = gradTV(varargin)
% gradTV(M,x) compute the gradient of the manifold total variation
%
% INPUT
%   M      :  a manifold
%   x      : data (size [manDims,dataDims])
%   
% OPTIONAL
%   'p'       : (p=1) compute TV with p-norm coupling in the dimensions of
%               the data, i.e. anisotropic TV for p=1 and isotropic for p=2
%  'epsilon'  : compute the gradient of the epsilon-relaxed TV
%  'weights'  : (ones(dataDims)) exclude certain data points from all
%               gradient terms
%
% OUTPUT
%   eta : the gradient
% ---
% MVIRT, R. Bergmann, 2017-12-08
ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addOptional(ip,'p',1);
addOptional(ip,'Epsilon',0);
addOptional(ip,'Weights',[]);
parse(ip, varargin{:});
vars = ip.Results;
sX = size(vars.x);
mD = length(vars.M.ItemSize);
dataDims = sX( (mD+1):end );
n = length(dataDims);
if isempty(vars.Weights)
    weights = ones([dataDims,1]);
else
    weights = vars.Weights;
end
eta = zeros(size(vars.x));
if vars.p>1
    prefactors1 = TV(vars.M,vars.x,'Sum',false,'p',vars.p,'Epsilon',vars.Epsilon);
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
        prefactors1 = sqrt(vars.M.dist(center,forward).^2+vars.Epsilon^2);
        prefactors2 = sqrt(vars.M.dist(center,backward).^2+vars.Epsilon^2);
    else
        prefactors2 = prefactors1(preFill{:},[1 1:(dataDims(i)-1)],postFill{:});
    end
    w1 = weights.*forwardweights.*double(prefactors1~=0)./(prefactors1+double(prefactors1==0));
    w1(preFill{:},dataDims(i),postFill{:}) = 0;
    w2 = weights.*backwardweights.*double(prefactors2~=0)./(prefactors2+double(prefactors2==0));
    w2(preFill{:},1,postFill{:}) = 0;
    eta = eta   - vars.M.log(center,forward).*shiftdim(w1,-mD)...
                - vars.M.log(center,backward).*shiftdim(w2,-mD);
end

