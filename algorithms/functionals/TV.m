function V = TV(varargin)
% TV(M,x) - compute all TV terms (and sum them
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
weights = vars.weights;
if isscalar(weights)
    weights = ones(1,n).*weights; %to vector
end
if isvector(weights) % is a vector
    weights = diag(weights); %to matrix
end
tv = zeros([dataDims,1]);
for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    center = vars.x;
    forward = vars.x(vars.M.allDims{:},preFill{:},[2:(dataDims(i)) (dataDims(i))],postFill{:});
    d = vars.M.dist(center,forward);
    if vars.Epsilon>0 && vars.p==1 %for anisotropic relax each summand
        d = sqrt(d.^2+vars.Epsilon^2);
    end
    tv = tv + weights(i,i)*d.^vars.p;
    for j=0:1 %even odd for easier access
        for i2=i+1:n % diagonals
            interfill = repelem({':'},i2-i-1);
            post2Fill = repelem({':'},n-i2);
            Nt2 = floor((dataDims(i2)-j)/2);
            if dataDims(i) >= 2+j && dataDims(i2) >= 2+j
                if weights(i,i2) > 0
                    % extract even/odd
                    subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:});
                    % +1,+1 diagonal
                    subXd1 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:});
                    d = vars.M.dist(subXo,subXd1);
                    if vars.Epsilon>0 && vars.p==1 %for anisotropic relax each summand
                        d = sqrt(d.^2+vars.Epsilon^2);
                    end
                    tv = tv + weights(i,i2)*d.^vars.p;
                end
            end
            if weights(i2,i) > 0
                subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:});
                subXd1 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:});
                d = vars.M.dist(subXo,subXd1);
                if vars.Epsilon>0 && vars.p==1 %for anisotropic relax each summand
                    d = sqrt(d.^2+vars.Epsilon^2);
                end
                tv = tv + weights(i,i2)*d.^vars.p;
            end
        end
    end
end
if vars.Epsilon>0 && vars.p>1 % for isotropic relax within coupling
    tv = tv+epsilon^2;
end
    tv = tv.^(1/vars.p);
if vars.Sum
    V = sum(tv(:));
else
    V = tv;
end