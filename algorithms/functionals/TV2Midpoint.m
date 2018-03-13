function d = TV2Midpoint(varargin)
% proxTV(M,x,lambda) compute the second order TV mid point model value
%
% INPUT
%     M       : a manifold
%     x       : data (manifold-valued)
%
% OPTIONAL
%   'p'      : (p=1) compute TV2 with p-norm coupling in the dimensions of
%              the data, i.e. anisotropic TV2 for p=1 and isotropic for p=2
%  'epsilon' : compute the gradient of the epsilon-relaxed TV2
%  'weights' : (ones(dataDims) exclude certain data points from all
%                   gradient terms
%  'Sum'     : (true) return a value (true) or a matrix of TV2 terms (false)
%
% OUTPUT
%   d          : TV_2(x)
% ---
% MVIRT | R. Bergmann | 2017-12-11
ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addParameter(ip,'Weights',[]);
addParameter(ip,'Epsilon',0);
addParameter(ip,'Sum',true);

parse(ip, varargin{:});
vars = ip.Results;
sX = size(vars.x);
dataDims = sX( (length(vars.M.ItemSize)+1):end );
n = length(dataDims);
if isempty(vars.Weights)
    vars.Weights = eye(n);
end
if isscalar(vars.Weights) %lambda is a number
    vars.Weights = diag(ones(1,n).*vars.Weights);
elseif isvector(vars.Weights) %array
    vars.Weights = diag(vars.Weights); %diagonal matrix -> no mixed
end
d = zeros([dataDims,1]);

for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    for j=0:2
        % all starting with mod j+1 indices
        Nt = floor(abs(dataDims(i)-j)/3);
        if dataDims(i) >= 3+j && vars.Weights(i,i) > 0
            subX1 = vars.x(vars.M.allDims{:},preFill{:}, 1+j:3:3*Nt+j, postFill{:}); % 1
            subX2 = vars.x(vars.M.allDims{:},preFill{:}, 2+j:3:3*Nt+j, postFill{:}); % 2
            subX3 = vars.x(vars.M.allDims{:},preFill{:}, 3+j:3:3*Nt+j, postFill{:}); % 3
            s = vars.Weights(i,i)*vars.M.dist(vars.M.midPoint(subX1,subX3),subX2);
            if vars.Epsilon > 0
                s = s.^2;
            end
            d(preFill{:}, 2+j:3:3*Nt+j, postFill{:}) = ...
                d(preFill{:}, 2+j:3:3*Nt+j, postFill{:}) + ...
                s;
        end
    end
    % Mixed Differences
    for i2=i+1:n
        interfill = repelem({':'},i2-i-1);
        post2Fill = repelem({':'},n-i2);
        for j=0:1 % even and odd
            for j2=0:1
                % all starting with odd (even) indices in i and i2
                Nt = floor((dataDims(i)-j)/2);
                Nt2 = floor((dataDims(i2)-j2)/2);
                if dataDims(i) >= 2+j && dataDims(i2) >= 2+j2 && vars.Weights(i,i2) > 0
                    % extract all 4 elements of each 2x2 matrix
                    subX1 = vars.x(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:});
                    subX2 = vars.x(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:});
                    subX3 = vars.x(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j2:2:2*Nt2+j2,post2Fill{:});
                    subX4 = vars.x(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j2:2:2*Nt2+j2,post2Fill{:});
                    s = vars.Weights(i,i2)*vars.M.dist(...
                        vars.M.midPoint(subX1,subX3),...
                        vars.M.midPoint(subX2,subX4));
                    if vars.Epsilon > 0
                        s = s.^2;
                    end
                    d(preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:}) = ...
                        d(preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:})...
                        + s;
                end
            end
        end
    end
end
if vars.Epsilon > 0
    d = sqrt(d+vars.Epsilon^2);
end
if vars.Sum
    d = sum(d(:));
end
end