function y = proxTV2(varargin)
% proxTV(M,x,lambda) compute all proxima of second order finite differences
% with a splitting modulo 3 for the mid point approach within the TV2
% term and a 2-2 splitting in the mixed terms in a cyclic manner
%
% INPUT
%     M       : a manifold
%     x       : data (manifold-valued)
%    lambda   : parameter of the prox
%
% OPTIONAL
% 'FixedMask' : binary mask the size of data items in x to indicate fixed
%   data items
%   'SecDiffProx' : (@(x1,x2,x3,lambda)
%                     proxAbsoluteSecondOrderDifference(M,x1,x2,x3,lambda))
%                      specify a prox for (mod 3) TV2 term proxes
%   'SecDiffMixProx' : (@(x1,x2,x3,x4,lambda)
%             proxAbsoluteSecondOrderMixedDifference(M,x1,x2,x3,x4,lambda))
%                      specify a prox for (mod 2,2) TV2 mixed term proxes
%
% OUTPUT
%   y          : result of applying the proximal maps in a cyclic manner.
% ---
% Manifold-valued Image Restoration Toolbox 1.2 | R. Bergmann | 2017-12-11
ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addRequired(ip,'lambda');
addParameter(ip,'FixedMask',[]);
addParameter(ip,'SecDiffProx', []);
addParameter(ip,'SecDiffMixProx', []);
parse(ip, varargin{:});
vars = ip.Results;
if isempty(vars.SecDiffProx)
    vars.SecDiffProx = @(x1,x2,x3,lambda) ...
        proxAbsoluteSecondOrderDifference(vars.M,x1,x2,x3,lambda);
end
if isempty(vars.SecDiffMixProx)
    vars.SecDiffMixProx = @(x1,x2,x3,x4,lambda) ...
        proxAbsoluteSecondOrderMixedDifference(vars.M,x1,x2,x3,,x4,lambda);
end
sX = size(vars.x);
dataDims = sX( (length(vars.M.ItemSize)+1):end );
n = length(dataDims);
if isscalar(vars.lambda) %lambda is a number
    vars.lambda = diag(ones(1,n).*vars.lambda);
elseif isvector(vars.lambda) %array
    vars.lambda = diag(vars.lambda); %diagonal matrix -> no mixed
end
y = vars.x;

for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    for j=0:2
        % all starting with mod j+1 indices
        Nt = floor(abs(dataDims(i)-j)/3);
        if dataDims(i) >= 3+j && vars.lambda(i,i) > 0
            subX1 = y(vars.M.allDims{:},preFill{:}, 1+j:3:3*Nt+j, postFill{:}); % 1
            subX2 = y(vars.M.allDims{:},preFill{:}, 2+j:3:3*Nt+j, postFill{:}); % 2
            subX3 = y(vars.M.allDims{:},preFill{:}, 3+j:3:3*Nt+j, postFill{:}); % 3
            [pX1,pX2,pX3] = vars.SecDiffProx(subX1,subX2,subX3,vars.lambda(i,i));
            % write back
            y(vars.M.allDims{:},preFill{:}, 1+j:3:3*Nt+j, postFill{:}) = pX1;
            y(vars.M.allDims{:},preFill{:}, 2+j:3:3*Nt+j, postFill{:}) = pX2;
            y(vars.M.allDims{:},preFill{:}, 3+j:3:3*Nt+j, postFill{:}) = pX3;
            if ~isempty(vars.FixedMask)
                y(vars.M.allDims{:},vars.FixedMask) = vars.x(vars.M.allDims{:},vars.FixedMask);
            end
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
                if dataDims(i) >= 2+j && dataDims(i2) >= 2+j2 && vars.lambda(i,i2) > 0
                    % extract all 4 elements of each 2x2 matrix
                    subX1 = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:});
                    subX2 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:});
                    subX3 = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j2:2:2*Nt2+j2,post2Fill{:});
                    subX4 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j2:2:2*Nt2+j2,post2Fill{:});
                    [pX1,pX2,pX3,pX4] = vars.SecDiffMixProx(subX1,subX2,subX3,subX4,vars.lambda(i,i2));
                    y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:}) = pX1;
                    y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j2:2:2*Nt2+j2,post2Fill{:}) = pX2;
                    y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j2:2:2*Nt2+j2,post2Fill{:}) = pX3;
                    y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j2:2:2*Nt2+j2,post2Fill{:}) = pX4;
                    if ~isempty(vars.FixedMask)
                        y(vars.M.allDims{:},vars.FixedMask) = vars.x(vars.M.allDims{:},vars.FixedMask);
                    end
                end
            end
        end
    end
end
end