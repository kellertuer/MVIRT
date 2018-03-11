function y = proxTV(varargin)
% proxTV(M,x,lambda) compute all proxima of TV in a cyclic manner.
% Items containing a NaN will be initalized with their first prox call.
%
% INPUT
%     M         : a manifold
%     x         : data (manifold-valued)
%   lambda      : parameter of the prox
%                   if given as a vector, an entry for each dimension,
%                   or a matrix, the diagonal for the TV terms,
%                   the offdiagonals for diagonal differences.
%
% OPTIONAL
%   'FixedMask' : binary mask the size of data items in x to indicate fixed
%   data items
%   'DifferenceProx' : (@(x1,x2,lambda)
%                                   proxAbsoluteDifference(M,x1,x2,lambda))
%                      specify a prox for the even/odd TV term proxes, i.e.
%                      switch the classical TV by Huber.
%
% OUTPUT
%     y     : result of the (iterated) proximal maps)
% ---
% Manifold-valued Image Restoration Toolbox 1.2 | R. Bergmann | 2018-02-09

ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addRequired(ip,'lambda');
addParameter(ip,'FixedMask',[]);
addParameter(ip,'DifferenceProx',[]);
parse(ip, varargin{:});
vars = ip.Results;
if isempty(vars.DifferenceProx)
    vars.DifferenceProx = @(x1,x2,lambda) proxAbsoluteDifference(vars.M,x1,x2,lambda);
end
sX = size(vars.x);
dataDims = sX( (length(vars.M.ItemSize)+1):end );
n = length(dataDims);
if sum(size(vars.lambda))==2 %lambda is a number
    vars.lambda = ones(1,n).*vars.lambda;
end
if isvector(vars.lambda)
    vars.lambda = diag(vars.lambda);
end
y = vars.x;
for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    for j=0:1
        % all starting with odd (even) indices
        Nt = floor((dataDims(i)-j)/2);
        if dataDims(i) >= 2+j && vars.lambda(i,i) > 0
            subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, postFill{:}); % odd
            subXe = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, postFill{:}); % even
            [pXo,pXe] = vars.DifferenceProx(subXo,subXe,vars.lambda(i,i));
            % write back
            y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, postFill{:}) = pXo;
            y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, postFill{:}) = pXe;
            if ~isempty(vars.FixedMask)
                y(vars.M.allDims{:},vars.FixedMask) = vars.x(vars.M.allDims{:},vars.FixedMask);
            end
        end
        for i2=i+1:n % diagonals
            interfill = repelem({':'},i2-i-1);
            post2Fill = repelem({':'},n-i2);
            Nt2 = floor((dataDims(i2)-j)/2);
            if dataDims(i) >= 2+j && dataDims(i2) >= 2+j
                if vars.lambda(i,i2) > 0
                % extract even/odd
                subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:});
                % +1,+1 diagonal
                subXd1 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:});
                [pXo,pXd1] = vars.DifferenceProx(subXo,subXd1,vars.lambda(i,i2));
                y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:}) = pXo;
                y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:}) = pXd1;
                if ~isempty(vars.FixedMask)
                    y(vars.M.allDims{:},vars.FixedMask) = vars.x(vars.M.allDims{:},vars.FixedMask);
                end
                end
                % +1 -1 diagonal
                % extract even/odd
                if vars.lambda(i2,i) > 0
                subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:});
                subXd1 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:});
                [pXo,pXd2] = vars.DifferenceProx(subXo,subXd1,vars.lambda(i2,i));
                y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:}) = pXo;
                y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:}) = pXd2;
                if ~isempty(vars.FixedMask)
                    y(vars.M.allDims{:},vars.FixedMask) = vars.x(vars.M.allDims{:},vars.FixedMask);
                end
                end
            end
        end
    end
end
end
