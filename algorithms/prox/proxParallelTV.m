function [y,w] = proxParallelTV(varargin)
% proxTV(M,x,lambda) - compute all proxima of the finite differences within
% the TV term in a cyclic manner (even/odd). Items containing a NaN will
% be initalized with their first prox call.
% and the results saved in indepentend arrays concatenated at the first
% signular dimension.
% INPUT
%     M         : a manifold
%     x         : an dataset from M^{m_1,m_2,...,m_n,K} of manifold-valued,
%                 data, where k equals 2*sum(lambda>0), i.e. the number of
%                  parallel proxes this term evaluates
%   lambda      : parameter of the proxes given as a vector, an entry for
%   each dimension, or a matrix, the diagonal for the TV terms, the
%   offdiagonals for diagonal differences.
%
% OPTIONAL
%   'FixedMask' : binary mask the size of data items in x to indicate fixed
%   data items
%   'DifferenceProx' : (@(x1,x2,lambda)
%   proxAbsoluteDifference(M,x1,x2,lambda))
%                      specify a prox for the even/odd TV term proxes, i.e.
%                      switch the classical TV by Huber.
%
% OUTPUT
%     y     : result of the (parallel) proximal maps
% ---
% MVIRT ? R. Bergmann ~ 2017-12-11

% TODO
% Diagonals?
ip = inputParser();
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'x');
addRequired(ip,'lambda');
addOptional(ip,'FixedMask',[]);
addOptional(ip,'DifferenceProx',[]);
parse(ip, varargin{:});
vars = ip.Results;
if isempty(vars.DifferenceProx)
    vars.DifferenceProx = @(x1,x2,lambda) proxAbsoluteDifference(vars.M,x1,x2,lambda);
end
sX = size(vars.x);
mL = length(vars.M.ItemSize);
if sum(size(vars.lambda))==2 %lambda is a number
    dataDims = sX( (length(vars.M.ItemSize)+1):end-1);
    n = length(dataDims);
    vars.lambda = ones(1,n).*vars.lambda;
end
if isvector(vars.lambda)
    vars.lambda = diag(vars.lambda);
end
numProxes = 2*sum(vars.lambda(:)>0);
if numProxes > 1
    dataDims = sX( (length(vars.M.ItemSize)+1):end-1);
    assert(numProxes== sX(end),['the number of parallel terms (last dimension, ',...
        num2str(sX(end)),' has to be twice the number of nonzero terms in lambda (',num2str(sum(vars.lambda>0)),...
        '.']);
else
    dataDims = sX( (length(vars.M.ItemSize)+1):end);
end
n = length(dataDims);
y = vars.x;
k = 0; %term to adress of the parallel terms
if nargout > 1
    w = ones([dataDims,numProxes]);
end
for i=1:n
    preFill = repelem({':'},i-1);
    postFill = repelem({':'},n-i);
    for j=0:1
        % all starting with odd (even) indices
        Nt = floor((dataDims(i)-j)/2);
        if dataDims(i) >= 2+j && vars.lambda(i,i) > 0
            k=k+1; %next prox
            subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, postFill{:},k); % odd
            subXe = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, postFill{:},k); % even
            [pXo,pXe] = vars.DifferenceProx(subXo,subXe,vars.lambda(i,i));
            % write back
            y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, postFill{:},k) = pXo;
            y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, postFill{:},k) = pXe;
            if ~isempty(vars.FixedMask)
                y(vars.M.allDims{:},vars.FixedMask,k) = vars.x(vars.M.allDims{:},vars.FixedMask,k);
            end
            if nargout > 1 %update w mask, exclude first if j==1
                if mod(dataDims(i),2)==0
                    if j==1 %divisible by 2 -> first and last not in odd terms
                        w(preFill{:},[1,dataDims(i)],postFill{:},k) = 0;
                    end
                else % for odd length:
                    if j==0% even: omit last
                        w(preFill{:},dataDims(i),postFill{:},k) = 0;
                    else % omit first
                        w(preFill{:},1,postFill{:},k) = 0;
                    end
                end
            end
        end
        for i2=i+1:n % diagonals
            interfill = repelem({':'},i2-i-1);
            post2Fill = repelem({':'},n-i2);
            Nt2 = floor((dataDims(i2)-j)/2);
            if dataDims(i) >= 2+j && dataDims(i2) >= 2+j
                if vars.lambda(i,i2) > 0
                    k=k+1; %next term
                    % extract even/odd
                    subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:},k);
                    % +1,+1 diagonal
                    subXd1 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:},k);
                    [pXo,pXd1] = vars.DifferenceProx(subXo,subXd1,vars.lambda(i,i2));
                    y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:},k) = pXo;
                    y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:},k) = pXd1;
                    if ~isempty(vars.FixedMask)
                        y(vars.M.allDims{:},vars.FixedMask,k) = vars.x(vars.M.allDims{:},vars.FixedMask,k);
                    end
                    if nargout > 1
                        warning('Mask w not jet updated for diagonal terms!');
                    end
                end
                % +1 -1 diagonal
                % extract even/odd
                if vars.lambda(i2,i) > 0
                    k=k+1;
                    subXo = y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:},k);
                    subXd1 = y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:},k);
                    [pXo,pXd2] = vars.DifferenceProx(subXo,subXd1,vars.lambda(i2,i));
                    y(vars.M.allDims{:},preFill{:}, 1+j:2:2*Nt+j, interfill{:}, 2+j:2:2*Nt2+j,post2Fill{:},k) = pXo;
                    y(vars.M.allDims{:},preFill{:}, 2+j:2:2*Nt+j, interfill{:}, 1+j:2:2*Nt2+j,post2Fill{:},k) = pXd2;
                    if ~isempty(vars.FixedMask)
                        y(vars.M.allDims{:},vars.FixedMask,k) = vars.x(vars.M.allDims{:},vars.FixedMask,k);
                    end
                    if nargout > 1
                        warning('Mask w not jet updated for diagonal terms!');
                    end
                end
            end
        end
    end
end
end