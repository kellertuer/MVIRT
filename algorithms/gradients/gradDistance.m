function eta = gradDistance(M,x,y,p)
% gradDistance(M,x,y,p) gradient of f(x) = 1/p d^p(x,y)for fixed y in M.
%
% INPUT
%  M   : a manfiold
% x,y  : two points on a manifold
%  p   : (2) exponent of the distance function
% ---
% MVIRT | R. Bergmann | 2018-01-22

if nargin < 4
    p=2;
end
mD = length(M.ItemSize);
if p ~= 2
    d = M.dist(y,x);
    eta = -M.log(x,y);
    eta = eta.*shiftdim(d.^(p-2),-mD);
    % securly avoid division by zero
    if any(d(:)<eps)
        eta(M.allDims{:},d<eps) = 0;
    end
else
    eta = -M.log(x,y);
end
end