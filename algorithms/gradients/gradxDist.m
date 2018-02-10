function xi = gradxDist(M,x,y)
% gradxDist(M,x,y) compute the (sub)gradient of d(x,y) with respect to x
%
% INPUT
%    M : a manifold
%  x,y : two points or sets of points on M
%
% OUTPUT
%   eta : the gradient
mL = length(M.ItemSize);
d = M.dist(x,y);
dZ = d < 10^(-7);
xi = -M.log(x,y).*shiftdim(~dZ,-mL)./shiftdim(d + double(dZ), -mL);
end

