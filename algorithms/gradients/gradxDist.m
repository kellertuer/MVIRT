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
xi = -M.log(x,y)./shiftdim(d + double(d==0), -mL);
end

