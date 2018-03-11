function [y,w] = proxDistance(M,x,f,lambda,FixedMask)
% proxDistanceSquared(M,x,f,lambda) prox of f(x) = d(x,y), y fixed,
% and prox paramater lambda. The fixedMask can be used to fix values, i.e.
% they are kept unchanged by hard reset.
%
% INPUT
%    M    : a manifold
%  x,f  : data point( sets/columns ) on the maifold M
%  lambda : stepsize towards f
%
% OPTIONAL
% FixedMask : fix a certain subset of x
%
% OUTPUT
%       x : result point( sets) of the proximal map
%       w : indicator which terms where proxed, here: all
% ---
% Manifold-valued Image Restoration Toolbox 1.0 ~ R. Bergmann, 2014-10-19
if nargin < 5
    FixedMask = [];
end
sizef = size(x);
mD = sizef(1:length(M.ItemSize));
dD = sizef(length(M.ItemSize)+1:end);
if any( sizef ~= size(f) )
    error('Input data has to be of same size')
end
%exclude equals and nans
mask = reshape(any(reshape(x~=f&~isnan(x)&~isnan(f) ,[prod(mD),dD,1]),1),[dD,1,1]);
% if only one term is unknown it can be intialized
% Init
y = x;
if sum(mask(:))>0
    v = M.log(x(M.allDims{:},mask),f(M.allDims{:},mask));
    d = M.dist(x(M.allDims{:},mask),f(M.allDims{:},mask));
    t = min(lambda./d,1);
    t = shiftdim(t,-length(size(M.ItemSize)));
    y(M.allDims{:},mask) = M.exp(x(M.allDims{:},mask), t.*v);
    if ~isempty(FixedMask) %hard reset
        y(M.allDims{:},FixedMask) = x(M.allDims{:},FixedMask);
    end
end
% if nargout > 1 we return the mask additionally
if nargout > 1
    w = ones([dD,1]);
end
end