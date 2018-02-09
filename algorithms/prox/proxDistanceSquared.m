function [x,w] = proxDistanceSquared(M,f1,f2,lambda,FixedMask)
% proxDist(Mg,f,lambda)
% Proximal step towards f from given data g with parameter
% lambda on an arbitrary manifold. This is the proximal map of
% the distance function squared with fixed g.
% INPUT
%    M    : a manifold
%  f1,f2    : data point( sets/columns ) on the maifold M
%  lambda : stepsize towards f
%
% OPTIONAL
% FixedMask : fix a certain subset of f1
%
% OUTPUT
%       x : result point( sets) of the proximal map
%       w : indicator which terms where proxed, here: all
% ---
% MVIRT 1.1 | R. Bergmann | 2014-10-19
if nargin < 5
    FixedMask = [];
end
sizef = size(f1);
mD = sizef(1:length(M.ItemSize));
dD = sizef(length(M.ItemSize)+1:end);
if any( sizef ~= size(f2) )
    error('Input data has to be of same size')
end
%exclude equals and nans
mask = reshape(any(reshape(f1~=f2&~isnan(f1)&~isnan(f2) ,[prod(mD),dD]),1),[dD,1]);
% if only one term is unknown it can be intialized
% Init
x = f1;
if sum(mask(:))>0
    v = M.log(f1(M.allDims{:},mask),f2(M.allDims{:},mask));
    if sum(size(lambda))==2 % a nubmer
        t = lambda/(1+lambda);
    else
        %                t = repmat(lambda./(1+lambda),[ones(1,length(size(lambda))),this.ItemSize]);
        t = lambda./(1+lambda); %avoid repmat
        % manifolds first with one-dim stuff to avoid repmat
        t = shiftdim(-length(size(M.ItemSize)),t);
    end
    x(M.allDims{:},mask) = M.exp(f1(M.allDims{:},mask), t.*v);
    if ~isempty(FixedMask)
        x(M.allDims{:},FixedMask) = f1(M.allDims{:},FixedMask);
    end
end
% if nargout > 1 we return the mask additionally
if nargout > 1
    w = ones([dD,1]);
end
end