function [x1,x2] = proxAbsoluteDifference(M,f1,f2,lambda)
% proxAbsoluteDifference(M,f1,f2,lambda)
% Proximal steps for the absolute difference (TV summand) of f1 and f2 with
% parameter lambda on an arbitrary manifold. This is
% the proximal map of the piointwise distance function d(f1,f2).
%
% INPUT
%   M         : A manifold
%  f1,f2   : data columns
%  lambda     : parameter within the proximal map
%
% OUTPUT
%  x1,x2   : resulting columns of the proximal map
% ---
% Manifold-valued Image Restoration Toolbox 1.0 ~ R. Bergmann, 2014-10-19
sizef = size(f1);
mD = sizef(1:length(M.ItemSize));
dD = sizef(length(M.ItemSize)+1:end);
if any(sizef ~= size(f2))
    error('Input data has to be of same size')
end
%exclude equals and nans
mask = reshape(any(reshape(f1~=f2&~isnan(f1)&~isnan(f2) ,[prod(mD),dD]),1),[dD,1]);
% if only one term is unknown it can be intialized
umaskf1 = reshape(any(reshape(isnan(f1)&~isnan(f2),[prod(mD),dD]),1),[dD,1]);
umaskf2 = reshape(any(reshape(isnan(f2)&~isnan(f1),[prod(mD),dD]),1),[dD,1]);
% Init
x1 = f1;
if any(umaskf1(:))
    x1(M.allDims{:},umaskf1) = f2(M.allDims{:},umaskf1); % initialized unknowns
end
x2 = f2;
if any(umaskf2(:))
    x2(M.allDims{:},umaskf2) = f1(M.allDims{:},umaskf2);
end
if any(mask(:)) % we have equal nonequal values -> compute proxes
    t1 = f1(M.allDims{:},mask);
    t2 = f2(M.allDims{:},mask);
    step = min(1/2,lambda./M.dist(t1,t2));
    step = permute(repmat(step(:),[1,M.ItemSize]), length(M.ItemSize)+1:-1:1);
    x1(M.allDims{:},mask) = M.exp(t1,step.*M.log(t1,t2));
    x2(M.allDims{:},mask) = M.exp(t2,step.*M.log(t2,t1));
end
end