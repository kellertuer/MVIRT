function [x1,x2] = proxAbsoluteDifferenceSquared(M,f1,f2,lambda)
% proxAbsoluteDifference(M,f1,f2,lambda) prox of d^2_M(f1,f2) in both args
% with parameter lambda on an arbitrary manifold. Values containing NaN
% are initialized to the other argument, which is the minimizer.
%
% INPUT
%   M       : A manifold
%  f1,f2    : data columns
%  lambda   : proxParameter
%
% OUTPUT
%  x1,x2    : resulting columns of the proximal map
% ---
% Manifold-valued Image Restoration Toolbox 1.0 ~ R. Bergmann, 2014-10-19
sizef = size(f1);
mD = sizef(1:length(M.ItemSize));
dD = sizef(length(M.ItemSize)+1:end);
if any(sizef ~= size(f2))
    error('Input data has to be of same size')
end
%exclude equals and nans
mask = reshape(any(reshape(f1~=f2&~isnan(f1)&~isnan(f2) ,[prod(mD),dD]),1),[ones(length(mD),1),dD]);
% if only one term is unknown it can be intialized
umaskf1 = reshape(any(reshape(isnan(f1)&~isnan(f2),[prod(mD),dD]),1),[ones(length(mD),1),dD]);
umaskf2 = reshape(any(reshape(isnan(f2)&~isnan(f1),[prod(mD),dD]),1),[ones(length(mD),1),dD]);
% Init
x1 = f1;
if sum(umaskf1(:)) > 0
    x1(M.allDims{:},umaskf1) = f2(M.allDims{:},umaskf1); % initialized unknowns
end
x2 = f2;
if sum(umaskf2(:)) > 0
    x2(M.allDims{:},umaskf2) = f1(M.allDims{:},umaskf2);
end
if any(mask(:)) % we have equal nonequal values -> compute proxes
    t1 = f1(M.allDims{:},mask);
    t2 = f2(M.allDims{:},mask);
    step = lambda/(1+2*lambda)*M.dist(t1,t2);
    step = permute(repmat(step(:),[1,M.ItemSize]), length(M.ItemSize)+1:-1:1);
    x1(M.allDims{:},mask) = M.exp(t1,step.*M.log(t1,t2));
    x2(M.allDims{:},mask) = M.exp(t2,step.*M.log(t2,t1));
end
end