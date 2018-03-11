function [x1,x2] = proxAbsoluteDifferenceHuber(M,f1,f2,lambda,tau,omega)
% proxAbsoluteDifferenceHuber(M,f1,f2,lambda,tau,omega) prox d(x,y) relaxed
% Given the Huber function
% h(t) = tau^2t^2 for t< omega/(sqrt(2)tau), omega*sqrt(2)*tau*t - omega^2/2
% i.e. for small t a suqared function (steered by tau) and for large t a
% linear part with ascent omega.
% lambda is the prox parameter, as for the absolute difference, this
% function inpaints NaN values and works on arbitrary manifolds
%
% INPUT
%   M         : A manifold
%  f1,f2   : data columns
%  lambda     : parameter within the proximal map
%  tau, omega : parameters of the Huber function
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
    d = M.dist(t1,t2);
    step = 2*lambda*tau^2/(1+4*lambda.*tau.^2)*ones(size(d)); % small case: squared (cf. ADSquared)
    % filter linear part
    step(d >= omega*(1+4*lambda.*tau.^2)./(sqrt(2)*tau) ) = ...
        min(1/2, sqrt(2)*lambda*omega*tau./d(d >= omega*(1+4*lambda.*tau.^2)./(sqrt(2)*tau) ) );
    step = permute(repmat(step(:),[1,M.ItemSize]), length(M.ItemSize)+1:-1:1);
    x1(M.allDims{:},mask) = M.exp(t1,step.*M.log(t1,t2));
    x2(M.allDims{:},mask) = M.exp(t2,step.*M.log(t2,t1));
end
end