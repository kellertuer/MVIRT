function [x1,x2,x3,x4] = proxAbsoluteSecondOrderMixedDifference(M,f1,f2,f3,f4,lambda)
% proxAbsoluteSecondOrderMixedDifference(M,f1,f2,f3,f4,lambda)
% Compute the second order mixed difference based on the mid point model
% employing a sub gradient descent
%
% INPUT
%    M          : a manifold
%   f1,f2,f3,f4 : Data items from a 2x2 matrix each
%   lambda      : parameter of the proximal map
%
% OUTPUT
%  x1,x2,x3,x4 : result of the proximal map
% ---
% MVIRT | R. Bergmann, 2018-01-22
sizef = size(f1);
mD = sizef(1:length(M.ItemSize));
mL = length(mD);
dD = sizef(length(M.ItemSize)+1:end);
if any(sizef ~= size(f2)) || any(sizef ~= size(f3)) || any(sizef ~= size(f4))
    error('Input data has to be of same size')
end
% which ones do we have to compute?
mask = reshape(any(reshape(~isnan(f1)&~isnan(f2)&~isnan(f3)&~isnan(f4),[prod(M.ItemSize),dD]),1)&...
        shiftdim(M.dist(M.midPoint(f1,f3),M.midPoint(f2,f4))>10^(-15),-1),[dD,1]);
% exactely one is unknown
umaskf1 = reshape(any(reshape(~isnan(f2)&~isnan(f3)&~isnan(f4)&isnan(f1),[prod(mD),dD]),1),[dD,1]);
umaskf2 = reshape(any(reshape(~isnan(f3)&~isnan(f1)&~isnan(f4)&isnan(f2),[prod(mD),dD]),1),[dD,1]);
umaskf3 = reshape(any(reshape(~isnan(f1)&~isnan(f2)&~isnan(f4)&isnan(f3),[prod(mD),dD]),1),[dD,1]);
umaskf4 = reshape(any(reshape(~isnan(f1)&~isnan(f2)&~isnan(f3)&isnan(f4),[prod(mD),dD]),1),[dD,1]);
% Init
x1 = f1;
if any(umaskf1(:))
    x1(M.allDims{:},umaskf1) = M.geopoint(f3(M.allDims{:},umaskf1),M.geopoint(f2(M.allDims{:},umaskf1),f4(M.allDims{:},umaskf1),0.5),2);
end
x2 = f2;
if any(umaskf2(:))
    x2(M.allDims{:},umaskf2) = M.geopoint(f4(M.allDims{:},umaskf2),M.geopoint(f1(M.allDims{:},umaskf2),f3(M.allDims{:},umaskf2),0.5),2);
end
x3 = f3;
if any(umaskf3(:))
    x3(M.allDims{:},umaskf3) = M.geopoint(f1(M.allDims{:},umaskf3),M.geopoint(f2(M.allDims{:},umaskf3),f4(M.allDims{:},umaskf3),0.5),2);
end
x4 = f4;
if any(umaskf4(:))
    x4(M.allDims{:},umaskf4) = M.geopoint(f2(M.allDims{:},umaskf4),M.geopoint(f1(M.allDims{:},umaskf4),f3(M.allDims{:},umaskf4),0.5),2);
end
if any(mask(:)) % we have equal nonequal values -> compute proxes
    t1 = permute(f1(M.allDims{:},mask),[1:mL,mL+2,mL+1]);
    t2 = permute(f2(M.allDims{:},mask),[1:mL,mL+2,mL+1]);
    t3 = permute(f3(M.allDims{:},mask),[1:mL,mL+2,mL+1]);
    t4 = permute(f4(M.allDims{:},mask),[1:mL,mL+2,mL+1]);
    if isa(M,'Rn') % closed forms known
        m = min(lambda,abs(t1-t2-t3+t4)/4);
        s = sign(t1-t2-t3+t4);
        xt = cat(mL+1, t1 - s.*m, t2 + s.*m, t3 + s.*m, t4 - s.*m);
    elseif isa(M,'S1')
        m = min(lambda,abs(symMod(t1-t2-t3+t4,2*pi))/4);
        s = sign(symMod(t1-t2-t3+t4,2*pi));
        xt = cat(mL+1, t1 - s.*m, t2 + s.*m, t3 + s.*m, t4 - s.*m);
        xt = symMod(xt,2*pi);
    elseif isa(M,'S1mRn')
        prem = t1-t2-t3+t4;
        prem(M.SComponents,:,:) = symMod(prem(M.SComponents,:,:),2*pi);
        m = min(lambda,abs(prem)/4);
        s = sign(prem);
        xt = cat(mL+1, t1 - s.*m, t2 + s.*m, t3 + s.*m, t4 - s.*m);
        xt(M.SComponents,:,:) = symMod(xt(M.SComponents,:,:),2*pi);
    else %gradient descent
        xInit = cat(mL+1,t1,t2,t3,t4);
        % the data is gicen by f1 f2
        %                      f3 f4
        % such that the functional is d(c(f1,f3),c(f2,f4))
        gradF = @(x) cat(mL+1,...%concatenate gradx1,2,3,4
            ... grad x1:
            -2*M.log(x(M.allDims{:},1,:),t1) + lambda*M.AdjDxGeo(...
            x(M.allDims{:},1,:),...
            x(M.allDims{:},3,:),1/2,...
            gradxDist(M,...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:)),...
            M.midPoint(x(M.allDims{:},2,:),x(M.allDims{:},4,:))...
            )),...
            ... grad x2:
            -2*M.log(x(M.allDims{:},2,:),t2) + lambda*M.AdjDxGeo(...
            x(M.allDims{:},2,:),...
            x(M.allDims{:},4,:),1/2,...
            gradxDist(M,...
            M.midPoint(x(M.allDims{:},2,:),x(M.allDims{:},4,:)),...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:))...
            )),...
            ... grad x3:
            -2*M.log(x(M.allDims{:},3,:),t3) + lambda*M.AdjDxGeo(...
            x(M.allDims{:},3,:),...
            x(M.allDims{:},1,:),1/2,...
            gradxDist(M,...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:)),...
            M.midPoint(x(M.allDims{:},2,:),x(M.allDims{:},4,:))...
            )),...
            ... grad x4:
            -2*M.log(x(M.allDims{:},4,:),t4) + lambda*M.AdjDxGeo(...
            x(M.allDims{:},4,:),...
            x(M.allDims{:},2,:),1/2,...
            gradxDist(M,...
            M.midPoint(x(M.allDims{:},2,:),x(M.allDims{:},4,:)),...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:))...
            ))...
            );
        F = @(x) M.dist(x(M.allDims{:},1,:),t1).^2 + ...
            M.dist(x(M.allDims{:},2,:),t2).^2 + ...
            M.dist(x(M.allDims{:},3,:),t3).^2 + ...
            M.dist(x(M.allDims{:},4,:),t4).^2 + ...
            lambda*M.dist(...
                M.midPoint(x(M.allDims{:},2,:),x(M.allDims{:},4,:)),...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:)) );
        stepSize = @(x,eta,iter,old_step) 1/iter;
        stopCrit = stopCritMaxIterEpsilonCreator(M,3,10^(-9));
        xt = subGradientDescent(M,xInit,F,gradF,stepSize,stopCrit);
    end
    x1(M.allDims{:},mask) = xt(M.allDims{:},1,:);
    x2(M.allDims{:},mask) = xt(M.allDims{:},2,:);
    x3(M.allDims{:},mask) = xt(M.allDims{:},3,:);
    x4(M.allDims{:},mask) = xt(M.allDims{:},4,:);
end
end

