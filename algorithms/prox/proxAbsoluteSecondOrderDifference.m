function [x1,x2,x3] = proxAbsoluteSecondOrderDifference(M,f1,f2,f3,lambda)
sizef = size(f1);
mD = sizef(1:length(M.ItemSize));
mL = length(mD);
dD = sizef(length(M.ItemSize)+1:end);
if any(sizef ~= size(f2)) || any(sizef ~= size(f3))
    error('Input data has to be of same size')
end
mask = reshape(any(reshape(~isnan(f1)&~isnan(f2)&~isnan(f3),[prod(mD),dD]),1)&...
    shiftdim(M.dist(M.midPoint(f1,f3),f2)>10^(-15),-1),[dD,1]);
%only 1 unknown
umaskf1 = reshape(any(reshape(~isnan(f2)&~isnan(f3)&isnan(f1),[prod(mD),dD]),1),[dD,1]);
umaskf2 = reshape(any(reshape(~isnan(f3)&~isnan(f1)&isnan(f2),[prod(mD),dD]),1),[dD,1]);
umaskf3 = reshape(any(reshape(~isnan(f1)&~isnan(f2)&isnan(f3),[prod(mD),dD]),1),[dD,1]);
% Init
x1 = f1;
if any(umaskf1(:))
    x1(M.allDims{:},umaskf1) = M.geopoint(f3(M.allDims{:},umaskf1),f2(M.allDims{:},umaskf1),2);
end
x2 = f2;
if any(umaskf2(:))
    x2(M.allDims{:},umaskf2) = M.geopoint(f1(M.allDims{:},umaskf2),f3(M.allDims{:},umaskf2),0.5);
end
x3 = f3;
if any(umaskf3(:))
    x3(M.allDims{:},umaskf3) = M.geopoint(f1(M.allDims{:},umaskf3),f2(M.allDims{:},umaskf3),2);
end
if any(mask(:)) % we have equal nonequal values -> compute proxes
    t1 = permute(f1(M.allDims{:},mask),[1:mL,mL+2,mL+1]);
    t2 = permute(f2(M.allDims{:},mask),[1:mL,mL+2,mL+1]);
    t3 = permute(f3(M.allDims{:},mask),[1:mL,mL+2,mL+1]);
    if isa(M,'Rn') % closed forms known
        m = min(lambda,abs(t1-2*t2+t3)/6);
        s = sign(t1-2*t2+t3);
        xt = cat(mL+1, t1 - s.*m, t2 + 2*s.*m, t3 - s.*m);
    elseif isa(M,'S1')
        m = min(lambda,abs(symMod(t1-2*t2+t3,2*pi))/6);
        s = sign(symMod(t1-2*t2+t3,2*pi));
        xt = cat(mL+1, t1 - s.*m, t2 + 2*s.*m, t3 - s.*m);
        xt = symMod(xt,2*pi);
    elseif isa(M,'S1mRn')
        prem = t1-2*t2+t3;
        prem(M.SComponents,:,:) = symMod(prem(M.SComponents,:,:),2*pi);
        m = min(lambda,abs(prem)/6);
        s = sign(prem);
        xt = cat(mL+1, t1 - s.*m, t2 + 2*s.*m, t3 - s.*m);
        xt(M.SComponents,:,:) = symMod(xt(M.SComponents,:,:),2*pi);
    else %gradient descent
        xInit = cat(mL+1,t1,t2,t3);
        gradF = @(x) cat(mL+1,...%concatenate gradx,grady,gradz
            ... grad X:
            -2*M.log(x(M.allDims{:},1,:),t1) + lambda*M.AdjDxGeo(...
            x(M.allDims{:},1,:),...
            x(M.allDims{:},3,:),1/2,...
            gradxDist(M,...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:)),...
            x(M.allDims{:},2,:)...
            )),...
            ...% grad y - just the log but opposite direction
            -2*M.log(x(M.allDims{:},2,:),t2) + lambda*gradxDist(M,...
            x(M.allDims{:},2,:),...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:))...
            ),...
            ... grad Z:
            -2*M.log(x(M.allDims{:},3,:),t3) + lambda*M.AdjDxGeo(...
            x(M.allDims{:},3,:),...
            x(M.allDims{:},1,:),1/2,...
            gradxDist(M,...
            M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:)),...
            x(M.allDims{:},2,:)...
            ))...
            );
        F = @(x) M.dist(x(M.allDims{:},1,:),t1).^2 + ...
            M.dist(x(M.allDims{:},2,:),t2).^2 + ...
            M.dist(x(M.allDims{:},3,:),t3).^2 + ...
            lambda*M.dist(x(M.allDims{:},2,:),M.midPoint(x(M.allDims{:},1,:),x(M.allDims{:},3,:)));
        stepSize = @(x,eta,iter,old_step) 1/iter;
        stopCrit = stopCritMaxIterEpsilonCreator(M,10,10^(-9)); %
        xt = subGradientDescent(M,xInit,F,gradF,stepSize,stopCrit);
    end
    x1(M.allDims{:},mask) = xt(M.allDims{:},1,:);
    x2(M.allDims{:},mask) = xt(M.allDims{:},2,:);
    x3(M.allDims{:},mask) = xt(M.allDims{:},3,:);
end
end

