%
% Illustrationg the gradient descent algorithm by 
% Computing the median of a set of points on the sphere
% =========================================================================
%
% This short demo illustrates the gradient descent algorithm to compute
% the Riemannian center of mass on the sphere
%
% ---
% Manifold-valued Image Restoration Toolbox 1.2
% Ronny Bergmann | 2018-01-28

%% Change folder, initialize
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
%%
run('../initMVIRT.m');
results = 'RiemannianMedianCyclicProx/';
mainC = [65, 161, 198]./255;
secondC = [255, 107, 76]/255;
%% We want to work on the 2-sphere
M = Sn(2);
pts = 30;
x = ArtificialS2SignalLemniscate(pts,'a',1);
figure(1);
plotS2(x,'o');
title('Given data on S2');
exportSpherePCT2Asy({x},{},{},{mainC},'File',[results,'InitialPoints.asy'],'Camera',[.3,.3,1]);

%% The function we would like to minize: p with smallest distance squared to all x
% and the proximal maps as well as our stopping criterion
F = @(p) 1/pts*sum( M.dist(repmat(p,[ones(size(M.ItemSize)),pts]),x) );
proximalMaps = cell(pts,1);
for i=1:pts
    proximalMaps{i} = @(y,lambda) proxDistance(M,y,x(M.allDims{:},i),lambda);
end
stoppingCriterion = stopCritMaxIterEpsilonCreator(M,2000,10^(-8));

%% the corresponding proximal point reads
[p1,recData1] = cyclicProximalPoint(M,x(M.allDims{:},1),proximalMaps,...
    stoppingCriterion,'lambda',pi/64,...
    'Record',@(x,xold,iter) [F(x);M.dist(x,xold)]);
[p2,recData2] = cyclicProximalPoint(M,x(M.allDims{:},1),proximalMaps,...
    stoppingCriterion,'Order','permute','lambda',pi/64,...
    'Record',@(x,xold,iter) [F(x);M.dist(x,xold)]);
%% Plot mean into figure
hold on
plotS2(p2,'o');
hold off
exportSpherePCT2Asy({x,p1},{},{},{mainC,secondC},'File',[results,'InitialAndMedian.asy'],'Camera',[.3,.3,1]);

%% Let's compare this to a subgradient descent
mD = length(M.ItemSize);
stoppingCriterion = stopCritMaxIterEpsilonCreator(M,2000,10^(-8));
gradF = @(p) sum(...
    -1./(shiftdim(... %the following line avoids division by zero
        pts*M.dist(repmat(p,[ones(1,mD),pts]),x) + ...
        double(M.dist(repmat(p,[ones(mD,1),pts]),x)==0),-mD)).*...
        M.log(repmat(p,[ones(mD,1),pts]),x),mD+1);
[p3,recData3] = subGradientDescent(M,x(M.allDims{:},1),F,gradF,...
    @(x,descentDir,iter,s) 1/iter,...
    stoppingCriterion,...
    'Record',@(x,xold,iter) [F(x);M.dist(x,xold)]);
%% Compare
hold on
plotS2(p3,'o');
hold off

%% PART II: Symmetric positive definite matrices

M = SymPosDef(3);

data = ArtificialSPDImage(pts,1.5); x = data(:,:,:,18);
figure(2);
[~,gD] = plotSPD(x,'EllipsoidPoints',30); %fine ellipses
title('Given SPD data');
% export and reuse the grid distance used abobe
exportSPD2Asymptote(x,'GridDistance',gD,'File',[results,'SPDdata.asy']);
%% Definitions look the same, they just have to save M in local scope again.
F = @(p) 1/pts*sum( M.dist(repmat(p,[ones(size(M.ItemSize)),pts]),x) );
proximalMaps = cell(pts,1);
for i=1:pts
    proximalMaps{i} = @(y,lambda) proxDistance(M,y,x(M.allDims{:},i),lambda);
end
stoppingCriterion = stopCritMaxIterEpsilonCreator(M,2000,10^(-8));

[p4,recData4] = cyclicProximalPoint(M,x(M.allDims{:},1),proximalMaps,...
    stoppingCriterion,'lambdaIterate',@(iter) pi/8/iter,...
    'Order','permute',...
    'Record',@(x,xold,iter) [F(x);M.dist(x,xold)]);

figure(3);
plotSPD(p4,'EllipsoidPoints',30);
exportSPD2Asymptote(p4,'GridDistance',gD,'File',[results,'SPDMedian.asy']);
