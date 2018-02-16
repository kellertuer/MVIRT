%
% Illustrationg the gradient descent algorithm by 
% Computing the mean of a set of points on the sphere
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
results = 'RiemannianMeanGradDesc/';
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
F = @(p) 1/pts*sum( M.dist(repmat(p,[ones(size(M.ItemSize)),pts]),x).^2 );
% The corresponding gradient reads
gradF = @(p) sum( -2/pts*M.log(repmat(p,[ones(size(M.ItemSize)),pts]),x),length(M.ItemSize)+1);
stoppingCriterion = stopCritMaxIterEpsilonCreator(M,100,10^(-9));
stepSizeRule = @(x,eta,iter,initial) 1/2;
[p1,recData1] = gradientDescent(M,x(M.allDims{:},1),gradF,stepSizeRule,...
    stoppingCriterion,'Record',@(x,xold,iter) [F(x);iter]);

%% Plot mean into figure
hold on
plotS2(p1,'o');
hold off
exportSpherePCT2Asy({x,p1},{},{},{mainC,secondC},'File',[results,'InitialAndMean.asy'],'Camera',[.3,.3,1]);

%% Let's improve on the choice of the stepsize
stepSizeRuleA = @(x,eta,iter,initial) ...
    stepSizeArmijo(M,F,x,eta,'InitialStepSize',initial,'rho',0.99,'c',0.5);
[p2,recData2] = gradientDescent(M,x(M.allDims{:},1),gradF,stepSizeRuleA,...
    stoppingCriterion,'Record',@(x,xold,iter) [F(x);iter]);

%% Compare
M.dist(p1,p2)

%% TEIL II: Symmetric positive definite matricesPD

M = SymPosDef(3);

data = ArtificialSPDImage(pts,1.5); x = data(:,:,:,18);
figure(2);
[~,gD] = plotSPD(x,'EllipsoidPoints',30); %fine ellipses
title('Given SPD data');
% export and reuse the grid distance used abobe
exportSPD2Asymptote(x,'GridDistance',gD,'File',[results,'SPDdata.asy']);
%% Definitions look the same, they just have to save M in local scope again.
F = @(p) 1/pts*sum( M.dist(repmat(p,[ones(size(M.ItemSize)),pts]),x).^2 );
gradF = @(p) sum( -2/pts*M.log(repmat(p,[ones(size(M.ItemSize)),pts]),x),length(M.ItemSize)+1);
stoppingCriterion = stopCritMaxIterEpsilonCreator(M,1000,10^(-9));

[p3,recData3] = gradientDescent(M,x(M.allDims{:},1),gradF,stepSizeRule,...
    stoppingCriterion,'Record',@(x,xold,iter) [F(x);iter]);

figure(3);
plotSPD(p3,'EllipsoidPoints',30);
exportSPD2Asymptote(p3,'GridDistance',gD,'File',[results,'SPDMean.asy']);
