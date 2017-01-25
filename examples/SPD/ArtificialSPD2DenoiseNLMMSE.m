%%       Denoising of an artificial SPD(2) valued image
%
%   Denoising of an artificial SPD(2) image with the NL MMSE denoising
%   method, with same and different parameters in both steps.
%
%
%   This example is published in Sect 5.4 of
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%       Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
%
%  This file can be started without any changes; it initializes the Toolbox
%  itself
%
% Manifold-valued Image Restoration Toolbox 1.2
%  J. Persch | R. Bergmann ~ 2017-01-25
% see LICENSE.txt
run('../../initMVIRT.m')
%% Debug
setDebugLevel('LevelMin',0);
setDebugLevel('LevelMax',100);
setDebugLevel('text',3);
setDebugLevel('time',3);
setDebugLevel('LoadData',1);
setDebugLevel('WriteImages',1);
setDebugLevel('Figures',1);
% Parameters for image and noise
sigma = 0.2;
pts = 65;
M = SymPosDef(2);
resultsFolder = ['ArtificialSPD2',filesep];
if getDebugLevel('NewNoise') == 1
    T_orig = ArtificialSPD2Image(pts);
    T = M.addNoise(T_orig,sigma,'Distribution','Gaussian');
else
    load([resultsFolder,'SPD2Test.mat'],'T','T_orig');
end

% Initialize problem
clear problem
problem.M = M;
problem.f = T;
problem.sigma = sigma;
problem.gamma = 1;

problem.patch_size_1 = 9;
problem.window_size_1 = 115;

problem.patch_size_2 = 9;
problem.window_size_2 = 115;

problem.K_1 = 1038;
problem.K_2 = 1038;
[u_final, u_oracle] = NL_MMSE_2D(problem);
problem.patch_size_2 = 7;
problem.window_size_2 = 41;
problem.K_2 = 193;
[u_final_sw, u_oracle_sw] = NL_MMSE_2D(problem);
% TV alpha in 1/100N
problem.lambda = 1;
problem.alpha = 0.25;
problem.beta = 0;
u_tv = cppa_ad_2D(problem);
% NL Mean
tic
u_nl = NL_Mean(M,T,33,9,81,2,.2);
toc
% calculate erros
SPD2_mean_err_tv = sum(sum(M.dist(T_orig,u_tv).^2))/pts^2;
SPD2_mean_err_nl = sum(sum(M.dist(T_orig,u_nl).^2))/pts^2;
SPD2_mean_err_noisy = sum(sum(M.dist(T_orig,T).^2))/pts^2;
SPD2_mean_err_oracle = sum(sum(M.dist(T_orig,u_oracle).^2))/pts^2;
SPD2_mean_err_final = sum(sum(M.dist(T_orig,u_final).^2))/pts^2;
SPD2_mean_err_oracle_sw = sum(sum(M.dist(T_orig,u_oracle_sw).^2))/pts^2;
SPD2_mean_err_final_sw = sum(sum(M.dist(T_orig,u_final_sw).^2))/pts^2;
%% Display Resutls
if getDebugLevel('Figures') == 1
    figure(1)
    plotSPD2(T,'GridDistance',4);
    figure(3)
    plotSPD2(u_final_sw,'GridDistance',4);
    title('Different Parameters');
    figure(4)
    plotSPD2(u_final,'GridDistance',4);
    title('Same Parameters');    
    figure(5)
    plotSPD2(u_nl,'GridDistance',4);
    title('NL Means');
end
%% Export Results
if getDebugLevel('WriteImages') ==1
   exportEllipses2Tikz(T_orig,'File',[resultsFolder,'SPD2orig.tex'],'ExportHeader',true,'GridDistance',4);
   exportEllipses2Tikz(T,'File',[resultsFolder,'SPD2noisy.tex'],'ExportHeader',true,'GridDistance',4);
   exportEllipses2Tikz(u_oracle,'File',[resultsFolder,'SPD2oracle.tex'],'ExportHeader',true,'GridDistance',4);
   exportEllipses2Tikz(u_final,'File',[resultsFolder,'SPD2final.tex'],'ExportHeader',true,'GridDistance',4);
   exportEllipses2Tikz(u_tv,'File',[resultsFolder,'SPD2tv.tex'],'ExportHeader',true,'GridDistance',4); 
   exportEllipses2Tikz(u_nl,'File',[resultsFolder,'SPD2nl.tex'],'ExportHeader',true,'GridDistance',4);
   exportEllipses2Tikz(u_oracle_sw,'File',[resultsFolder,'SPD2oracle_c.tex'],'ExportHeader',true,'GridDistance',4);
   exportEllipses2Tikz(u_final_sw,'File',[resultsFolder,'SPD2final_c.tex'],'ExportHeader',true,'GridDistance',4);     
end
