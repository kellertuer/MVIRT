%%       Denoising of an artificial S2 image
%
%   Denoising of an artificial S2 image with the NL MMSE denoising method
%   in comparison to not using a second order update, a simple NL Mean 
%   approach, TV and TV&TV2 methods
%
%
%   This example is published in Sect 5.6 of
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%       Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
%
%
%  This file can be started without any changes; it initializes the Toolbox
%  itself
%
% Manifold-valued Image Restoration Toolbox 1.2
%  J. Persch  ~ 2016-07-05 | R. Bergmann ~ 2017-01-25
% see LICENSE.txt
run('../../initMVIRT.m')
%% Debug
setDebugLevel('LevelMin',0);
setDebugLevel('LevelMax',100);
setDebugLevel('text',3);
setDebugLevel('time',3);
setDebugLevel('LoadData',1);
setDebugLevel('WriteImages',1);
%
M = Sn(2);
sigma = 0.3;
pts = 64;
%
% Get data
resultsFolder = ['S2WhirlImage',filesep];
if getDebugLevel('LoadData')
    load([resultsFolder,'S2Whirl.mat'],'u_noisy_S2','u_0_S2');
else    
    u_0_S2 = S2WhirlImage;
    u_noisy_S2 = M.addNoise(u_0_S2,sigma);
end
% Parameters

% Initialize problem struct
clear problem
problem.M = M;
problem.f = u_noisy_S2;
problem.sigma =sigma;
% optimized gamma on 1/20 N, K on N and p = 3,5,7
problem.gamma = .8;
problem.patch_size_1 = 3;
problem.window_size_1 = 127;
problem.K_1 = 65;

problem.patch_size_2 = 5;
problem.window_size_2 = 127;
problem.K_2 = 54;
% Parameters for TV1&2 CPPA TV 1/100 N TV2 1/10N
problem.lambda = pi/2;
problem.alpha = 0.18;
problem.beta = 2.6;
%%
[u_final,u_oracle] = NL_MMSE_2D(problem);
tic
% patch+neighbor on N, sigmas on 1/10N
u_nl = NL_Mean(M,u_noisy_S2,23,127,104,1.5,.2);
toc
u_tv12 = cppa_ad_2D(problem);
problem.beta = 0;
problem.alpha = 0.24;
u_tv1 = cppa_ad_2D(problem);

% Parameters for no update
problem.K_1 = 6;
problem.patch_size_1 = 5;
problem.window_size_1 = 127;
problem.oracle = false;
problem.simple = true;
u_no_update = NL_MMSE_2D(problem);       
        
arts2_mean_err_noisy = sum(sum(M.dist(u_0_S2,u_noisy_S2).^2))/pts^2;
arts2_mean_err_oracle = sum(sum(M.dist(u_0_S2,u_oracle).^2))/pts^2;
arts2_mean_err_final = sum(sum(M.dist(u_0_S2,u_final).^2))/pts^2;
arts2_mean_err_nl = sum(sum(M.dist(u_0_S2,u_nl).^2))/pts^2;
arts2_mean_err_tv12 = sum(sum(M.dist(u_0_S2,u_tv12).^2))/pts^2;
arts2_mean_err_tv1 = sum(sum(M.dist(u_0_S2,u_tv1).^2))/pts^2;
arts2_mean_err_no_update = sum(sum(M.dist(u_0_S2,u_no_update).^2))/pts^2;
%%
if getDebugLevel('WriteImages') == 1
    exportSphere2Asymptote(u_0_S2,'ExportHeader',true,'File',[resultsFolder,'arts2_orig.asy']);
    exportSphere2Asymptote(u_noisy_S2,'ExportHeader',true,'File',[resultsFolder,'arts2_noisy.asy']);
    exportSphere2Asymptote(u_oracle,'ExportHeader',true,'File',[resultsFolder,'arts2_oracle.asy']);
    exportSphere2Asymptote(u_final,'ExportHeader',true,'File',[resultsFolder,'arts2_final.asy']);
    exportSphere2Asymptote(u_nl,'ExportHeader',true,'File',[resultsFolder,'arts2_nl.asy']);
    exportSphere2Asymptote(u_tv12,'ExportHeader',true,'File',[resultsFolder,'arts2_tv12.asy']);
    exportSphere2Asymptote(u_tv1,'ExportHeader',true,'File',[resultsFolder,'arts2_tv1.asy']);
    exportSphere2Asymptote(u_no_update,'ExportHeader',true,'File',[resultsFolder,'arts2_no_update.asy']);
    save([resultsFolder,'arts_s2_me.mat'],'arts2_mean_err_noisy','arts2_mean_err_oracle','arts2_mean_err_final','arts2_mean_err_nl','arts2_mean_err_tv12','arts2_mean_err_tv1','arts2_mean_err_no_update');
end

