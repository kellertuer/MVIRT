%%       Denoising of an artificial S1 image
%
%   Denoising of an artificial S1 image with the NL MMSE denoising method
%   in comparison to a simple NL Mean approach, for the first/second order
%   TV restorations run the S1 Example "syntheticSARImage"
%
%
%   This example is published in Sect 5.5 of
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%       Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
%
%  This file can be started without any changes; it initializes the Toolbox
%  itself
%
% ---
% Manifold-valued Image Restoration Toolbox 1.2
%  J. Persch  ~ 2016-07-05 | R. Bergmann ~ 2017-01-07
% see LICENSE.txt
run('../../initMVIRT.m')
resultsFolder = ['syntheticSARImage',filesep];
%% Debug
setDebugLevel('LevelMin',0);
setDebugLevel('LevelMax',100);
setDebugLevel('text',3);
setDebugLevel('time',3);
setDebugLevel('LoadData',1);
setDebugLevel('WriteImages',1);
setDebugLevel('Figures',1);
% Parameters for image and noise
sigma = 0.3;
pts = 256;
[u_noisy_SAR,u_0_SAR] = ArtificalSARImage(pts,sigma);
if getDebugLevel('LoadData')
    load([resultsFolder,'S1-syntheticSARImageN256.mat'],'Zn');
    u_noisy_SAR = Zn;
end
u_noisy_Sn1 = zeros(2,size(u_noisy_SAR,1),size(u_noisy_SAR,2));
u_orig = zeros(2,size(u_noisy_SAR,1),size(u_noisy_SAR,2));
u_noisy_Sn1(1,:,:) = cos(u_noisy_SAR);
u_noisy_Sn1(2,:,:) = sin(u_noisy_SAR);
u_orig(1,:,:) = cos(u_0_SAR);
u_orig(2,:,:) = sin(u_0_SAR);

M = Sn(1);
% Initialize problem
clear problemImage
problem.M = M;
problem.f = u_noisy_Sn1;
problem.sigma =sigma;
problem.gamma = 1.1;
problem.patch_size_1 = 9;
problem.window_size_1 = 119;
problem.patch_size_2 = 7;
problem.window_size_2 = 123;
problem.K_1 = 186;
problem.K_2 = 86;
%% Start Denoising
[u_final,u_oracle] = NL_MMSE_2D(problem);
u_nl = NL_Mean(M,u_noisy_Sn1,11,23,33,46,.2);
arts1_mean_err_noisy = sum(sum(M.dist(u_orig,u_noisy_Sn1).^2))/pts^2;
arts1_mean_err_oracle = sum(sum(M.dist(u_orig,u_oracle).^2))/pts^2;
arts1_mean_err_final = sum(sum(M.dist(u_orig,u_final).^2))/pts^2;
arts1_mean_err_nl = sum(sum(M.dist(u_orig,u_nl).^2))/pts^2;
%% Show Figures
if getDebugLevel('Figures') == 1
    figure;
    subplot(2,2,1)
    imagesc(u_0_SAR),colormap hsv, title('original');axis equal, axis tight;axis off
    subplot(2,2,3)
    imagesc(squeeze(atan2(u_final(2,:,:),u_final(1,:,:)))), colormap hsv, title('Final');axis equal, axis tight;axis off
    subplot(2,2,4)
    imagesc(squeeze(atan2(u_oracle(2,:,:),u_oracle(1,:,:)))),colormap hsv, title('Oracle');axis equal, axis tight;axis off
    subplot(2,2,2)
    imagesc(u_noisy_SAR),colormap hsv, title('noisy');axis equal, axis tight;axis off
    figure
    imagesc(squeeze(atan2(u_nl(2,:,:),u_nl(1,:,:))))
    title('NL Mean');colormap hsv;axis equal, axis tight;axis off
end
%% Write results
if getDebugLevel('WriteImages') == 1    
    map = colormap(hsv(256));
    imwrite(255*(squeeze(atan2(u_final(2,:,:),u_final(1,:,:)))+pi)/2/pi,map,[resultsFolder,'artInSAR_MMSE.png']);
    imwrite(255*(squeeze(atan2(u_oracle(2,:,:),u_oracle(1,:,:)))+pi)/2/pi,map,[resultsFolder,'artInSAR_MMSE_Oracle.png']);
    imwrite(255*(squeeze(atan2(u_nl(2,:,:),u_nl(1,:,:)))+pi)/2/pi,map,[resultsFolder,'artInSAR_MMSE_NL.png']);
    imwrite(255*(u_0_SAR+pi)/2/pi,map,[resultsFolder,'artInSAR_MMSE_Orig.png']);
    imwrite(255*(u_noisy_SAR+pi)/2/pi,map,[resultsFolder,'artInSAR_MSSE_Noisy.png']);
end
