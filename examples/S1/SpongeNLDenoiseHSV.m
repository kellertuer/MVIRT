%%       Denoising of an Hue valued image
%
%   Denoising of the hue of a color image, with the NL MMSE method.
%
%
%   This example is published in Sect 5.3 of
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%     Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
%
%  This file can be started without any changes; it initializes the Toolbox
%  itself
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
%  J. Persch  ~ 2017-01-06 | R. Bergmann 2016-01-07
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
%
resultsFolder = ['Sponge',filesep];
% Parameter for noise
sigma = 0.6;
% Load image and transform it into HSV colorspace
M = Sn(1);
sponges = double(imread('sponges.png'))/255;
transformed_corals = rgb2hsv(sponges);
h = 2*pi*transformed_corals(:,:,1)-pi;
u_noisy_Sn1 = zeros(2,size(h,1),size(h,2));
u_orig = zeros(2,size(h,1),size(h,2));
u_orig(1,:,:) = cos(h);
u_orig(2,:,:) = sin(h);
if ~(getDebugLevel('LoadData')==1)
noise = sigma*randn(size(h));
u_noisy_Sn1(1,:,:) = cos(h+noise);
u_noisy_Sn1(2,:,:) = sin(h+noise);
else
    load([resultsFolder,'NoisySponge.mat'],'u_noisy_Sn1');
end
% Initialize problem
clear problem
problem.M = M;
problem.f = u_noisy_Sn1;
problem.sigma =sigma;
problem.gamma = 1;

problem.patch_size_1 = 7;
problem.window_size_1 = 81;

problem.patch_size_2 = 7;
problem.window_size_2 = 81;

problem.K_1 = 70;
problem.K_2 = 70;
% TV parameters
problem.lambda = pi/2;
problem.alpha = 0.45;
problem.beta = 0;
%% Run Algorithm
[u_final,u_oracle] = NL_MMSE_2D(problem);
u_tv = cppa_ad_2D(problem);
%% get Images
noisy = transformed_corals;
noisy(:,:,1) = mod(squeeze(atan2(u_noisy_Sn1(2,:,:),u_noisy_Sn1(1,:,:)))/2/pi+0.5,1);
sponge_s1psnr_noisy =  psnr(hsv2rgb(noisy),sponges);
sponge_s1mse_noisy = sum(sum(M.dist(u_noisy_Sn1,u_orig).^2))/(numel(sponges)/3);
denoised_oracle = transformed_corals;
denoised_oracle(:,:,1) = mod(squeeze(atan2(u_oracle(2,:,:),u_oracle(1,:,:)))/2/pi+0.5,1);
sponge_s1psnr_oracle = psnr(hsv2rgb(denoised_oracle),sponges);
sponge_s1mse_oracle = sum(sum(M.dist(u_oracle,u_orig).^2))/(numel(sponges)/3);
denoised = transformed_corals;
denoised(:,:,1) = mod(squeeze(atan2(u_final(2,:,:),u_final(1,:,:)))/2/pi+0.5,1);
sponge_s1psnr_final = psnr(hsv2rgb(denoised),sponges);
sponge_s1mse_final = sum(sum(M.dist(u_final,u_orig).^2))/(numel(sponges)/3);
tv_denoised = transformed_corals;
tv_denoised(:,:,1) = mod(squeeze(atan2(u_tv(2,:,:),u_tv(1,:,:)))/2/pi+0.5,1);
sponge_s1psnr_tv = psnr(hsv2rgb(tv_denoised),sponges);
sponge_s1mse_tv = sum(sum(M.dist(u_tv,u_orig).^2))/(numel(sponges)/3);
%% Show Images
if getDebugLevel('Figures')
figure(1);
imagesc(sponges)
figure
imagesc(hsv2rgb(noisy));
figure
imagesc(hsv2rgb(denoised))
title('Final HSV')
figure
imagesc(hsv2rgb(tv_denoised))
title('TV HSV')
end
%% Export Images
if getDebugLevel('WriteImages')
   imwrite(sponges,[resultsFolder,'orig_corals.png']);
   imwrite(hsv2rgb(noisy),[resultsFolder,'noisy_corals.png']);
   imwrite(hsv2rgb(denoised),[resultsFolder,'denoised_corals.png']);
   imwrite(hsv2rgb(tv_denoised),[resultsFolder,'tv_denoised_corals.png']);
   imwrite(uint8(transformed_corals(:,:,1)*255),hsv(256),[resultsFolder,'hue.png']);
   imwrite(uint8(noisy(:,:,1)*255),hsv(256),[resultsFolder,'hue_noisy.png']);
   imwrite(uint8(denoised(:,:,1)*255),hsv(256),[resultsFolder,'hue_denoised_corals.png']);
   imwrite(uint8(tv_denoised(:,:,1)*255),hsv(256),[resultsFolder,'hue_tv_denoised_corals.png']);
end