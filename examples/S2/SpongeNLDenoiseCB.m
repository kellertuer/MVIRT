%%       Denoising of the Chromaticity
%
%   Denoising of an the chromaticity channel of an RGB image of a sponge
%
%
%   This example is published in Sect 5.3 of
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%       Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
%
%
%  This file can be started without any changes; it initializes the Toolbox
%  itself
%
% Manifold-valued Image Restoration Toolbox 1.2
%  J. Persch  ~ 2016-07-05 | R. Bergmann ~ 2017-01-24
% see LICENSE.txt
run('../../initMVIRT.m')
%% Debug
setDebugLevel('LevelMin',0);
setDebugLevel('LevelMax',100);
setDebugLevel('text',3);
setDebugLevel('time',3);
setDebugLevel('LoadData',1);
setDebugLevel('WriteImages',1);
setDebugLevel('Figures',true);

resultsFolder = ['Sponge',filesep];
% Parameters for noise and import image
sigma = 0.2;
M = Sn(2);
sponge = double(imread('sponges.png'))/255;
B = sqrt(sum(sponge.^2,3));
C = sponge./repmat(B,1,1,3);
if getDebugLevel('LoadData')
    C_noisy_S2 = M.addNoise(permute(C,[3,1,2]),sigma);
else
    load([resultsFolder,'S2Sponge.mat'],'u_noisy_Sn2') ;
    C_noisy_S2 = u_noisy_Sn2;
end

u_orig = permute(C,[3,1,2]);
% Initialize problem
clear problem
problem.M = M;
problem.f = C_noisy_S2;
problem.sigma =sigma;
problem.gamma = 1;

problem.patch_size_1 = 5;
problem.window_size_1 = 37;

problem.patch_size_2 = 5;
problem.window_size_2 = 37;

problem.K_1 = 110;
problem.K_2 = 110;
% TV Parameter
problem.lambda = pi/2;
problem.alpha = 0.21;
problem.beta = 0;
%% NL MMSE Denoising
[u_final,u_oracle] = NL_MMSE_2D(problem);
u_tv = cppa_ad_2D(problem);
%% Get Images
noisy_corals = permute(C_noisy_S2,[2,3,1]).*repmat(B,1,1,3);
sponge_s2psnr_noisy = psnr(noisy_corals,sponge);
sponge_s2mse_noisy = sum(sum(M.dist(C_noisy_S2,u_orig).^2))/(numel(sponge)/3);
denoised_corals = permute(u_final,[2,3,1]).*repmat(B,1,1,3);
sponge_s2psnr_final = psnr(denoised_corals,sponge);
sponge_s2mse_final = sum(sum(M.dist(u_final,u_orig).^2))/(numel(sponge)/3);
oracle_corals = permute(u_oracle,[2,3,1]).*repmat(B,1,1,3);
sponge_s2psnr_oracle = psnr(oracle_corals,sponge);
sponge_s2mse_oracle = sum(sum(M.dist(u_oracle,u_orig).^2))/(numel(sponge)/3);
tv_corals = permute(u_tv,[2,3,1]).*repmat(B,1,1,3);
sponge_s2psnr_tv = psnr(tv_corals,sponge);
sponge_s2mse_tv = sum(sum(M.dist(u_tv,u_orig).^2))/(numel(sponge)/3);
%% Show Results
if getDebugLevel('Figures') == 1
    figure
    imagesc(sponge)
    title('Original sponges');
    figure
    imagesc(noisy_corals)
    title(['Noisy Corals: ',num2str(sponge_s2psnr_noisy)]);
    figure
    imagesc(denoised_corals)
    title(['Denoised Corals Final: ',num2str(sponge_s2psnr_final)]);
    figure
    imagesc(oracle_corals)
    title(['Denoised Corals Oracle: ',num2str(sponge_s2psnr_oracle)]);
    figure
    imagesc(tv_corals)
    title(['Denoised Corals TV: ',num2str(sponge_s2psnr_tv)]);
end
%% Export Results
if getDebugLevel('WriteImages') == 1
    imwrite(sponge,[resultsFolder,'S2corals_orig.png']);
    imwrite(noisy_corals,[resultsFolder,'S2corals_noisy.png']);
    imwrite(denoised_corals,[resultsFolder,'S2corals_ora.png']);
    imwrite(oracle_corals,[resultsFolder,'S2corals_denoised.png']);    
    imwrite(tv_corals,[resultsFolder,'S2corals_tv.png']);
    imwrite(C,[resultsFolder,'chromaticity.png']);
    imwrite(permute(C_noisy_S2,[2,3,1]),[resultsFolder,'chromaticity_noisy.png']);
    imwrite(permute(u_final,[2,3,1]),[resultsFolder,'chromaticity_denoised.png']);    
    imwrite(permute(u_tv,[2,3,1]),[resultsFolder,'chromaticity_tv.png']);
end