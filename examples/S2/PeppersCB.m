%
% Denoising the Peppers Image with an CB channel wise -TV
% ==========================================
%
% The peppers image is obstructed with sigma = 0.1 Gaussian noise and
% then denoised on several color space models.
% 
% Here we employ a channel wise TV model on the CB color space
%
% This Matlab script is used to generate Figure 5.3(e) and (f) in
%
% Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann, M. Bacak, G. Steidl, A. Weinmann: A second order non-smooth
%   variational model for restoring manifold-valued images, 2015.
% 
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-02-13 | 2015-05-05

start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
    run('../../initMVIRT.m')
    
%% Settings & Variables
writeImages = true;
loadData = true;
writeData = false;
showFigures = true;
useLogfile = true;

format compact

folder = ['examples',filesep,'S2',filesep];
results = ['Peppers',filesep];
dataFolder = ['..',filesep,'..',filesep,'data',filesep];
dataName = 'peppers';
name = 'DenoisePeppers_CB';

if useLogfile
    clc
    if exist([results,name,'.log'],'file')
        delete([results,name,'.log'])
    end
    diary([results,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
Vrgb = double(imread([dataFolder,dataName,'.jpg']))/255;
if loadData && exist([results,dataName,'-data.mat'],'file')
    load([dataFolder,dataName,'-data.mat']);
    metaData = dir([results,dataName,'-data.mat']);
    disp(['Using File Data generated ',datestr(metaData.date),'.']);
else
    %Create Data
    sigma = 0.1; %Paper Grohs
    Vrgbn = max(min(Vrgb + sigma*randn(size(Vrgb)),1),0); %How does Grohs handle this?
    disp(['Using new Data generated ',datestr(datetime),'.']);
    if writeData %Write this version to file
        save([dataFolder,dataName,'-data.mat'],'Vrgbn','sigma')
    end
end
    
if showFigures
    figure(1); imagesc(Vrgb); axis image; axis off; title('Original');
    figure(2); imagesc(Vrgbn); axis image; axis off; title(['noisy, \sigma=',num2str(sigma),'.']);
end
if writeImages
    imwrite(Vrgb,[results,name,'-original.png'],'png');
    imwrite(Vrgbn,[results,name,'-noisy-sigma',num2str(sigma),'.png'],'png');
end

% VB is on S2, VC is in R
VB = sqrt(sum(Vrgbn.^2,3));
VC = Vrgbn./repmat(VB,[1,1,3]);
% Zero -> Northpole
VC1 = VC(:,:,1); VC1(VB==0) = 0; VC(:,:,1) = VC1;
VC2 = VC(:,:,2); VC2(VB==0) = 0; VC(:,:,2) = VC2;
VC3 = VC(:,:,3); VC3(VB==0) = 1; VC(:,:,3) = VC3;

%% Parameters
alphas = 0:1/4:1;
betas = 0:1/4:1;

problem.lambda = pi/2;
problem.Debug = 100;
problemR = problem;
%
problem.M = Sn(2);
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,400,0);
problem.f = permute(VC,[3,1,2]); 
%
problemR.M = Rn(1);
problemR.stoppingCriterion = stopCritMaxIterEpsilonCreator(problemR.M,800,0);
problemR.f = permute(VB,[3,1,2]);

% Search for best
bestPSNR = 0;
bestParams = [0,0];

bestPSNRTV1 = 0;
bestParamTV1 = 0;

if useLogfile %display parameters
    problem
    problem.M
end

for beta = betas
    for alpha = alphas
        [alpha beta]
        if (alpha+beta)>0 % not both zero
            problem.alpha = alpha*[1,1,0,0]/2;
            problem.beta= beta/2;
            tic;
            RC = CPP_AdditiveTV12(problem);
            %
            problemR.alpha = alpha*[1,1,0,0];
            problemR.beta= beta;
            RB = CPP_AdditiveTV12(problemR); %using same values for lambda alpha and beta
            toc
            %
            Rrgb = permute(RC,[2,3,1]).*repmat(permute(RB,[2,3,1]),[1,1,3]);

            psnrV = psnr(Rrgb,Vrgb);
            if psnrV > bestPSNR
                bestPSNR = psnrV;
                bestParams = [alpha,beta];
                disp(['New optimal PSNR: ',num2str(psnrV),'.']);
            else
                disp(['PSNR: ',num2str(psnrV),'.']);
            end
            if (beta==0) && (psnrV > bestPSNRTV1)
                bestPSNRTV1 = psnrV;
                bestParamTV1 = alpha;
                disp(['New optimal PSNR for TV1: ',num2str(psnrV),'.']);
            end
            if showFigures
                figure(3); imagesc(Rrgb); axis image; axis off; title(['Denoised, \alpha=',num2str(alpha),' \beta=',num2str(beta),'.']);
                pause(0.02);
            end
            if writeImages
                imwrite(Rrgb,[results,name,'-denoised-alpha',num2str(alpha),'-beta',num2str(beta),'.png'],'png');
            end
        end
    end
end
disp(['Best PSNR ',num2str(bestPSNR),' obtained for ',num2str(bestParams),'.']); 
disp(['Best PSNR for TV1 ',num2str(bestPSNRTV1),' obtained for ',num2str(bestParamTV1),'.']); 
if useLogfile
    diary off;
end
cd(start)