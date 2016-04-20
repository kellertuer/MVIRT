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
setDebugLevel('LevelMin',0);        % Minimal Value of all global Values
setDebugLevel('LevelMax',1000);     % Maximal Value ...
setDebugLevel('text',3);            % quite verbose text, but...
setDebugLevel('IterationStep',100); % only every Xth iteration
setDebugLevel('WriteImages',1);     % 0: no file writing, 1: file writing
setDebugLevel('LoadData',1);        % 0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0);       % 0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('time',3);            % verbose time
setDebugLevel('Figures',1);         % 0: no figure display, 1: figures are displayed
setDebugLevel('logfile',1);         % 0: no logfile 1: logfile
format compact

folder = ['examples',filesep,'S2',filesep];
results = ['Peppers',filesep];
dataFolder = ['..',filesep,'..',filesep,'data',filesep];
dataName = 'peppers';
name = 'DenoisePeppers_CB';

if getDebugLevel('logfile')
    clc
    if exist([results,name,'.log'],'file')
        delete([results,name,'.log'])
    end
    diary([results,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
Vrgb = double(imread([dataFolder,dataName,'.jpg']))/255;
if getDebugLevel('LoadData') && exist([results,dataName,'-data.mat'],'file')
    load([results,dataName,'-data.mat']);
    metaData = dir([results,dataName,'-data.mat']);
    debug('text',2,'Text',['Using File Data generated ',datestr(metaData.date),'.']);
else
    %Create Data
    sigma = 0.1; %Paper Grohs
    Vrgbn = max(min(Vrgb + sigma*randn(size(Vrgb)),1),0); %How does Grohs handle this?
    debug('text',2,'Text',['Using new Data generated ',datestr(datetime),'.']);
    if getDebugLevel('WriteData')>0 %Write this version to file
        save([results,dataName,'-data.mat'],'Vrgbn','sigma')
    end
end
    
if getDebugLevel('figures')
    figure(1); imagesc(Vrgb); axis image; axis off; title('Original');
    figure(2); imagesc(Vrgbn); axis image; axis off; title(['noisy, \sigma=',num2str(sigma),'.']);
end
if getDebugLevel('WriteImageFiles')
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
alphas = [0 0.0938];
betas = [0 0.469];

problem.Epsilon = 0;
problem.MaxIterations = 400;
problem.tau=1;
problem.lambda = pi/2;
problemR = problem;
%
problem.M = Sn(2);
problem.M.steps = 5;
problem.f = permute(VC,[3,1,2]); 
%
problemR.M = S1mRn(0,1);
problemR.f = permute(VB,[3,1,2]);

% Search for best
bestPSNR = 0;
bestParams = [0,0];

bestPSNRTV1 = 0;
bestParamTV1 = 0;

if getDebugLevel('logfile') %display parameters
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
            RC = cppa_ad_2D(problem);
            %
            problemR.alpha = alpha*[1,1,0,0];
            problemR.beta= beta;
            RB = cppa_ad_2D(problemR); %using same values for lambda alpha and beta
            toc
            %
            Rrgb = permute(RC,[2,3,1]).*repmat(permute(RB,[2,3,1]),[1,1,3]);

            psnrV = psnr(Rrgb,Vrgb);
            if psnrV > bestPSNR
                bestPSNR = psnrV;
                bestParams = [alpha,beta];
                disp(['New optimal PSNR:',num2str(psnrV),'.']);
            else
                disp(['PSNR:',num2str(psnrV),'.']);
            end
            if (beta==0) && (psnrV > bestPSNRTV1)
                bestPSNRTV1 = psnrV;
                bestParamTV1 = alpha;
                disp(['New optimal PSNR for TV1:',num2str(psnrV),'.']);
            end
            if getDebugLevel('figures')>0
                figure(3); imagesc(Rrgb); axis image; axis off; title(['Denoised, \alpha=',num2str(alpha),' \beta=',num2str(beta),'.']);
                pause(0.02);
            end
            if getDebugLevel('WriteImageFiles')>0
                imwrite(Rrgb,[results,name,'-denoised-alpha',num2str(alpha),'-beta',num2str(beta),'.png'],'png');
            end
        end
    end
end
disp(['Best PSNR',num2str(bestPSNR),' obtained for ',num2str(bestParams),'.']); 
disp(['Best PSNR for TV1 ',num2str(bestPSNRTV1),' obtained for ',num2str(bestParamTV1),'.']); 
if getDebugLevel('logfile');
    diary off;
end
cd(start)