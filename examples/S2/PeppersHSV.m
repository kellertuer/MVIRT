%
% Denoising the Peppers Image with an HSV-TV
% ==========================================
%
% The peppers image is obstructed with sigma = 0.1 Gaussian noise and
% then denoised on several color space models.
% 
% Here we employ a vectorial TV model on the HSV color space
%
% This Matlab script is used to generate Figure 5.3(c) in
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann,
%     A second order non-smooth variational model for restoring
%     manifold-valued images,
%     SIAM Journal on Scientific Computing, 38, (1), A567?A597, 2016.
% 
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-02-13 | 2016-03-01
% see LICENSE.txt

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
name = 'DenoisePeppers_HSV1';

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
    debug('text',2,'Text',['Using _new_ Data generated ',datestr(datetime),'.']);
    % never write here
    if getDebugLevel('WriteData')>0 %Write this version to file
        save([results,dataName,'-data.mat'],'Vrgbn','sigma')
    end
end
    
if getDebugLevel('Figures')
    figure(1); imagesc(Vrgb); axis image; axis off; title('Original');
    figure(2); imagesc(Vrgbn); axis image; axis off; title(['noisy, \sigma=',num2str(sigma),'.']);
end
if getDebugLevel('WriteImages')
    imwrite(Vrgb,[results,name,'-original.png'],'png');
    imwrite(Vrgbn,[results,name,'-noisy-sigma',num2str(sigma),'.png'],'png');
end

Vhsvn = rgb2hsv(Vrgbn);
%Scale to 2 pi
Vhsvs = 2*pi*Vhsvn;
Vhsvs(:,:,1) = symMod(Vhsvs(:,:,1),2 * pi);

%% Parameters
alphas = [1/16 3/32] * 2*pi;
betas = [0 1/32] * 2*pi;
problem.Epsilon = 10^(-9);
problem.MaxIterations = 400;
problem.lambda = pi/2;
%
problem.M = S1mRn(1,2);
problem.f = permute(Vhsvs,[3,1,2]); 

% Search for best

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
            problem.alpha = alpha;
            problem.beta= beta;
            tic;
            Vhsvr = permute(cppa_ad_2D(problem),[2,3,1])/(2*pi);
            Vhsvr(:,:,1) = mod(Vhsvr(:,:,1),1); %back from (-.5,.5) to (0,1)
            toc
            %
            Vrgbr = hsv2rgb(Vhsvr);
            psnrV = psnr(Vrgb,Vrgbr);
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
            if getDebugLevel('Figures')>0
                figure(3); imagesc(Vrgbr); axis image; axis off; title(['Denoised, \alpha=',num2str(alpha),' \beta=',num2str(beta),'.']);
                pause(0.02);
            end
            if getDebugLevel('WriteImages')>0
                imwrite(Vrgbr,[results,name,'-denoised-alpha',num2str(alpha),'-beta',num2str(beta),'.png'],'png');
            end
        end
    end
end
disp(['Best PSNR',num2str(bestPSNR),' obtained for ',num2str(bestParams),'.']); 
if bestPSNRTV1 > 0
    disp(['Best PSNR for TV1 ',num2str(bestPSNRTV1),' obtained for ',num2str(bestParamTV1),'.']); 
end
if getDebugLevel('logfile');
    diary off;
end
cd(start)