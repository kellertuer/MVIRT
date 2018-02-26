%
% Denoising a synthetic SAR Image
% ===============================

% Demonstrating the Cyclic Proximal Point Algorithm for absolute Differences
% of first and second order on a 2D phase valued image.
%
% In InSAR data, the phase valued images are often noisy. Denoising them is
% a task performed on S1 valued data and demonstrated in this example
%
% This example is published in
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014).
%       Second order differences of cyclic data and applications in variational denoising.
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
%
% This file can be started without any changes; it initialized the Toolbox
% itself
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-04-06 | 2016-02-24
% see LICENSE.txt

% Init Toolbox
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')
% Global settings
useLogfile= true;
writeData = false;
loadData = true;
showFigures = true;
writeImages = true;

dataFolder = ['syntheticSARImage',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['syntheticSARImage',filesep];

N = 256;
% Algorithm parameters
iter = 4000;
epsilon = 0;
% Details for the data - these are ignored if data is loaded
sigma = 0.3;
name = ['S1-syntheticSARImageN',num2str(N)];
dataName = ['S1-syntheticSARImageN',num2str(N)];


%% Create Image
%
%
% Logfile?
if useLogfile
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
if writeData % Create new data
    [Zn,Z] = ArtificalSARImage(N,sigma);
    save([resultsFolder,dataName,'.mat'],'Zn','Z','sigma');
elseif loadData
    load([resultsFolder,dataName,'.mat'],'Zn','Z','sigma');
    metaData = dir([resultsFolder,dataName,'.mat']);
    disp(['Using File Data generated ',datestr(metaData.date),'.']);
    N = size(Zn,1);
else
    error('Either Loading or Creating(Wirting) Data must be set');
end

%% Initial Export of Starting Values
if showFigures
    figure(1); imagesc(Z+pi); title('Original Data'); colormap hsv;
    figure(2); imagesc(Zn+pi); title('noisy measurements'); colormap hsv;
    pause(0.05); %show all 3
end
if writeImages
    map = colormap(hsv(256));
    frs = 255*(Z+pi)/(2*pi);
    imwrite(frs,map,[resultsFolder,name,'wrapped.png'],'png');
    frs = 255*(Zn+pi)/(2*pi);
    imwrite(frs,map,[resultsFolder,name,'measured.png'],'png');
end

%% Parameters
alpha1s = [3/8, 0, 1/4];
alpha2s = [1/4,0,1/8];
beta1s = [0,1/8,1/8];
beta2s = [0,1/8,1/8];
gammas = [0,1/8,0];

clear problem
problem.M = S1;
problem.lambda = pi;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,iter,epsilon);
problem.f = permute(Zn,[3,1,2]); %the (unseen, third) manifold dimension has to be the first one

for i=1:length(alpha1s)
    problem.alpha = diag([alpha1s(i),alpha2s(i)]);
    problem.beta = [beta1s(i),gammas(i);0,beta2s(i)];
    problem.Debug = 1000;
    tic
    Zr = permute(CPP_AdditiveTV12(problem),[2,3,1]); %permute back for imaging
    toc
    ME = sum(sum(problem.M.dist(Zr,Z).^2))/length(Z(:));
    disp(['For ',num2str([problem.alpha(:).' problem.beta(:).']),...
        ' the error is ME:',sprintf('%3.7G',ME)]);
    if writeImages
        map = colormap(hsv(256));
        frs = 255*(Zr+pi)/(2*pi);
        imwrite(frs,map,[resultsFolder,name,'-P',num2str(i),'-denoised.png'],'png');
    end
    if showFigures
        figure(i+2); imagesc(Zr+pi); title(['Parameters ',...
        num2str([problem.alpha(:).' problem.beta(:).']),' ME:',sprintf('%3.7G',ME)]); colormap hsv;
    end
end
% End logfile
if useLogfile
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)
