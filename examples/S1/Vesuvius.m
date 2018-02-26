%
% Denoising the Vesuvius InSAR Data
% ===============================
%
% In InSAR data, the phase-valued images are often noisy. Denoising them is
% a task performed on S1 valued data and demonstrated in this example.
%
% The Mount Vesuvius Data is obtained from
% https://earth.esa.int/workshops/ers97/program-details/speeches/rocca-et-al/
%
% This example is published in
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
%
% This file can be started without any changes; it initialized the Toolbox
% itself
% ---
% Manifold-valued Image Restoration 1.0
% R. Bergmann ~ 2014-04-06 | 2016-02-24
% see LICENSE.txt

% Init Toolbox
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')
% Global settings

writeImages = true;
showFigures = true;
useLogfile  = true;

dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['Vesuvius',filesep];

% Algorithm parameters
iter = 4000;
epsilon = 0;
% Details for the data - these are ignored if data is loaded
name = 'S1-Vesuvius';
dataName = 'Vesuvius';


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

img = imread([dataFolder,dataName,'.png']);
imgg = (double(rgb2gray(img))/256)*(2*pi)-pi;

%% Initial Export of Starting Values
if showFigures
    figure(1); imagesc(imgg+pi); title('Original Data'); colormap hsv;
    axis image; axis off
    pause(0.05); %show all 3
end
if writeImages
    map = colormap(hsv(256));
    frs = 255*(imgg+pi)/(2*pi);
    imwrite(frs,map,[resultsFolder,name,'-orig.png'],'png');
end

%% Parameters
alpha1s = 1/4;
alpha2s = 1/4;
beta1s = 3/4;
beta2s = 3/4;
gammas = 3/4;

problem.M = S1();
problem.lambda = pi;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,iter,epsilon);
problem.f = permute(imgg,[3,1,2]);
problem.Debug = 1000;

for i=1:length(alpha1s)
    problem.alpha = diag( [alpha1s(i),alpha2s(i)] );
    problem.beta = [beta1s(i),gammas(i);0,beta2s(i)];
    tic
        imgr = permute(CPP_AdditiveTV12(problem),[2,3,1]); %permute back for imaging
    toc
        
    if writeImages
        map = colormap(hsv(256));
        frs = 255*(imgr+pi)/(2*pi);
        imwrite(frs,map,[resultsFolder,name,'-P',num2str(i),'-denoised.png'],'png');
    end
    if showFigures
        figure(i+1); imagesc(imgr+pi); title(['Reconstruction',...
        num2str([problem.alpha(:).', problem.beta(:).']),' .']); colormap hsv;
        axis image; axis off
    end
end
% End logfile
if useLogfile
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)