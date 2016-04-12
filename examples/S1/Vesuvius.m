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
setDebugLevel('LevelMin',0);        % Minimal Value of all global Values
setDebugLevel('LevelMax',1000);     % Maximal Value ...
setDebugLevel('text',3);            % quite verbose text, but...
setDebugLevel('IterationStep',1000);% only every 100th iteration
setDebugLevel('WriteImages',1);     % 0: no file writing, 1: file writing
setDebugLevel('time',3);            % verbose time
setDebugLevel('Figures',1);         % 0: no figure display, 1: figures are displayed
setDebugLevel('logfile',1);         % 0: no logfile 1: logfile

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
if getDebugLevel('logfile')
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
if getDebugLevel('Figures')
    figure(1); imagesc(imgg+pi); title('Original Data'); colormap hsv;
    axis image; axis off
    pause(0.05); %show all 3
end
if getDebugLevel('WriteImages')
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

problem.MaxIterations = iter;
problem.Epsilon = epsilon;
problem.lambda=pi;
problem.M = S1;
problem.f = permute(imgg,[3,1,2]); %the (unseen, third) manifold dimension has to be the first one

for i=1:length(alpha1s)
    problem.alpha = [alpha1s(i),alpha2s(i) 0 0];
    problem.beta = [beta1s(i),beta2s(i),gammas(i)];
    imgr = permute(cppa_ad_2D(problem),[2,3,1]); %permute back for imaging
    if getDebugLevel('WriteImages')
        map = colormap(hsv(256));
        frs = 255*(imgr+pi)/(2*pi);
        imwrite(frs,map,[resultsFolder,name,'-P',num2str(i),'-denoised.png'],'png');
    end
    if getDebugLevel('Figures')
        figure(i+1); imagesc(imgr+pi); title(['Reconstruction',...
        num2str([problem.alpha problem.beta]),' .']); colormap hsv;
        axis image; axis off
    end
end
% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)