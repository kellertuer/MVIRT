%
% Denoising a noisy color channel in the HSV model of a real world image
% ======================================================================
%
% An image might be abstructed by noise on the color channel. Denoising
% such an image does not work on RGB color space, but on the HSV space it
% does, though it requires phase-valued denoising.
%
% The image ?sailboat on lake? is taken from the USC-SIPI Image Database,
% namely http://sipi.usc.edu/database/database php?volume=misc&image=14.
%
% This example is in a similar version (less noise, less iterations)
% published in
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
%
% This file can be started without any changes; it initialized the Toolbox
% itself
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-06-20 | 2016-02-24
% see LICENSE.txt
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')
% Global settings
setDebugLevel('LevelMin',0);        % Minimal Value of all global Values
setDebugLevel('LevelMax',100);      % Maximal Value ...
setDebugLevel('text',3);            % quite verbose text, but...
setDebugLevel('IterationStep',100); % only every 100th iteration
setDebugLevel('WriteImages',1);     % 0: no file writing, 1: file writing
setDebugLevel('time',3);            % verbose time
setDebugLevel('LoadData',1);        % 0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0);       % 0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('Figures',1);         % 0: no figure display, 1: figures are displayed
setDebugLevel('logfile',1);         % 0: no logfile 1: logfile

dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['SailboatHueDenoising',filesep];

% Algorithm parameters
iter = 4000;
epsilon = 0;
% Details for the data - these are ignored if data is loaded
sigma = 0.2;
name = 'sailboat on lake';
dataName = 'S1-sailboatHueDenoising';
%% Create Image
%
%
% use a Logfile?
if getDebugLevel('logfile')
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
if getDebugLevel('WriteData') % Create new data
    img = double(imread([dataFolder,name,'.tiff']))/255;
    Z = rgb2hsv(img);
    Zn = Z;
    Zn(:,:,1) = symMod(2*pi*Z(:,:,1) - pi + 2*pi*sigma*randn(size(Z(:,:,1))),2*pi); %noisy signal
    save([resultsFolder,dataName,'.mat'],'img','Zn','Z','sigma');
elseif getDebugLevel('LoadData')
    load([resultsFolder,dataName,'.mat'],'img','Zn','Z','sigma');
    metaData = dir([resultsFolder,dataName,'.mat']);
    debug('text',3,'Text',['Using File Data generated ',datestr(metaData.date),'.']);
else
    error('Either Loading or Creating (Writing) Data must be set');
end
numel = size(Zn,1)*size(Zn,2);
%
% image output
if getDebugLevel('Figures')
    figure(1); imagesc(img); title('Original image ''sailboat on lake''');
    axis image; axis off;
    imgn = Zn;
    imgn(:,:,1) = mod(imgn(:,:,1)+pi,2*pi)/(2*pi);
    figure(2); imagesc(hsv2rgb(imgn)); title('''sailboat on lake'' with noisy hue');
    axis image; axis off;
    pause(0.01);
end
if getDebugLevel('WriteImages')
    imgn = Zn;
    imgn(:,:,1) = mod(imgn(:,:,1)+pi,2*pi)/(2*pi);
    imwrite(img,[resultsFolder,name,'-orig.png'],'png');
    imwrite(hsv2rgb(imgn),[resultsFolder,name,'-noisy.png'],'png');
end
% Number 1: H phase valued denoising
clear problem
problem.MaxIterations = iter; problem.Epsilon = epsilon;
problem.f = permute(Zn(:,:,1),[3,1,2]); problem.M = S1;
problem.alpha = 2*pi*1/4*[1,1,0,0];
problem.beta = 2*pi*[1,1,0];
problem.lambda=pi;
ZnH = Zn;
ZnH(:,:,1) = mod( permute(cppa_ad_2D(problem),[2,3,1])+pi,2*pi )/(2*pi);
% 2: Same as real valued
problem.M = S1mRn(0,1);
ZnR = Zn;
ZnR(:,:,1) = mod( permute(cppa_ad_2D(problem),[2,3,1])+pi,2*pi )/(2*pi);

if getDebugLevel('Figures')
    figure(3); imagesc(hsv2rgb(ZnH)); title('Denoised Hue with S1-valued CPPA');
    axis image; axis off;
    figure(4); imagesc(hsv2rgb(ZnR)); title('Denoised Hue with R-valued CPPA');
    axis image; axis off;
    pause(0.01)
end
if getDebugLevel('WriteImages')
    imwrite(hsv2rgb(ZnH),[resultsFolder,name,'-S1_Denoised.png'],'png');
    imwrite(hsv2rgb(ZnR),[resultsFolder,name,'-R-Denoised.png'],'png');
end

%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)