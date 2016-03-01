%
% Denoising and Inpainting the Vesuvius InSAR Data
% ===============================
%
% In InSAR data, the phase valued images are often noisy and sometimes lossy.
% Inpainting and Denoising them is a task performed on S1 valued data and
% demonstrated in this example.
%
% The Mount Vesuvius Data is obtained from
% https://earth.esa.int/workshops/ers97/program-details/speeches/rocca-et-al/
%
% This example is published in
%    R. Bergmann, A. Weinmann: Inpainting of Cyclic Data Using First and
%      Second Order Differences Second Order Differences, in: EMMCVPR 2015.
%
% This file can be started without any changes; it initialized the Toolbox
% itself
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-10-12 | 2016-02-24
% see LICENSE.txt

% Init Toolbox
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')

%% Settings
setDebugLevel('LevelMin',0);        % Minimal Value of all global Values
setDebugLevel('LevelMax',1000);     % Maximal Value ...
setDebugLevel('text',3);            % quite verbose text, but...
setDebugLevel('IterationStep',1000);% only every 100th iteration
setDebugLevel('WriteImages',1);     % 0: no file writing, 1: file writing
setDebugLevel('time',3);            % verbose time
setDebugLevel('LoadData',1);        % 0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0);       % 0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('Figures',1);         % 0: no figure display, 1: figures are displayed
setDebugLevel('logfile',1);         % 0: no logfile 1: logfile

%
% ------
dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['Vesuvius',filesep];

% Details for the data - these are ignored if data is loaded
name = 'S1-Vesuvius-Inpaint';
dataName = 'Vesuvius';
lossy = .2; %loose 1/5th of data.

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
% ------
% Data generation / Loading
%
% Segments, similar to Strekalovskyi et al., but rotated here and with variable
% number of segments
%% Data generation (Just run this one once to be comparable
if getDebugLevel('WriteData') % Create new data
    imgl = imread([dataFolder,dataName,'.png']);
    V = (double(rgb2gray(imgl))/256)*(2*pi)-pi;
    Mask = (randi(ceil(1/lossy),size(V))~=1);
    Vl = V.*Mask; %Just loose data
    save([resultsFolder,dataName,'.mat'],'V','lossy','Vl','Mask');
elseif getDebugLevel('LoadData')
    load([resultsFolder,dataName,'.mat'],'V','lossy','Vl','Mask');
    metaData = dir([resultsFolder,dataName,'.mat']);
    debug('text',3,'Text',['Using File Data generated ',datestr(metaData.date),'.']);
else
    error('Either Loading or Creating (Writing) Data must be set');
end
% ------
if getDebugLevel('WriteImages')
    Vexp = uint8((Vl+pi)/(2*pi)*255);
    % white mask of cut out (rest transparent
    imwrite(uint8(255*(1-Mask)),[resultsFolder,name,'-mask-w.png'],'Alpha',255*(1-Mask));
    % black mask of cut out (rest transparent
    imwrite(uint8(255*Mask),[resultsFolder,name,'-mask-b.png'],'Alpha',255*(1-Mask));
    % Display mask
    pr = 100 * sum(sum(Mask==0))/length(Mask(:));
end
if getDebugLevel('Figures')
    figure(1); imagesc(V,[-pi,pi]); colormap(hsv(1024));
    title('Original Data V');
    axis image; axis off
    figure(2); imagesc(Mask,[0,1]); colormap(gray);
    axis image; axis off
    title(['Mask: White pixels are kept, black ones are destroyed. (',num2str(pr,3),'% destroyed).']);
    Vl = V.*Mask; %Just loose data
    figure(3); imagesc(Vl,[-pi,pi]);  colormap(hsv(1024));
    axis image; axis off
    title('input data with destroyed points');
end

%% Parameters for the CPPA
clear problem
problem.alpha = 1/4*[1,1,0,0];
problem.beta=3/4*[1,1,1];
problem.lambda = pi;
problem.f = permute(Vl,[3,1,2]);
problem.MaxIterations = 4000;
problem.Epsilon = 10^(-9);
problem.M = S1();
problem.UnknownMask = ~Mask;
%
% Datamask (1=existent, 0 nonexistent data)
% RegMask (1=pixel affected by surrounding TV terms, 0 not affected (just
% data applied if that's 1)
Vres = permute(cppa_ad_2D(problem),[2,3,1]);
if getDebugLevel('Figures')
    figure(4); imagesc(Vres,[-pi,pi]); colormap(hsv(1024));
    axis image; axis off
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=',num2str(problem.beta(1),4),' on R.']);
end
if getDebugLevel('WriteImages')
    Vresexp = uint8((Vres+pi)/(2*pi)*255);
    % Export result
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impainted-',num2str(problem.alpha(1),4),'-',num2str(problem.beta(1),4),'.png'],'png');
    % Compare to usual denoise
end
problem.alpha = problem.alpha/2;
problem.alpha = problem.beta/2;
problem.f = permute(V,[3,1,2]);
Vcomp = permute(cppa_ad_2D(problem),[2,3,1]);
VDist = problem.M.dist(Vcomp,Vres);
%
if getDebugLevel('Figures')
    figure(5); imagesc(VDist); colormap(flipud(gray(1024)));
    axis image; axis off
    ME = sum(sum(VDist.^2))/length(V(:));
    title(['Reconstruction Error. MSE = ',num2str(ME),'.']);
end
if getDebugLevel('WriteImages')
    VDistexp = uint8((VDist-min(min(VDist)))/(max(max(VDist))-min(min(VDist)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-error-',num2str(problem.alpha(1),4),'-',num2str(problem.beta(1),4),'.png'],'png');
end
if getDebugLevel('Figures')
    figure(6); imagesc(Vcomp); colormap hsv;
    axis image; axis off
    title('Denoised without missing data');
end
if getDebugLevel('WriteImages')
    Vcompexp = uint8((Vcomp+pi)/(2*pi)*255);
    imwrite(Vcompexp,colormap(hsv(255)),[resultsFolder,name,'-denoised-',num2str(problem.alpha(1),4),'-',num2str(problem.beta(1),4),'.png'],'png');
end
%
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)