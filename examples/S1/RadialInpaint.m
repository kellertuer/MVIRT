%
% Inpainting an artificial example first and second oder differences.
% =======================================================================
%
% Demonstrating the Cyclic Proximal Point Algorithm for absolute Differences
% of first and second order on a phase valued radial image and alost hole
% in the middle.
%
% The comparison is published as the third numerical example (Sec. 5) of
%
%    R. Bergmann, A. Weinmann: Inpainting of Cyclic Data Using First and
%      Second Order Differences Second Order Differences, in: EMMCVPR 2015.
%
% This file can be started without any changes; it initializes the Toolbox
% itself
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-10-05

% Init Toolbox
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')
%
%
%% Settings
setDebugLevel('LevelMin',0);
setDebugLevel('LevelMax',1000);
setDebugLevel('text',3); %verbose, but...
setDebugLevel('IterationStep',10000); %only every 50th iteration
setDebugLevel('WriteImages',1); %0: no file writing, 1: file writing
setDebugLevel('time',3); %many time measurements
setDebugLevel('Figures',1); %0: no figure display, 1: figures are displayed (disable e.g. for cluster/console work)
setDebugLevel('logfile',1); %0: no logfile 1: logfile

dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['RadialInpaint',filesep];
name = 'PhaseRadial';

%% Logfile?
if getDebugLevel('logfile')
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end

pts = 128;
[X,Y] = meshgrid(-.5:1/(pts-1):.5,-.5:1/(pts-1):.5);
% Linear segment:
V = atan2(Y,X);

if getDebugLevel('Figures')
    figure(1); imagesc(V,[-pi,pi]); colormap(hsv(1024));
    title(['Original Data ''',name,'''.']);
    axis image
    axis off
end
if getDebugLevel('WriteImages')
    Vexp = uint8((V+pi)/(2*pi)*255);
    %Export original image
    imwrite(Vexp,colormap(hsv(255)),[resultsFolder,name,'.png'],'png');
end
%% ------
% Mask generation
% Destroy two thirds of all rows/cols -> just 1/9 kept
% here, 1 means keeping data, 0 means destroyd
M = (Y.^2+X.^2) >= (0.255)^2;
img = V.*M; %Just loose data
if getDebugLevel('Figures')
    figure(2); imagesc(M,[0,1]); colormap(gray); axis square
    pr = 100 * sum(sum(M==0))/numel(M);
    title(['Mask: White pixels are kept, black ones are destroyed. (',num2str(pr,3),'% destroyed).']);
    figure(3); imagesc(img,[-pi,pi]);  colormap(hsv(1024)); axis square
    title('input data with destroyed points');
end
if getDebugLevel('WriteImages')
    % white mask of cut out (rest transparent
    imwrite(uint8(255*(1-M)),[resultsFolder,name,'-mask-w.png'],'Alpha',255*(1-M));
    % black mask of cut out (rest transparent
    imwrite(uint8(255*M),[resultsFolder,name,'-mask-b.png'],'Alpha',255*(1-M));
end
%% Parameters for the cyclic proximal point algorithm

problem.alpha = 1/2*[1,1,1,1];
problem.beta=1/4*[1,1,1];
problem.lambda = pi/2;
problem.f = permute(img,[3,1,2]);
problem.MaxIterations = 4000;
problem.Epsilon = 10^(-9);
problem.M = S1();
problem.UnknownMask = ~M;
problem.RegMask = ~M;
VresTV2 = permute(cppa_ad_2D(problem),[2,3,1]);
if getDebugLevel('Figures')
    figure(4); imagesc(VresTV2,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=',num2str(problem.beta(1),4),'.']);
end
problem.beta = [0,0,0];
VresTV = permute(cppa_ad_2D(problem),[2,3,1]);
if getDebugLevel('Figures')
    figure(5); imagesc(VresTV,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=0.']);
end
if getDebugLevel('WriteImages')
    Vresexp = uint8((VresTV2+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV12.png'],'png');
    Vresexp = uint8((VresTV+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV.png'],'png');
end
% Evaluation
VDistTV12 = problem.M.dist(permute(V,[3,1,2]),permute(VresTV2,[3,1,2]));
METV12 = sum(sum(VDistTV12.^2))/numel(V);
VDistTV = problem.M.dist(permute(V,[3,1,2]),permute(VresTV,[3,1,2]));
METV = sum(sum(VDistTV.^2))/numel(V);
if getDebugLevel('Figures')
    figure(6); imagesc(VDistTV12); colormap(flipud(gray(1024))); axis image; axis off;    title(['Reconstruction Error TV1&2. MSE = ',num2str(METV12),'.']);
    figure(7); imagesc(VDistTV); colormap(flipud(gray(1024))); axis image; axis off;    title(['Reconstruction Error TV. MSE = ',num2str(METV12),'.']);
end
if getDebugLevel('WriteImages')
    VDistexp = uint8((VDistTV12-min(min(VDistTV12)))/(max(max(VDistTV12))-min(min(VDistTV12)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV12-error.png'],'png');
    VDistexp = uint8((VDistTV-min(min(VDistTV)))/(max(max(VDistTV))-min(min(VDistTV)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV-error.png'],'png');
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)