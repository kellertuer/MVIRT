%
% Inpainting an artificial example first and second oder differences.
% =======================================================================
%
% Demonstrating the Cyclic Proximal Point Algorithm for absolute Differences
% of first and second order on a phase valued image containing a constant
% block and a linear increase inside one block. This extends an idea from
%
% K. Papafitsoros, C. B. Schönlieb: A Combined First and Second Order
%   Variational Approach for Image Reconstruciton,
%       J. Math. Imaging. Vis. 48, pp. 308-338, 2014.
%
% to linear blocks and the phase valued setting.
%
% The comparison is published as the second numerical example (Sec. 5) of
%
%    R. Bergmann, A. Weinmann: Inpainting of Cyclic Data Using First and
%      Second Order Differences Second Order Differences, in: EMMCVPR 2015.
%
% This file can be started without any changes; it initializes the Toolbox
% itself
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-10-05 | 2016-02-24
% see LICENSE.txt

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
setDebugLevel('IterationStep',1000); %only every 1000th iteration
setDebugLevel('WriteImages',1); %0: no file writing, 1: file writing
setDebugLevel('time',3); %many time measurements
setDebugLevel('Figures',1); %0: no figure display, 1: figures are displayed (disable e.g. for cluster/console work)
setDebugLevel('logfile',1); %0: no logfile 1: logfile

dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['BlocksInpaint',filesep];
name = 'PhaseBlocks';

%% Logfile?
if getDebugLevel('logfile')
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end

%% ---------
% Data generation

pts = 128;
[X,Y] = meshgrid(-.5:1/(pts-1):.5,-.5:1/(pts-1):.5);
%
% Blocks
V = 2*pi*(2*X/0.8 + 1/2).*((abs(X)<=0.4).*(abs(Y-0.25)<=0.125));
V = V + pi/2*((abs(X)<=0.4).*(abs(Y+0.25)<=0.125));
V = symMod(V,2*pi);
%
if getDebugLevel('Figures')
    figure(1); imagesc(V,[-pi,pi]); colormap(hsv(1024));
    title(['Original Data ''',name,'''.']);
    axis image
    axis off
end
if getDebugLevel('WriteImages')
    Vexp = uint8((V+pi)/(2*pi)*255);
    imwrite(Vexp,colormap(hsv(255)),[resultsFolder,name,'.png'],'png');
end
%% Mask
% On 128^2: 0.03125 -> 26 pixel
pxt = 0.03125;
ext=4;
Mint = ( ((abs(Y+0.25)<=0.125)|(abs(Y-0.25)<=0.125))&(abs(X)>=pxt) )...
    & ~ ( (abs(Y-0.25)<=pxt) & (abs(X) <= 2*ext*pxt) );
Mexthalf = (((abs(Y+0.25)>0.125)&(abs(Y+0.25)<=0.2))...
    | ((abs(Y-0.25)>0.125)&(abs(Y-0.25)<=0.2)))...
    & (abs(X)>=ext*pxt);
Mext = (abs(Y+0.25)>0.2) & (abs(Y-0.25)>0.2);
Mask = Mexthalf|Mext|Mint;
% here, 1 means keeping data, 0 means destroyd
img = V.*Mask; %Just loose data
if getDebugLevel('Figures')
    figure(2); imagesc(Mask,[0,1]); colormap(gray); axis square
    title(['Mask: White pixels are kept, black ones are destroyed.']);
    figure(3); imagesc(img,[-pi,pi]);  colormap(hsv(1024)); axis image; axis off;
    title('input data with destroyed points');
end
if getDebugLevel('WriteImages')
    % white mask of cut out (rest transparent
    imwrite(uint8(255*(1-Mask)),[resultsFolder,name,'-mask-w.png'],'Alpha',255*(1-Mask));
    % black mask of cut out (rest transparent
    imwrite(uint8(255*Mask),[resultsFolder,name,'-mask-b.png'],'Alpha',255*(1-Mask));
end
%% Parameters for the cyclic proximal point algorithm

problem.alpha = 2*[1,1,1,1];
problem.beta=[1,1,1];
problem.lambda = pi;
problem.f = permute(img,[3,1,2]);
problem.MaxIterations = 4000;
problem.Epsilon = 10^(-9);
problem.M = S1mRn(0,1);
problem.UnknownMask = ~Mask;
problem.RegMask = ~Mask;
disp(' --- CPPA TV1&2 on R ---');
VresTV12R = permute(cppa_ad_2D(problem),[2,3,1]);
if getDebugLevel('Figures')
    figure(4); imagesc(VresTV12R,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=',num2str(problem.beta(1),4),' on R.']);
end
problem.M = S1();
disp(' --- CPPA TV1&2 ---');
VresTV12 = permute(cppa_ad_2D(problem),[2,3,1]);
if getDebugLevel('Figures')
    figure(5); imagesc(VresTV12,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=',num2str(problem.beta(1),4),'.']);
end
problem.beta = [0,0,0];
disp(' --- CPPA TV ---');
VresTV = permute(cppa_ad_2D(problem),[2,3,1]);
if getDebugLevel('Figures')
    figure(6); imagesc(VresTV,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=0.']);
end
%use a small value of alpha for initialization
problem.alpha = 0*[1,1,1,1]; problem.beta = [1,1,1];
disp(' --- CPPA TV2 ---');
VresTV2 = permute(cppa_ad_2D(problem),[2,3,1]);
if getDebugLevel('Figures')
    figure(7); imagesc(VresTV2,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=0 and \beta=',num2str(problem.beta(1)),'.']);
end
if getDebugLevel('WriteImages')
    Vresexp = uint8((VresTV12R+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV12R.png'],'png');
    Vresexp = uint8((VresTV12+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV12.png'],'png');
    Vresexp = uint8((VresTV+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV.png'],'png');
    Vresexp = uint8((VresTV2+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV2.png'],'png');
end
% Evaluation
VDistTV12 = problem.M.dist(permute(V,[3,1,2]),permute(VresTV12,[3,1,2]));
METV12 = sum(sum(VDistTV12.^2))/numel(V);
VDistTV12R = problem.M.dist(permute(V,[3,1,2]),permute(VresTV12R,[3,1,2]));
METV12R = sum(sum(VDistTV12R.^2))/numel(V);
VDistTV = problem.M.dist(permute(V,[3,1,2]),permute(VresTV,[3,1,2]));
METV = sum(sum(VDistTV.^2))/numel(V);
VDistTV2 = problem.M.dist(permute(V,[3,1,2]),permute(VresTV2,[3,1,2]));
METV2 = sum(sum(VDistTV2.^2))/numel(V);
if getDebugLevel('Figures')
    figure(8); imagesc(VDistTV12); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV1&2. MSE = ',num2str(METV12),'.']);
    figure(9); imagesc(VDistTV12R); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV1&2 on R. MSE = ',num2str(METV12R),'.']);
    figure(10); imagesc(VDistTV); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV. MSE = ',num2str(METV),'.']);
    figure(11); imagesc(VDistTV2); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV2. MSE = ',num2str(METV2),'.']);
end
if getDebugLevel('WriteImages')
    VDistexp = uint8((VDistTV12-min(min(VDistTV12)))/(max(max(VDistTV12))-min(min(VDistTV12)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV12-error.png'],'png');
    VDistexp = uint8((VDistTV12R-min(min(VDistTV12R)))/(max(max(VDistTV12R))-min(min(VDistTV12R)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV12R-error.png'],'png');
    VDistexp = uint8((VDistTV-min(min(VDistTV12)))/(max(max(VDistTV))-min(min(VDistTV12)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV-error.png'],'png');
    VDistexp = uint8((VDistTV2-min(min(VDistTV2)))/(max(max(VDistTV2))-min(min(VDistTV2)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV2-error.png'],'png');
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)
