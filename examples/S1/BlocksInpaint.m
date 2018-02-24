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
% R. Bergmann ~ 2015-10-05 | 2018-02-24
% see LICENSE.txt

% Changelog
% 2018-02-24 RB streamline file

start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')
%
%
%% Settings for this script
writeLogfile = true; %write a logfile into the results folder
writeImages  = true; % write result images
showFigures  = true; %activate (1) or deactivate (0) figures
% Folders
dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['BlocksInpaint',filesep];
name = 'PhaseBlocks';

%% Logfile?
if writeLogfile
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end

%% Data generation
pts = 128;
[X,Y] = meshgrid(-.5:1/(pts-1):.5,-.5:1/(pts-1):.5);
%
% Blocks
V = 2*pi*(2*X/0.8 + 1/2).*((abs(X)<=0.4).*(abs(Y-0.25)<=0.125));
V = V + pi/2*((abs(X)<=0.4).*(abs(Y+0.25)<=0.125));
V = symMod(V,2*pi);
%
if showFigures
    figure(1); imagesc(V,[-pi,pi]); colormap(hsv(1024));
    title(['Original Data ''',name,'''.']);
    axis image
    axis off
end
if writeImages
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
if showFigures
    figure(2); imagesc(Mask,[0,1]); colormap(gray); axis square
    title(['Mask: White pixels are kept, black ones are destroyed.']);
    figure(3); imagesc(img,[-pi,pi]);  colormap(hsv(1024)); axis image; axis off;
    title('input data with destroyed points');
end
if writeImages
    % white mask of cut out (rest transparent
    imwrite(uint8(255*(1-Mask)),[resultsFolder,name,'-mask-w.png'],'Alpha',255*(1-Mask));
    % black mask of cut out (rest transparent
    imwrite(uint8(255*Mask),[resultsFolder,name,'-mask-b.png'],'Alpha',255*(1-Mask));
end
%% Parameters for the cyclic proximal point algorithm

%% Init problem structure
problem.M = Rn(1);
problem.alpha = 2*ones(2);
problem.beta = [1,1;0,1];
problem.f = permute(img,[3,1,2]);
problem.UnknownMask = ~Mask;
problem.FixedMask = Mask;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M, 8000,10^(-9));
tic
VresTV12R = permute(CPP_AdditiveTV12(problem),[2,3,1]);
toc
%%
if showFigures
    figure(4); imagesc(VresTV12R,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1,1),4),' and \beta=',num2str(problem.beta(1,1),4),' on R.']);
end
problem.M = S1();
disp(' --- CPPA TV1&2 ---');
tic
VresTV12 = permute(CPP_AdditiveTV12(problem),[2,3,1]);
toc
if showFigures
    figure(5); imagesc(VresTV12,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=',num2str(problem.beta(1),4),'.']);
end
problem.beta = zeros(2);
disp(' --- CPPA TV ---');
tic
VresTV = permute(CPP_AdditiveTV12(problem),[2,3,1]);
toc
if showFigures
    figure(6); imagesc(VresTV,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=0.']);
end
%use a small value of alpha for initialization
problem.alpha = zeros(2); problem.beta = [1,1;0,1];
disp(' --- CPPA TV2 ---');
tic
VresTV2 = permute(CPP_AdditiveTV12(problem),[2,3,1]);
toc
if showFigures
    figure(7); imagesc(VresTV2,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=0 and \beta=',num2str(problem.beta(1)),'.']);
end
if writeImages
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
if showFigures
    figure(8); imagesc(VDistTV12); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV1&2. MSE = ',num2str(METV12),'.']);
    figure(9); imagesc(VDistTV12R); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV1&2 on R. MSE = ',num2str(METV12R),'.']);
    figure(10); imagesc(VDistTV); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV. MSE = ',num2str(METV),'.']);
    figure(11); imagesc(VDistTV2); colormap(flipud(gray(1024))); axis image; axis off;
    title(['Reconstruction Error TV2. MSE = ',num2str(METV2),'.']);
end
if writeImages
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
if writeLogfile
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)