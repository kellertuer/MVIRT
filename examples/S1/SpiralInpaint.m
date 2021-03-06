%
% Inpainting an artificial example first and second oder differences.
% =======================================================================
%
% Demonstrating the Cyclic Proximal Point Algorithm for absolute Differences
% of first and second order on a phase valued image (a spiral adapted from
% Ghiglia, D.C., Pritt, M.D.: Two-dimensional phase unwrapping:
%   theory, algorithms, and software. Wiley, 1998.
% This exaple also illustrates the advantage of a second order phase valued
% TV when the last data was previsouly linear.
%
% The comparison is published as the first numerical example (Sec. 5) of
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
useLogfile  = true;
showFigures = true;
writeImages = true;

dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['SpiralInpaint',filesep];
name = 'PhaseSpiral';
load([dataFolder,name,'.mat']);

%% Logfile?
if useLogfile
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end

if showFigures
    figure(1); imagesc(spiral,[-pi,pi]); colormap(hsv(1024));
    title(['Original Data ''',name,'''.']);
    axis image
    axis off
end
if writeImages
    Vexp = uint8((spiral+pi)/(2*pi)*255);
    %Export original image
    imwrite(Vexp,colormap(hsv(255)),[resultsFolder,name,'.png'],'png');
end
%% ------
% Mask generation
% Destroy two thirds of all rows/cols -> just 1/9 kept
% here, 1 means keeping data, 0 means destroyd
M = ones(size(spiral));
M(1:3:size(spiral,1),:) = 0;
M(2:3:size(spiral,1),:) = 0;
%
M(:,1:3:size(spiral,2)) = 0;
M(:,2:3:size(spiral,2)) = 0;
M = M>0; %binarize
img = spiral.*M; %Just loose data
if showFigures
    figure(2); imagesc(M,[0,1]); colormap(gray); axis square
    pr = 100 * sum(~M(:))/length(M(:));
    title(['Mask: White pixels are kept, black ones are destroyed. (',num2str(pr,3),'% destroyed).']);
    figure(3); imagesc(img,[-pi,pi]);  colormap(hsv(1024)); axis square
    title('input data with destroyed points');
end
if writeImages
    % white mask of cut out (rest transparent
    imwrite(uint8(255*(1-M)),[resultsFolder,name,'-mask-w.png'],'Alpha',255*(1-M));
    % black mask of cut out (rest transparent
    imwrite(uint8(255*M),[resultsFolder,name,'-mask-b.png'],'Alpha',255*(1-M));
end
%% Parameters for the cyclic proximal point algorithm

problem.M = S1();
problem.alpha = eye(2);
problem.beta = [1,1;0,1];
problem.lambdaIterate = @(iter) pi/2/iter;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,4000,0);
problem.f = permute(img,[3,1,2]);
problem.UnknownMask = ~M;
problem.FixedMask = M;
tic
VresTV2 = permute(CPP_AdditiveTV12(problem),[2,3,1]);
toc

if showFigures
    figure(4); imagesc(VresTV2,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=',num2str(problem.beta(1),4),'.']);
end
problem.beta = zeros(2);
tic
VresTV = permute(CPP_AdditiveTV12(problem),[2,3,1]);
toc
if showFigures
    figure(5); imagesc(VresTV,[-pi,pi]); colormap(hsv(1024));axis image; axis off;
    title(['Result of the reconstruction using \alpha=',num2str(problem.alpha(1),4),' and \beta=0.']);
end
if writeImages
    Vresexp = uint8((VresTV2+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV12.png'],'png');
    Vresexp = uint8((VresTV+pi)/(2*pi)*255);
    imwrite(Vresexp,colormap(hsv(255)),[resultsFolder,name,'-impaintedTV.png'],'png');
end
% Evaluation
VDistTV12 = problem.M.dist(permute(spiral,[3,1,2]),permute(VresTV2,[3,1,2]));
METV12 = sum(sum(VDistTV12.^2))/length(spiral(:));
VDistTV = problem.M.dist(permute(spiral,[3,1,2]),permute(VresTV,[3,1,2]));
METV = sum(sum(VDistTV.^2))/length(spiral(:));
if showFigures
    figure(6); imagesc(VDistTV12); colormap(flipud(gray(1024))); axis image; axis off;    title(['Reconstruction Error TV1&2. MSE = ',num2str(METV12),'.']);
    figure(7); imagesc(VDistTV); colormap(flipud(gray(1024))); axis image; axis off;    title(['Reconstruction Error TV. MSE = ',num2str(METV12),'.']);
end
if writeImages
    VDistexp = uint8((VDistTV12-min(min(VDistTV12)))/(max(max(VDistTV12))-min(min(VDistTV12)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV12-error.png'],'png');
    VDistexp = uint8((VDistTV-min(min(VDistTV)))/(max(max(VDistTV))-min(min(VDistTV)))*255);
    imwrite(VDistexp,colormap(flipud(gray(255))),[resultsFolder,name,'-impainting-TV-error.png'],'png');
end
%% End logfile
if useLogfile
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)