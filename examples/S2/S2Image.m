%
% S2 Image
% ========
%
% An example of a (unit) vector field, i.e. an image having data items on
% S2, which are obstructed by noise.
%
% See also Section 5.1 and Figure 5.2 of
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann: A second order non-smooth
%   variational model for restoring manifold-valued images, 2015.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-03-09 | 2016-03-01
% see LICENSE.txt
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')
%
%
%% Settings
writeImages = true;
loadData = true;
writeData = false;
showFigures = true;
useLogfile = true;
format compact
folder = ['examples',filesep,'S2',filesep];
results = ['S2Image',filesep];
name = 'S2Image';
dataName = 'S2Image';
% The following Values are ignored, when loading data.
sigma = 16/360*2*pi; %10 degrees of Gaussian
pts = 64; %image size in x and y
surroundings = 2.5; % while in tY we have [0,2*pi] in tX we have [0,2*pi*surroundings]
stepSize = pi/4; % Size of each step
steps = [3,3];
%
%
% Logfile
if useLogfile
    clc
    if exist([results,name,'.log'],'file')
        delete([results,name,'.log'])
    end
    diary([results,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
M = Sn(2);
%
%
% Loading or creating and saving data
if loadData && exist([results,dataName,'-data.mat'],'file')
    load([results,dataName,'-data.mat']); %loads f and fn
    metaData = dir([results,dataName,'-data.mat']);
    disp(['Using loaded Data generated ',datestr(metaData.date),'.']);
    imgSize = [size(f,2),size(f,3)];
else
    f = ArtificialS2Data(pts,surroundings,steps,stepSize);
    fn = M.addNoise(f,sigma);
    disp(['Using _new_ Data generated ',datestr(datetime),'.']);
    if writeData %Write this version to file
        save([results,dataName,'-data.mat'],'f','fn','sigma','pts');
    end
end
%
%
% Plot original and noisy data
if showFigures % init meshgrids
    % Show in quiver (no color coding)
    [Xpl,Ypl] = meshgrid(1:pts);
    Zpl = zeros(size(Xpl));
end
if showFigures
    % Show in quiver (no color coding)
    figure(1);
    quiver3(Xpl,Ypl,Zpl,...
        permute(f(1,:,:),[2,3,1]),...
        permute(f(2,:,:),[2,3,1]),...
        permute(f(3,:,:),[2,3,1])...'AutoScale','off'
    ); title('Original Data'); axis image; daspect([1,1,20/pts]); axis off
    figure(2); 
    quiver3(Xpl,Ypl,Zpl,...
        permute(fn(1,:,:),[2,3,1]),...
        permute(fn(2,:,:),[2,3,1]),...
        permute(fn(3,:,:),[2,3,1])... 'AutoScale','off'
    ); title('Noisy data'); axis image; daspect([1,1,20/pts]); axis off
end
%Export original and noisy
if writeImages
    exportSphere2Asymptote(f,'File',[results,name,'-original.asy'],'ExportHeader',true);
    exportSphere2Asymptote(fn,'File',[results,name,'-noisy.asy'],'ExportHeader',true);
end
%
%
%% General Parameters
problem.M = M;
problem.f = fn;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,800,0);
problem.Debug = 100;
problem.lambda = pi/2;

alpha = [0 0.22];
beta = [0 29.5];
disp(['Parameter range alpha (',num2str(length(alpha)),' values): ',regexprep(num2str(alpha,5), '\s*', ','),'.']);
disp(['Parameter range beta  (',num2str(length(beta)),' values): ',regexprep(num2str(beta,5), '\s*', ','),'.']);

mDist = zeros(length(alpha)*length(beta),1);
mResults = zeros([size(f), length(alpha)*length(beta)]);
minTV = Inf;
minTVp = [];
minTV12 = Inf;
minTV12p = [];
for i=1:length(alpha)*length(beta)
    [j1,j2] = ind2sub([length(alpha),length(beta)],i);    
    problem.alpha = alpha(j1);
    problem.beta = beta(j2);
    Lfr = CPP_AdditiveTV12(problem);
    mResults(:,:,:,i) = Lfr;
    mDist(i) = 1/pts^2* sum(sum(problem.M.dist(f,Lfr)));
    if minDist(i) < minTV12
        minTV12p = [problem.alpha,problem.beta];
        minTV2 = mDist(i);
        minTV12F = Lfr;
    end
    if mDist(i) < minTV && problem.beta==0
        minTVp = problem.alpha;
        minTV = minDist(i);
        minTVF = Lfr;
    end
    if writeImages
        fileStr = [results,name,'-p-',num2str(alpha(j1)),'-',num2str(beta(j2))];
        fileStr(fileStr=='.') = [];
        exportSphere2Asymptote(Lfr,'File',[fileStr,'.asy'],'ExportHeader',true);
    end
    disp(['Parameters: \alpha=',num2str(alpha(j1)),' \beta=',num2str(beta(j2)),' yield ',num2str(mDist(i)),'.']);
end
disp(['Minimum: Parameters: \alpha=',num2str(minTV12p(1)),' \beta=',num2str(minTV12p(2)),' yields minimal value ',num2str(minTV12),'.']);
disp(['Minimum TV: Parameter: \alpha=',num2str(minTVp),' yields minimal value ',num2str(minTV),'.']);

%%
if showFigures
    % Show in quiver (no color coding)
    figure;
    quiver3(Xpl,Ypl,Zpl,...
        permute(minTV12F(1,:,:),[2,3,1]),...
        permute(minTV12F(2,:,:),[2,3,1]),...
        permute(minTV12F(3,:,:),[2,3,1])...'AutoScale','off'
    ); title(['''Optimal'' Reconstruction \alpha=',num2str(minTV12p(1)),' \beta=',num2str(minTV12p(2)),' yields minimal value ',num2str(minTV12),'.']);
    axis image; daspect([1,1,20/pts]); axis off
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,'; ended ',datestr(datetime),' ---']);
    diary off;
end