%
% SPDImage
% ========
%
% Illustrating the denoising capabilities of the first and second
% differences order on the manifold of symmetric postive definite matrices
%
% This Matlab script is used to generate Figure 5.4 in
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann: A second order non-smooth
%   variational model for restoring manifold-valued images, 2015.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-04-02 | 2015-03-02
% see LICENSE.txt

clear problem; close all; clc
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')
%
%
%% Settings & Variables
writeImages = true;
loadData = true;
showFigures = true;
useLogfile = true;
format compact
%
folder = ['examples',filesep,'SPD',filesep];
results = ['SPDImage',filesep];
name = 'SPDImage';
dataName = 'SPDImage';
% The following Values might get ignored, when loading data!
sigma = 0.03; %6 degrees of Gaussian standard deviation
pts = 25; %image size in x and y
jumpSize = 1.5;
gS = 4; %GridDistance for all plots to have the same presentation
%
%
%% Initialization
if useLogfile
    if exist([results,name,'.log'],'file')
        delete([results,name,'.log'])
    end
    diary([results,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
M = SymPosDef(3);
%
%
% Loading or creating and saving data
if loadData && exist([results,dataName,'-data.mat'],'file')
    load([results,dataName,'-data.mat']); %loads f and fn, sigma and pts
    metaData = dir([results,dataName,'-data.mat']);
    disp(['Using File Data generated ',datestr(metaData.date),'.']);
else
    %% Create Data
    f = ArtificialSPDImage(pts,jumpSize);
    fn = M.addNoise(f,sigma);
    %%
    disp(['Using Data generated ',datestr(datetime),'.']);
    if writeData %Write this version to file
        save([results,dataName,'-data.mat'],'f','fn','sigma','pts');
    end
end
%
% Plotinitial data
if showFigures
    figure(1);
    plotSPD(f,'GridDistance',gS,'EllipsoidPoints',20);
    title('Original Data');
    figure(2);
    plotSPD(fn,'GridDistance',gS,'EllipsoidPoints',20)
    title('Noisy Data as Input');
end
if writeImages
    exportSPD2Asymptote(fn,'GridDistance',gS,'ExportHeader',true,'File',[results,name,'-noisy.asy']);
    exportSPD2Asymptote(f,'GridDistance',gS,'ExportHeader',true,'File',[results,name,'-original.asy']);
end
%%Denoising with several parameters
M = SymPosDef(3);

problem.M = M;
problem.f = fn;
problem.lambda = .1;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,800,0);
problem.Debug = 100;
if useLogfile
    problem %#ok<NOPTS>
    M       %#ok<NOPTS>
end
alpha = 0:1/10:1;
beta = 0:1/10:1;
disp(['Parameter range alpha (',num2str(length(alpha)),' values): ',regexprep(num2str(alpha,5), '\s*', ','),'.']);
disp(['Parameter range beta  (',num2str(length(beta)),' values): ',regexprep(num2str(beta,5), '\s*', ','),'.']);

mResults = zeros(3,3,pts,pts,length(alpha)*length(beta));
mDists = zeros(1,length(alpha)*length(beta));

%% Run Tests
% extract debug parameters for parallel kernels
minTV = Inf;
minTVa = Inf;
minTV12 = Inf;
minTV12a = Inf;
minTV12b = Inf;
for i=1:length(alpha)*length(beta)
    [j1,j2] = ind2sub([length(alpha),length(beta)],i);    
    problem.alpha = alpha(j1)*[1,1,0,0];
    problem.beta = beta(j2);
    fr = CPP_AdditiveTV12(problem);
    mResults(:,:,:,:,i) = fr;
    mDists(i) = 1/(pts*pts)*sum(sum(problem.M.dist(f,fr)));
    t = mDists(i);
    if problem.beta==0 && t < minTV
        minTV = t;
        minTVa = problem.alpha;
        disp(['Min TV ',num2str(t),'.']);
        TVsig = fr;
    end
    if t < minTV12
        minTV12 = t;
        minTV12a = problem.alpha;
        minTV121b = problem.beta;
        disp(['Min TV12 ',num2str(t),'.']);
        TV12sig = fr;
    end
    if writeImages
        fileStr = [results,name,'-p-',num2str(alpha(j1)),'-',num2str(beta(j2))];
        fileStr(fileStr=='.') = [];
        exportSPD2Asymptote(fr,'GridDistance',gS,'ExportHeader',true,'File',[fileStr,'.asy']);
    end
    disp(['Parameters: \alpha=',num2str(alpha(j1)),' \beta=',num2str(beta(j2)),' yield ',num2str(mDists(i)),'.']);
end

if writeImages %Write this version to file
    save([results,dataName,'-results.mat'],'mResults','mDists','f','fn','sigma','pts');
end
    
%% Evaluate results
disp(['Minimum: Parameters: \alpha=',num2str(minTV12a),' \beta=',num2str(minTV12b),' yields minimal value ',num2str(minTV12),'.']);
disp(['Minimum: Parameter: \alpha=',num2str(minTVa),' yields minimal value ',num2str(minTV),'.']);
%% Plot Results
if showFigures
    figure(3);
    plotSPD(TVsig,'GridDistance',gS,'EllipsoidPoints',12);
    title(['Best TV1&2 Result having E=',num2str(minValue)]);
    figure(4);
    plotSPD(TV12sig,'GridDistance',gS,'EllipsoidPoints',12)
    title(['Best TV Result having E=',num2str(minValueTV)]);
end
%% End logfile
if useLogfile
    disp([' --- Logfile of Experiment ',name,'; ended ',datestr(datetime),' ---']);
    diary off;
end