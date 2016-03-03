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
setDebugLevel('LevelMin',0);        % Minimal Value of all global Values
setDebugLevel('LevelMax',1000);     % Maximal Value ...
setDebugLevel('text',3);            % quite verbose text, but...
setDebugLevel('IterationStep',100); % only every Xth iteration
setDebugLevel('WriteImages',1);     % 0: no file writing, 1: file writing
setDebugLevel('LoadData',1);        % 0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0);       % 0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('time',3);            % verbose time
setDebugLevel('Figures',1);         % 0: no figure display, 1: figures are displayed
setDebugLevel('logfile',1);         % 0: no logfile 1: logfile
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
if getDebugLevel('logfile')
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
if getDebugLevel('LoadData') && exist([results,dataName,'-data.mat'],'file')
    load([results,dataName,'-data.mat']); %loads f and fn, sigma and pts
    metaData = dir([results,dataName,'-data.mat']);
    disp(['Using File Data generated ',datestr(metaData.date),'.']);
else
    %% Create Data
    f = ArtificialSPDImage(pts,jumpSize);
    fn = M.addNoise(f,sigma);
    %%
    debug('text',2,'Text',['Using Data generated ',datestr(datetime),'.']);
    if getDebugLevel('WriteData') %Write this version to file
        save([results,dataName,'-data.mat'],'f','fn','sigma','pts');
    end
end
%
% Plotinitial data
if getDebugLevel('Figures')
    figure(1);
    plotSPD(f,'GridDistance',gS,'EllipsoidPoints',20);
    title('Original Data');
    figure(2);
    plotSPD(fn,'GridDistance',gS,'EllipsoidPoints',20)
    title('Noisy Data as Input');
end
if getDebugLevel('WriteImages')
    exportSPD2Asymptote(fn,'GridDistance',gS,'ExportHeader',true,'File',[results,name,'-noisy.asy']);
    exportSPD2Asymptote(f,'GridDistance',gS,'ExportHeader',true,'File',[results,name,'-original.asy']);
end
%%Denoising with several parameters
M = SymPosDef(3);

problem.M = M;
problem.f = fn;
problem.MaxIterations = 400;
problem.Epsilon = 0;
problem.lambda = .1;
problem.M.tau = 1;
problem.M.steps = 2; %just do 5 gradient descent steps
if getDebugLevel('logfile');
    problem
    M
end
alpha = [0 5/16 3/8];
beta = [0 1/8 5/8];
debug('text',2,'Text',['Parameter range alpha (',num2str(length(alpha)),' values): ',regexprep(num2str(alpha,5), '\s*', ','),'.']);
debug('text',2,'Text',['Parameter range beta  (',num2str(length(beta)),' values): ',regexprep(num2str(beta,5), '\s*', ','),'.']);

mResults = zeros(3,3,pts,pts,length(alpha)*length(beta));
mDists = zeros(1,length(alpha)*length(beta));

%% Run Tests
% extract debug parameters for parallel kernels
for i=1:length(alpha)*length(beta)
    [j1,j2] = ind2sub([length(alpha),length(beta)],i);    
    problem.alpha = alpha(j1)*[1,1,0,0];
    problem.beta = beta(j2);
    fr = cppa_ad_2D(problem);
    mResults(:,:,:,:,i) = fr;
    mDists(i) = 1/(pts*pts)*sum(sum(problem.M.dist(f,fr)));
    if getDebugLevel('WriteImages') 
        fileStr = [results,name,'-p-',num2str(alpha(j1)),'-',num2str(beta(j2))];
        fileStr(fileStr=='.') = [];
        exportSPD2Asymptote(fr,'GridDistance',gS,'ExportHeader',true,'File',[fileStr,'.asy']);
    end
    debug('text',2,'Text',['Parameters: \alpha=',num2str(alpha(j1)),' \beta=',num2str(beta(j2)),' yield ',num2str(mDists(i)),'.']);
end

if getDebugLevel('WriteImages') %Write this version to file
    save([results,dataName,'-results.mat'],'mResults','mDists','f','fn','sigma','pts');
end
    
%% Evaluate results
[b,a] = (meshgrid(beta,alpha)); a = a(:); b = b(:);
% general minimum
[minValue,minIndex] = min(mDists);
[minJ1,minJ2] = ind2sub([length(alpha),length(beta)],minIndex);
    debug('text',1,'Text',['Minimum: Parameters: \alpha=',num2str(alpha(minJ1)),' \beta=',num2str(beta(minJ2)),' yields minimal value ',num2str(mDists(minIndex)),'.']);
minF = squeeze(mResults(:,:,:,:,minIndex));
% TV1 minimum
[minValueTV,minIndexTV] = min(mDists(b==0));
[minJ1TV,~] = ind2sub([length(alpha),length(beta)],minIndexTV);
    debug('text',1,'Text',['Minimum: Parameter: \alpha=',num2str(alpha(minJ1TV)),' yields minimal value ',num2str(minValueTV),'.']);
minFTV = squeeze(mResults(:,:,:,:,minIndex));
% just TV2 minimum
[minValueTV2,minIndexTV2] = min(mDists(a==0));
    debug('text',1,'Text',['Minimum TV2: Parameter: \beta=',num2str(beta(minIndexTV2)),' yields minimal value ',num2str(minValueTV2),'.']);
minFTV2 = squeeze(mResults(:,:,:,minIndexTV2));
%% Plot Results
if getDebugLevel('Figures')
    figure(3);
    plotSPD(minF,'GridDistance',gS,'EllipsoidPoints',12);
    title(['Best TV1&2 Result having E=',num2str(minValue)]);
    figure(4);
    plotSPD(minFTV,'GridDistance',gS,'EllipsoidPoints',12)
    title(['Best TV Result having E=',num2str(minValueTV)]);
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,'; ended ',datestr(datetime),' ---']);
    diary off;
end