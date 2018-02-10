%
% Camino Slice 28
% ===============
%
% Denoising the Slice 28 of the Camino Data Set available at
% http://cmic.cs.ucl.ac.uk/camino//index.php?n=Tutorials.DTI
% (the file dt.Bdouble from step 5 used here as Caminodt.Bdouble)
%
% This file requires either a previous run, such that CaminoSlice28.mat is
% prtesent in the Camino/ sub folder or that the Caminodt.Bdouble is provided
% in the data/ folder to generate the CaminoSlice28.mat file.
%
% See also Section 5.2 and Figure 5.5 of
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann: A second order non-smooth
%   variational model for restoring manifold-valued images,
%   SIAM Journal on Scientific Computing. 38, (1), A567ÐA597, 2016.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-04-01 | 2016-03-01
% see LICENSE.txt

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
setDebugLevel('IterationStep',1);   % only every Xth iteration
setDebugLevel('WriteImages',1);     % 0: no file writing, 1: file writing
setDebugLevel('LoadData',1);        % 0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0);       % 0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('time',3);            % verbose time
setDebugLevel('Figures',1);         % 0: no figure display, 1: figures are displayed
setDebugLevel('logfile',1);         % 0: no logfile 1: logfile
format compact
%% Strings
dataFolder = ['..',filesep,'..',filesep,'data',filesep];
folder = ['examples',filesep,'SPD',filesep];
% relative to folder
resultsFolder = ['Camino',filesep];

% Data details - ignored if only loading data
sliceNum = 28;
name = ['CaminoSlice',num2str(sliceNum)];
dataName = ['CaminoSlice',num2str(sliceNum)];
rangeX = 1:112;
rangeY = 1:112;
subRangeX = 28:87;
subRangeY = 34:73;
fn = [];

%% Initialization
if getDebugLevel('logfile')
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
if getDebugLevel('WriteData') % recreate data from original data and save
    %% Load and reformate Data
    % Open Camino Data file from the tutorial at http://cmic.cs.ucl.ac.uk/camino/index.php?n=Tutorials.DTI
    if ~(exist([dataFolder,'Caminodt.Bdouble'],'file')==2) %Original file missing
        error('Camino:FileNotFound',['The File ''Caminodt.Bdouble'' is missing in the folder ''data/''',...
            'of the MVDP project.\n Please follow the tutorial at\n',...
            'http://cmic.cs.ucl.ac.uk/camino//index.php?n=Tutorials.DTI\n',...
            'and save the file fro step 5 as ''Caminodt.Bdouble'' in ''data/''']);
    end
    fid = fopen([dataFolder,'Caminodt.Bdouble'], 'r', 'b');
    dataFile = fread(fid, 'double');
    fclose(fid);
    % Reshape into 8blocks of data
    data=reshape(dataFile,8,size(dataFile,1)/8);
    % remove each first 2 to keep just 6 values and reshape to known data size
    data=data(3:8,:);
    l=112;
    m=112;
    n=50;
    data = reshape(data,[6,l,m,n]);
    % create tensor image from data
    dataDTI  = zeros(3,3,l,m,n);
    for i=1:l             % 1 2 3    | Ordering of the elements for the matrix
        for j=1:m         % 2 4 5
            for k=1:n     % 3 5 6
                dataDTI(1,1,i,j,k) = data(1,i,j,k);
                dataDTI(1,2,i,j,k) = data(2,i,j,k);
                dataDTI(1,3,i,j,k) = data(3,i,j,k);
                dataDTI(2,1,i,j,k) = data(2,i,j,k);
                dataDTI(2,2,i,j,k) = data(4,i,j,k);
                dataDTI(2,3,i,j,k) = data(5,i,j,k);
                dataDTI(3,1,i,j,k) = data(3,i,j,k);
                dataDTI(3,2,i,j,k) = data(5,i,j,k);
                dataDTI(3,3,i,j,k) = data(6,i,j,k);
            end
        end
    end
    fn = dataDTI(:,:,rangeX,rangeY,sliceNum);
    l = length(rangeX); m = length(rangeY);
    save([resultsFolder,dataName,'.mat'],'fn','l','m');
end
if getDebugLevel('LoadData')
    load([resultsFolder,dataName,'.mat'],'fn','l','m');
end
if sum(size(fn))==0
    error('Both Writing and Loading Data are disabled, no data was loaded.');
end
rM = ~permute(all(reshape(fn,9,l,m)==0,1),[2,3,1]);
if getDebugLevel('Figures')
    figure(1);
    plotSPD(fn/max(fn(:)),'GridDistance',0.25,'EllipsoidPoints',20);
    title('Original data.');
    figure(2);
    plotSPD(fn(:,:,subRangeX,subRangeY)/max(fn(:)),'GridDistance',0.25,'EllipsoidPoints',20);
    title('Part of the original data.');
    pause(0.02);
end
%Export original data
if getDebugLevel('WriteImages')
    exportSPD2Asymptote(fn/max(fn(:)),'File',[resultsFolder,name,'-original.asy'],'GridDistance',0.25,'ExportHeader',true);
    exportSPD2Asymptote(fn(:,:,subRangeX,subRangeY)/max(fn(:)),'File',[resultsFolder,name,'-original-sub.asy'],'GridDistance',0.25,'ExportHeader',true);
end
%
%
%% General Parameters for all problems
M = SymPosDef(3);
problem.M = M;
problem.f = fn;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,800,10^(-5));
problem.lambda = pi/2;
problem.alpha = 0.05;
problem.beta = 0.1;
if getDebugLevel('logfile') %display parameters in logfile
    problem     %#ok<NOPTS>
    problem.M   
end
%% CPPA
fr = CPP_AdditiveTV12(problem);
%% Export Results
if getDebugLevel('WriteImages')>0
    fileStr = [resultsFolder,name,'-p-',num2str(problem.alpha),'-',num2str(problem.beta)];
    fileStr(fileStr=='.') = [];
    exportSPD2Asymptote(fr/max(fr(:)),'File',[fileStr,'.asy'],'GridDistance',0.25,'ExportHeader',true);
    exportSPD2Asymptote(fr(:,:,subRangeX,subRangeY)/max(fr(:)),'File',[fileStr,'sub.asy'],'GridDistance',0.25,'ExportHeader',true);
    save([fileStr,'.mat'],'fr')
end
if getDebugLevel('Figures')
    figure(3);
    plotSPD(fr/max(fr(:)),'GridDistance',0.25,'EllipsoidPoints',20);
    title(['Denoised Data. Parameters: \alpha=',num2str(problem.alpha),' \beta=',num2str(problem.beta),'.']);
    figure(4)
    plotSPD(fr(:,:,subRangeX,subRangeY)/max(fr(:)),'GridDistance',0.25,'EllipsoidPoints',20);
    title(['Denoised part. Parameters: \alpha=',num2str(problem.alpha),' \beta=',num2str(problem.beta),'.']);
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,'; ended ',datestr(datetime),' ---']);
    diary off;
end
