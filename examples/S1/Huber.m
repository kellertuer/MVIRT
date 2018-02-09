 %
% Comparing the CPPA-prox of absolute differences with a Huber relaxation
% =======================================================================
%
% This example compares the denoising with respect to an Huber type
% regularization to the first and second order TV-type regularization for
% phase valued data.
%
% This example is published in
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-06-19 | 2016-02-24
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
setDebugLevel('IterationStep',1000); %only every 50th iteration
setDebugLevel('time',2);
% The following two variables do the following
% LoadData: If this is set to 1,
setDebugLevel('LoadData',1); %0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0); %0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('Figures',1); %0: no figure display, 1: figures are displayed (disable e.g. for cluster/console work)
setDebugLevel('logfile',1); %0: no logfile 1: logfile
format compact
%
%
% Strings (only change these if you're knowing what you do
dataFolder = ['Huber',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['Huber',filesep];
N = 500;
sigma = 0.3;
name = ['S1-HuberN',num2str(N)];
dataName = ['S1-HuberN',num2str(N)];
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
if getDebugLevel('WriteData') % Create new data
    [fn,fo,t] = CyclicPiecewiseSignal(N,sigma);
    save([resultsFolder,dataName,'.mat'],'fn','fo','t','sigma');
elseif getDebugLevel('LoadData')
    load([resultsFolder,dataName,'.mat'],'fn','fo','t','sigma');
    metaData = dir([resultsFolder,dataName,'.mat']);
    debug('text',3,'Text',['Using File Data generated ',datestr(metaData.date),'.']);
    N = length(fn);
else
    error('Either Loading or Creating(Wirting) Data must be set');
end
M = S1;
%%
format compact
figure(1); plot(t,fo,'.b',t,fn,'.k');
ylim([-pi,pi]);
title(['Original & Noise (\sigma=',num2str(sigma),').']);

% Set up parameters
clear problem;
problem.M = S1();
problem.alpha = 1/2;
problem.beta = 1;
problem.f = fn;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,4000,0);
problem.lambdaIterate = @(iter) pi/iter;

% Algorithm parameters
iter = 4000;
epsilon = 0;

alpha=1/2; beta=1;
tic
frTV12 = CPP_AdditiveTV12(problem);
toc

MSETV = sum(M.dist(fo,frTV12).^2)*1/N;
if getDebugLevel('Figures')
    figure(2); plot(t,fo,'.b',t,frTV12,'.k');
    ylim([-pi,pi]);
    title(['TV1&2 Minimization by CPPA, alpha=',num2str(alpha),', beta=',num2str(beta),' MSE=',num2str(MSETV)]);
end

problemH.M = S1();
problemH.alpha = 1/2;
problemH.tau = 4*sqrt(2);
problemH.omega = pi/16;
problemH.f = fn;
problemH.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,4000,0);
problemH.lambdaIterate = @(iter) pi/iter;
frH = CPP_HuberTV(problemH);

MSEH = sum(M.dist(fo,frH).^2)*1/N;
if getDebugLevel('Figures')
    figure(3); plot(t,fo,'.b',t,frH,'.k');
    ylim([-pi,pi]);
    title(['Huber-type TV Minimization by CPPA, alpha=',num2str(problemH.alpha),', tau=',num2str(problemH.tau),', omega=',num2str(problemH.omega),' MSE=',num2str(MSETV)]);
end