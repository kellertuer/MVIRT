%
% Denoising a phase-valued signal with first and second oder differences.
% =======================================================================
%
% Demonstrating the Cyclic Proximal Point Algorithm for absolute Differences
% of first and second order on a 1D signal
%
% The comparison is published as the first numerical example (Sec. 5) of
%
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann: Second Order Differences
%        of S1-valued Data and Applications in Variational Image Denoising.
%        J. SIAM Imaging Sci. 7(4), 2916?2953,
%
% This file can be started without any changes; it initializes the Toolbox
% itself
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-03-24 | 2016-02-24
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

format compact
%
%
% Strings (only change these if you're knowing what you do
dataFolder = ['syntheticSignal',filesep];
folder = ['examples',filesep,'S1',filesep];
resultsFolder = ['syntheticSignal',filesep];
N = 500;
% Algorithm parameters
iter = 4000;
epsilon = 0;
alphas = [0.75,0,0.25];
betas= [0,1.5,0.5];
% Details for the data - these are ignored if data is loaded
sigma = 0.2;
name = ['S1-syntheticSignalN',num2str(N)];
dataName = ['S1-syntheticSignalN',num2str(N)];
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
if getDebugLevel('WriteData') % Create new data
    [fn, fo, t] = cyclicSignal(N,0.2);
    save([resultsFolder,dataName,'.mat'],'fn','fo','t','sigma');
elseif getDebugLevel('LoadData')
    load([resultsFolder,dataName,'.mat'],'fn','fo','t','sigma');
    metaData = dir([resultsFolder,dataName,'.mat']);
    debug('text',3,'Text',['Using File Data generated ',datestr(metaData.date),'.']);
    N = length(fn);
else
    error('Either Loading or Creating(Wirting) Data must be set');
end
%%
if getDebugLevel('Figures')
    figure(1); plot(t,fo,'.b',t,fn,'.k');
    ylim([-pi,pi]);
    title(['Original & Noise (\sigma=',num2str(sigma),').']);
end
M=S1;
fresults = zeros(length(alphas),N);
problem.M = S1; problem.f = fn; problem.lambda = pi;
problem.MaxIterations = iter; problem.Epsilon = epsilon;
for i=1:length(alphas)
    debug('text',2,'Text',['- TV1&2 Minimization by CPPA, alpha=',num2str(alphas(i)),', beta=',num2str(betas(i)),' -']);
    problem.alpha = alphas(i);
    problem.beta = betas(i);
    fresults(i,:) = cppa_ad_1D(problem);
end
if getDebugLevel('Figures')
    for i=1:length(alphas)
        figure(i+1);
        plot(t,fo,'.b',t,fresults(i,:),'.k');
        ylim([-pi,pi]);
        title(['Result of parameters \alpha=',num2str(alphas(i)),', \beta=',num2str(betas(i)),...
            ' MSE ',num2str(sum(M.dist(fo,fresults(i,:)).^2)*1/length(fo(:)),'%6.5f')]);
    end
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)