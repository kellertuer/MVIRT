%
% Denoising a phase-valued signal with first and second oder differences.
% =======================================================================
%
% Demonstrating the Cyclic Proximal Point Algorithm for absolute Differences
% of first and second order on a 1D signal
%
% The comparison is published as the first numerical example (Sec. 5) of
%
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
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
writeData = true;
loadData = true;
writeImages=true;
showFigures = true;
useLogfile = true;

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
betas= [0,1.5,.75];
% Details for the data - these are ignored if data is loaded
sigma = 0.2;
name = ['S1-syntheticSignalN',num2str(N)];
dataName = ['S1-syntheticSignalN',num2str(N)];
%
%
% Logfile?
if useLogfile
    clc
    if exist([resultsFolder,name,'.log'],'file')
        delete([resultsFolder,name,'.log'])
    end
    diary([resultsFolder,name,'.log']);
    disp([' --- Logfile of Experiment ',name,' started ',datestr(datetime),' ---']);
end
if writeData
    [fn, fo, t] = cyclicSignal(N,0.2);
    save([resultsFolder,dataName,'.mat'],'fn','fo','t','sigma');
elseif loadData
    load([resultsFolder,dataName,'.mat'],'fn','fo','t','sigma');
    metaData = dir([resultsFolder,dataName,'.mat']);
    doisp(['Using File Data generated ',datestr(metaData.date),'.']);
    N = length(fn);
else
    error('Either Loading or Creating(Wirting) Data must be set');
end
%%
if showFigures
    figure(1); plot(t,fo,'.b',t,fn,'.k');
    ylim([-pi,pi]);
    title(['Original & Noise (\sigma=',num2str(sigma),').']);
end
M=S1;
fresults = zeros(length(alphas),N);

problem.M = S1();
problem.lambdaIterate = @(iter) pi/iter;
problem.stoppingCriterion = stopCritMaxIterEpsilonCreator(problem.M,iter,epsilon);
problem.f = fn;

for i=1:length(alphas)
    disp(['- TV1&2 Minimization by CPPA, alpha=',num2str(alphas(i)),', beta=',num2str(betas(i)),' -']);    
    problem.alpha = alphas(i);
    problem.beta = betas(i);
    tic
    fresults(i,:) = CPP_AdditiveTV12(problem);
    toc
end
    disp(['- TV1 Minimization by CPPA on R, alpha=',num2str(alphas(1)),', beta=0 -']);
    M2 = Rn(1);    
    problem.M = M2;
    problem.alpha=alphas(1);
    problem.beta = 0;
    tic
    fResult = CPP_AdditiveTV12(problem);
    toc;
    
    
if showFigures
    for i=1:length(alphas)
        figure(i+1);
        plot(t,fo,'.b',t,fresults(i,:),'.k');
        ylim([-pi,pi]);
        title(['Result of parameters \alpha=',num2str(alphas(i)),', \beta=',num2str(betas(i)),...
            ' MSE ',num2str(sum(M.dist(fo,fresults(i,:)).^2)*1/length(fo(:)),'%6.5f')]);
    end
    figure(length(alphas)+2)
    plot(t,fo,'.b',t,fResult,'.k');
    ylim([-pi,pi]);
    title(['Result of real-valued TV parameters \alpha=',num2str(alphas(i)),', \beta=0',...
        ' MSE (on R) ',num2str(sum(M2.dist(fo,fResult).^2)*1/length(fo(:)),'%6.5f')]);
end
if writeImages
    fo2 = fo; fo2(157:344) = fo2(157:344)+2*pi; %adapt to general N?
    T = table(...
      t.',fo.',fo2.',fn.',fresults(1,:).',fresults(2,:).',fresults(3,:).',fResult.',...
      'VariableNames',{'x','orig','origS','noisy','TV','TV2','TV1a2','TV1R'});
    writetable(T,[resultsFolder,name,'.dat']);
end
%% End logfile
if useLogfile
    disp([' --- Logfile of Experiment ',name,' ended ',datestr(datetime),' ---']);
    diary off;
end
cd(start)