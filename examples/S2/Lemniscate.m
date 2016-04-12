%
% S2-Lemniscate
% =============
%
% The Curve as a onedimensional xample on the Sphere S2. Based on an
%     Lemniscate of Bernoulli.
% in the plane projected onto the sphere, we take a look at its
% reconstruction from noisy measurements.
%
% Cf. also Section 5.1 of, where the parameters are different due to a
% different parametrization of the proximal maps
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann,
%     A second order non-smooth variational model for restoring
%     manifold-valued images,
%     SIAM Journal on Scientific Computing, 38, (1), A567?A597, 2016.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-04-22 | 2016-03-01
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
setDebugLevel('LevelMin',0);     % Lower bound
setDebugLevel('LevelMax',3);     % Upper bound
setDebugLevel('text',2);    % Not so much text debug
setDebugLevel('time',3);    % Times
setDebugLevel('WriteImages',1); %0: no file writing, 1: file writing
setDebugLevel('LoadData',1); %0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0); %0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('figures',1); %0: no figure display, 1: figures are displayed (disable e.g. for cluster/console work)
setDebugLevel('logfile',1); %0: no logfile 1: logfile
folder = ['examples',filesep,'S2',filesep]; %this files folder
results = ['Lemniscate',filesep]; % subfolder for reulsts
name = 'Lemniscate'; %fildename
dataName = 'Lemniscate'; % Data filename
%
colors = [0,0,0;.3,.3,1]'; % colors for the signals
% The following Values might get ignored, when loading data!
sigma = 6/360*2*pi; %6 degrees of Gaussian standard deviation
pts = 512; %image size in x and y
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
M = Sn(2);
%
%
% Loading or creating and saving data
if getDebugLevel('LoadData') && exist([results,dataName,'-data.mat'],'file')
    load([results,dataName,'-data.mat']); %loads f and fn
    metaData = dir([results,dataName,'-data.mat']);
    disp(['Using File Data generated ',datestr(metaData.date),'.']);
    imgSize = [size(f,2),size(f,3)];
else
    %Create Data
    t = linspace(0,2*pi,pts);
    a = pi/(2*sqrt(2)); %Lemniscate of Bernoulli
    xc = a*sqrt(2)*(cos(t)./(sin(t).^2+1));
    yc = a*sqrt(2)*(cos(t).*sin(t)./(sin(t).^2+1));
    zc = ones(size(xc));
    %f = [xc;yc;zc]./repmat(sqrt(sum([xc;yc;zc].^2,1)),[3,1]); %old
    f = M.exp(repmat([0;0;1],[1,pts]),[xc;yc;zeros(size(xc))]); %length-keeping 
    fn = M.addNoise(f,sigma);
    debug('text',2,'Text',['Using Data generated ',datestr(datetime),'.']);
    if getDebugLevel('WriteData') %Write this version to file
        save([results,dataName,'-data.mat'],'f','fn','sigma','pts');
        debug('text',2,'Text',['New Data saved to ',results,dataName,'-data.mat'.']);
    end
end
%
%
% Plot original and noisy data
X = cell(1,3);
[X{:}] = sphere(40);
lightGrey = 0.8*[1 1  1]; % It looks better if the lines are lighter
if getDebugLevel('figures') % init meshgrids
    % Show in quiver (no color coding)
    surface(X{:},'FaceColor', 'none','EdgeColor',lightGrey,'LineWidth',.1)
    hold on
    plot3(f(1,:),f(2,:),f(3,:),'o','Markersize',3) 
    plot3(fn(1,:),fn(2,:),fn(3,:),'o','Markersize',3) 
    hold off
    title('Original and Noisy Data');
    axis off;
    pause(0.02);
end
if getDebugLevel('WriteImages')>0
     fileStr = [results,name,'-noisy'];
     fileStr(fileStr=='.') = [];
     exportSphereSignals2Asymptote(cat(3,f,fn),colors,'File',[fileStr,'.asy'],'ExportHeader',true);
     fileStr = [results,name,'-original'];
     fileStr(fileStr=='.') = [];
     exportSphereSignals2Asymptote(cat(3,f,f),colors,'File',[fileStr,'.asy'],'ExportHeader',true);
end

% Generate general problem statement and settings mentioned in the paper
problem.M = M;
problem.f = fn;
problem.MaxIterations = 1000;
problem.Epsilon = 0;
problem.lambda = pi/2;
problem.M.steps = 10; %just do 10 gradient descent steps
problem.M.tau = 1;
if getDebugLevel('logfile') % for the logfile print manifold and problem.
    M        %#ok<NOPTS>
    problem  %#ok<NOPTS>
end
alpha = [0 0.55 0.7];
beta =  [0 21 29];
debug('text',2,'Text',['Parameter range alpha (',num2str(length(alpha)),' values): ',regexprep(num2str(alpha,5), '\s*', ','),'.']);
debug('text',2,'Text',['Parameter range beta  (',num2str(length(beta)),' values): ',regexprep(num2str(beta,5), '\s*', ','),'.']);

%% Iterate parameters
mDist = zeros(length(alpha)*length(beta),1);
mResults = zeros([length(alpha)*length(beta), size(f)]);
for i=1:length(alpha)*length(beta)
    % the global structure is not passed to the workers, hence we have
    % to distribute the needed DebugValues ourselves 
    [j1,j2] = ind2sub([length(alpha),length(beta)],i);    
    problem.alpha = alpha(j1);
    problem.beta = beta(j2);
    Lfr = cppa_ad_1D(problem);
    mResults(i,:,:,:) = Lfr;
    mDist(i) = 1/pts*sum(sum( problem.M.dist(f,Lfr))); %SD
    if getDebugLevel('WriteImages')>0 
        fileStr = [results,name,'-p-',num2str(alpha(j1)),'-',num2str(beta(j2))];
        fileStr(fileStr=='.') = [];
        exportSphereSignals2Asymptote(cat(3,f,Lfr),colors,'File',[fileStr,'.asy'],'ExportHeader',true);
    end
    debug('text',1,'Text',...
        ['Parameters: \alpha=',num2str(alpha(j1)),' \beta=',num2str(beta(j2)),' yield ',num2str(mDist(i)),'.']);
end

%% Evaluate results
[b,a] = (meshgrid(alpha,beta)); a = a(:); b = b(:);
% general minimum
[minValue,minIndex] = min(mDist);
[minJ1,minJ2] = ind2sub([length(alpha),length(beta)],minIndex);
    debug('text',1,'Text',['Minimum: Parameters: \alpha=',num2str(alpha(minJ1)),' \beta=',num2str(beta(minJ2)),' yields minimal value ',num2str(mDist(minIndex)),'.']);
minF = squeeze(mResults(minIndex,:,:));
% TV1 minimum
[minValueTV,minIndexTV] = min(mDist(b==0));
[minJ1TV,minJ2TV] = ind2sub([length(alpha),length(beta)],minIndexTV);
    debug('text',1,'Text',['Minimum: Parameter: \alpha=',num2str(alpha(minJ1TV)),' yields minimal value ',num2str(minValueTV),'.']);
minFTV = squeeze(mResults(minIndex,:,:));

%% Plot result
if getDebugLevel('figures')
    % Show in quiver (no color coding)
    figure;
    surface(X{:},'FaceColor', 'none','EdgeColor',lightGrey,'LineWidth',.1)
    hold on
    plot3(f(1,:),f(2,:),f(3,:),'o','Markersize',3) 
    plot3(minF(1,:),minF(2,:),minF(3,:),'o','Markersize',3) 
    hold off
    title(['Reconstruction \alpha=',num2str(alpha(minJ1)),' \beta=',num2str(beta(minJ2)),' yields ',num2str(mDist(minIndex)),'.']);
    axis off;

    figure;
    surface(X{:},'FaceColor', 'none','EdgeColor',lightGrey,'LineWidth',.1)
    hold on
    plot3(f(1,:),f(2,:),f(3,:),'o','Markersize',3) 
    plot3(minFTV(1,:),minFTV(2,:),minFTV(3,:),'o','Markersize',3) 
    hold off
    title([' Reconstruction TV: \alpha=',num2str(alpha(minJ1)),' yields ',num2str(mDist(minIndexTV)),'.']);
    axis off;
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment >',name,'< ended ',datestr(datetime),' ---']);
    diary off;
end