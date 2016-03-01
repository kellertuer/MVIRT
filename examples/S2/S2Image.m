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
clear problem; close all; clc
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
setDebugLevel('IterationStep',100);   % only every Xth iteration
setDebugLevel('WriteImages',1);     % 0: no file writing, 1: file writing
setDebugLevel('LoadData',1);        % 0: generate new data 1: load existing data (if it exists), (other wise it is generated)
setDebugLevel('WriteData',0);       % 0: do not write data to file 1: write data to file (overwrites last experiment data!)
setDebugLevel('time',3);            % verbose time
setDebugLevel('Figures',1);         % 0: no figure display, 1: figures are displayed
setDebugLevel('logfile',1);         % 0: no logfile 1: logfile
format compact
folder = ['examples',filesep,'S2',filesep];
results = ['S2Image',filesep];
name = 'S2Image';
dataName = 'S2Image';
% The following Values might get ignored, when loading data!
sigma = 16/360*2*pi; %10 degrees of Gaussian
pts = 64; %image size in x and y
surroundings = 2.5; % while in tY we have [0,2*pi] in tX we have [0,2*pi*surroundings]
stepSize = pi/4; % Size of each step
steps = [3,3];
%
%
% Logfile
if getDebugLevel('logfile')
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
if getDebugLevel('LoadData') && exist([results,dataName,'-data.mat'],'file')
    load([results,dataName,'-data.mat']); %loads f and fn
    metaData = dir([results,dataName,'-data.mat']);
    debug('text',2,'Text',['Using loaded Data generated ',datestr(metaData.date),'.']);
    imgSize = [size(f,2),size(f,3)];
else
    f = ArtificialS2Data(pts,surroundings,steps,stepSize);
    fn = M.addNoise(f,sigma);
    debug('text',2,'Text',['Using _new_ Data generated ',datestr(datetime),'.']);
    if getDebugLevel('WriteData') %Write this version to file
        save([results,dataName,'-data.mat'],'f','fn','sigma','pts');
    end
end
%
%
% Plot original and noisy data
if getDebugLevel('Figures') % init meshgrids
    % Show in quiver (no color coding)
    [Xpl,Ypl] = meshgrid(1:pts);
    Zpl = zeros(size(Xpl));
end
if getDebugLevel('Figures')
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
if getDebugLevel('WriteImages')
    exportSphere2Asymptote(f,'File',[results,name,'-original.asy'],'ExportHeader',true);
    exportSphere2Asymptote(fn,'File',[results,name,'-noisy.asy'],'ExportHeader',true);
end
%
%
%% General Parameters
problem.M = M;
problem.f = fn;
problem.MaxIterations = 400;
problem.Epsilon = 0;
problem.lambda = pi/2;
problem.M.steps = 10; %10 gradient descent steps

alpha = [0 0.22];
beta = [0 29.5];
debug('text',2,'Text',['Parameter range alpha (',num2str(length(alpha)),' values): ',regexprep(num2str(alpha,5), '\s*', ','),'.']);
debug('text',2,'Text',['Parameter range beta  (',num2str(length(beta)),' values): ',regexprep(num2str(beta,5), '\s*', ','),'.']);

mDist = zeros(length(alpha)*length(beta),1);
mResults = zeros([size(f), length(alpha)*length(beta)]);
for i=1:length(alpha)*length(beta)
    [j1,j2] = ind2sub([length(alpha),length(beta)],i);    
    problem.alpha = alpha(j1);
    problem.beta = beta(j2);
    Lfr = cppa_ad_2D(problem);
    mResults(:,:,:,i) = Lfr;
    mDist(i) = 1/pts^2* sum(sum(problem.M.dist(f,Lfr)));
    if getDebugLevel('WriteImages')>0
        fileStr = [results,name,'-p-',num2str(alpha(j1)),'-',num2str(beta(j2))];
        fileStr(fileStr=='.') = [];
        exportSphere2Asymptote(Lfr,'File',[fileStr,'.asy'],'ExportHeader',true);
    end
    debug('text',1,'Text',...
        ['Parameters: \alpha=',num2str(alpha(j1)),' \beta=',num2str(beta(j2)),' yield ',num2str(mDist(i)),'.']);
end
%% Evaluate Results
[b,a] = (meshgrid(beta,alpha)); a = a(:); b = b(:);
% general minimum
[minValue,minIndex] = min(mDist);
[minJ1,minJ2] = ind2sub([length(alpha),length(beta)],minIndex);
    debug('text',1,'Text',['Minimum: Parameters: \alpha=',num2str(alpha(minJ1)),' \beta=',num2str(beta(minJ2)),' yields minimal value ',num2str(mDist(minIndex)),'.']);
minF = squeeze(mResults(:,:,:,minIndex));
% TV1 minimum
[minValueTV,minIndexTV] = min(mDist(b==0));
    debug('text',1,'Text',['Minimum TV: Parameter: \alpha=',num2str(alpha(minIndexTV)),' yields minimal value ',num2str(minValueTV),'.']);
minFTV = squeeze(mResults(:,:,:,minIndexTV));

% just TV2 minimum
[minValueTV2,minIndexTV2] = min(mDist(a==0));
    debug('text',1,'Text',['Minimum TV2: Parameter: \beta=',num2str(beta(minIndexTV2)),' yields minimal value ',num2str(minValueTV2),'.']);
minFTV2 = squeeze(mResults(:,:,:,minIndexTV2));

%%
if getDebugLevel('Figures')
    % Show in quiver (no color coding)
    figure;
    quiver3(Xpl,Ypl,Zpl,...
        permute(minF(1,:,:),[2,3,1]),...
        permute(minF(2,:,:),[2,3,1]),...
        permute(minF(3,:,:),[2,3,1])...'AutoScale','off'
    ); title(['''Optimal'' Reconstruction \alpha=',num2str(alpha(minJ1)),' \beta=',num2str(beta(minJ2)),' yields ',num2str(mDist(minIndex)),'.']); axis image; daspect([1,1,20/pts]); axis off
end
%% End logfile
if getDebugLevel('logfile')
    disp([' --- Logfile of Experiment ',name,'; ended ',datestr(datetime),' ---']);
    diary off;
end