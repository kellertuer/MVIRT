%
% S2-Lemniscate - Median and mean
% =============
%
% Illustrates the functions of the median and mean on a manifold for a
% set of S2Signals
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2016-10-07
% see LICENSE.txt

close all; clc
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../../initMVIRT.m')

setDebugLevel('LevelMin',0);     % Lower bound
setDebugLevel('LevelMax',3);     % Upper bound
setDebugLevel('text',2);    % Not so much text debug
setDebugLevel('time',3);    % Times

s1 = ArtificialS2Signal(129,'Interval',[0,pi/2]);
s2 = ArtificialS2Signal(129,'Interval',[pi/2,pi]);
s3 = ArtificialS2Signal(129,'Interval',[pi,3*pi/2]);
s4 = ArtificialS2Signal(129,'Interval',[3*pi/2,2*pi]);

%visual sphere data
X = cell(1,3);
[X{:}] = sphere(40);
lightGrey = 0.8*[1 1  1]; % It looks better if the lines are lighter

hS=surface(X{:},'FaceColor', 'none','EdgeColor',lightGrey,'LineWidth',.1)
hold on
h1 = plot3(s1(1,:),s1(2,:),s1(3,:),'o','Markersize',3) 
h2 = plot3(s2(1,:),s2(2,:),s2(3,:),'o','Markersize',3) 
h3 = plot3(s3(1,:),s3(2,:),s3(3,:),'o','Markersize',3) 
h4 = plot3(s4(1,:),s4(2,:),s4(3,:),'o','Markersize',3) 


M = Sn(2);
data = cat(2,permute(s1,[1,3,2]), permute(s2,[1,3,2]), ...
    permute(s3,[1,3,2]), permute(s4,[1,3,2]));
means = M.mean(data);
hMe=plot3(means(1,:),means(2,:),means(3,:),'x','Markersize',5) 

medians = M.median(data);
hMe2=plot3(medians(1,:),medians(2,:),medians(3,:),'x','Markersize',5) 

hold off
legend([h1,h2,h3,h4,hMe,hMe2],{'data 1', 'data 2', 'data 3', 'data 4', 'means', 'medians'});
title('4 Signal parts on S2.');


