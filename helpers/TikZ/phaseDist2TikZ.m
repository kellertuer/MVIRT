function s = phaseDist2TikZ(f,r,varargin)
% phaseDist2TikZ(f,R)
% Draw given data f on a circle and the complex number R as a vector.
% If mor than one data (column) f is given, they are plotted on different
% circles
%
%   INPUT
%       f : column(s) of phase data
%       R : a complex number for each column (length(f)) of f, representing
%               a certain mean.
%
%   OPTIONAL ARGUMENTS (standard value)
%       'Scale'       : (1) scale of the TikZ-Picture (1)
%       'Radius'      : (0.15) increase of the radii for data f starting
%                           from 'Raduis'/2+1. Can also be specified as a 
%                           vector of length(f).
%       'Markersize'  : (0.125) Size of the Marker denoted at each data
%                           point of f
%       'File'        : ('') specify a file to save the TikZ figure to.
%
%   OUTPUT
%       s : a string containing the TikZ-figure;
%           (formated string with escaped \, \n, and \t which are _not_ escaped)
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-04-24 
% see LICENSE.txt

[n,m] = size(f); %m is the number of data sets we have, n their length
if length(r)~=m
    error('There has to be an r for each data f');
end
ip = inputParser;
addParamValue(ip, 'Scale', 1);
addParamValue(ip, 'Radius', 0.125);
addParamValue(ip, 'MarkerSize', 1.25);
addParamValue(ip, 'File', '');
parse(ip, varargin{:});
par = ip.Results;
if sum(size(par.Radius))==2 %number -> increase from 1+r/2 to that number
    radius = 1+par.Radius/2+linspace(0,(m-1)*par.Radius,m);
elseif isvector(par.Radius) && (length(par.Radius)==m)
    radius = par.Radius;
else
    error('Wrong format of the optional radius parameter.');
end
    s = '\tikzstyle{unitcircle}=[thin,gray]';
    for i=1:m
        s = [s,'% Data #',num2str(i),'\n'];        
        s = [s,'\tikzstyle{circle',num2str(i),'}=[black!',num2str(100*(1-(i-1)/m)),']\n'];
        s = [s, '\tikzstyle{plot',num2str(i),'}=[mark=o,only marks, mark size=',num2str(par.MarkerSize,'%8.7f'),', mark options={color=black!',num2str(100*(1-(i-1)/m)),'}]\n'];
        s = [s, '\tikzstyle{dataR',num2str(i),'}=[thick,draw=red!',num2str(100*(1-(i-1)/m)),']\n'];
    end
    s = [s,'\begin{tikzpicture}[scale=',num2str(par.Scale,'%8.7f'),']\n',...
		'\t\draw[unitcircle] (0,0) circle (1);\n'];    
    for i=1:m
        %Draw circles of the radii
        s = [s,'\t\draw[circle',num2str(i),'] (0,0) circle (',num2str(radius(i),'%8.7f'),');\n'];
        s = [s,'\t\draw plot[plot',num2str(i),'] coordinates {'];
        for j=1:n
            s = [s,'(',num2str(180/pi*(f(j,i)+pi),'%8.7f'),':',num2str(radius(i),'%8.7f'),') '];
        end
        s = [s, '};\n'];
        s = [s, '\draw[dataR',num2str(i),'] (0,0) -- (',...
            num2str(180/pi*(angle(r(i))+pi),'%8.7f'),':',...
            num2str(abs(r(i)),'%8.7f'),') node [pos=0.5, above,font=\tiny,rotate=',...
            num2str(180/pi*(angle(r(i))+pi), '%8.7f'),'] {',...
            num2str(abs(r(i)),'%8.7f'),'};\n'];  
    end
    s = [s,'\end{tikzpicture}'];    
    if numel(par.File)>0
        fi = fopen(par.File,'w'); 
        fs = strrep(s,'\','\\');
        fs = strrep(fs,'\\t\','\t\');
        fs = strrep(fs,'\\n\','\n\');
        fs = strrep(fs,'\\n%','\n%');
        fs = strrep(fs,'%','%%');
        fprintf(fi,fs); 
        fclose(fi);
    end
end

