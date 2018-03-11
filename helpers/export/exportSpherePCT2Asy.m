function fileStr = exportSpherePCT2Asy(varargin)
% fileStr = exportSphereSignals2Asymptote(pts,curves,xi,colors)
% Export a onedimensional signal of spherical data to points on the sphere
%
% INPUT
%   pts    : data of 3xl in s cells, l data points on the sphere embedded in R3
%   curves : data cell array of points, t cells, on S2 to be plotted as curves
%   xi     : 3x2xm arrays in cells, u cells 
%   colors : 3x(s+t+u), i.e. columns of color values for the sets
% OUTPUT
%   fileStr : (if specified) string containing the file contents
%
% OPTIONAL PARAMETERS
%   'DotSize': ('1pt') dot size for the data points
%   'File' : ('') file name of the file to export to: If not specified,
%             only the string will be generated.
%   'ExportHeader' : (true) whether to export the Header of the .asy-File
%                    or just the drawing commands.
%   'OpacityVector' : ([0.5,1,...,1], s+t+u entries Make some signals opaque.
%                     The standard sets all to 1.
%
% ---
% Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-03-26 | 2015-04-08

% Lofgile
% 2015-04-08 Added the optional parameter OpacityVector
% 2015-03-27 Added the optional parameter ExportHeader

ip = inputParser;
addRequired(ip,'pts');
addRequired(ip,'curves');
addRequired(ip,'xi');
addRequired(ip,'colors');
addParameter(ip,'File','');
addParameter(ip,'DotSize',6);
addParameter(ip,'ArrowHead',6);
addParameter(ip,'OpacityVector',[]);
addParameter(ip,'Camera',[1,1,0]);
addParameter(ip,'ExportHeader', true);
parse(ip, varargin{:});
vars = ip.Results;
s = length(vars.pts);
t = length(vars.curves);
u = length(vars.xi);
m = s+t+u;
if isscalar(vars.DotSize)
    dS = vars.DotSize.*ones([1,m]);
else
    dS = vars.DotSize;    
end
if length(vars.colors) ~= m
    error(['Wrong color dimensions ',num2str(size(vars.colors,2)),' (data requires',num2str(m),')']);
end
if ~isempty(vars.File)
    fID = fopen(vars.File,'w');
end
if isempty(vars.OpacityVector)
    vars.OpacityVector = ones(1,s+t+u);
end
lfileStr = ''; %Header?
if vars.ExportHeader
    header = sprintf(['import settings;\nimport three;\nimport solids;',...
        'unitsize(4cm);\n\n',...
        'currentprojection=perspective( camera=(',num2str(vars.Camera(1)),...
        ',',num2str(vars.Camera(2)),',',num2str(vars.Camera(3)),'), target = (0,0,0));',...
        'currentlight=nolight;\n\n',...
        'revolution S=sphere(O,1);\n',...
        'draw(surface(S),surfacepen=lightgrey+opacity(.6), meshpen=0.6*white+linewidth(.5pt));',...
    ]);
    lfileStr = [lfileStr,header];
    if ~isempty(vars.File)
        fprintf(fID,header);
    end
end
for i=1:m %linestyles
    lColor = vars.colors{i};
    line = sprintf(['\npen LineStyle',num2str(i),' = ',...
        'rgb(',num2str(lColor(1)),',',num2str(lColor(2)),',',num2str(lColor(3)),')+',...
        'linewidth(',num2str(dS(i)),')+opacity(',num2str(vars.OpacityVector(i)),');']);
    if ~isempty(vars.File)
        fprintf(fID,line);
    end
    lfileStr = [lfileStr,line]; %#ok<*AGROW>
end
% for the signals, 1,...,s
for i=1:s %for each signal
    lSig = vars.pts{i};
    l = size(lSig,2);
    for j=1:l %for all dots each
        line = sprintf(['dot((',num2str(lSig(1,j)),',',num2str(lSig(2,j)),',',num2str(lSig(3,j)),'),LineStyle',num2str(i),');\n']);
        if ~isempty(vars.File)
            fprintf(fID,line);
        end
        lfileStr = [lfileStr,line];
    end
end
% for curves
for i=1:t
    lCurve = vars.curves{i};
    l = size(lCurve,2);
    line = ['path3 p',num2str(i),' = '];
    for j=1:l %for all dots each
        if j > 1
            line = sprintf([line,' .. ']);
        end
        line = sprintf([line,'(',num2str(lCurve(1,j)),',',num2str(lCurve(2,j)),',',num2str(lCurve(3,j)),')']);
    end
    line = sprintf([line,';\ndraw(p',num2str(i),',LineStyle',num2str(s+i),');\n']);
    if ~isempty(vars.File)
        fprintf(fID,line);
    end
    lfileStr = [lfileStr,line];
end
% for tangent vectors
for i=1:u
    lxi = vars.xi{i};
    l = size(lxi,3);
    for j=1:l
        lineS = sprintf(['draw((',num2str(lxi(1,1,j)),',',num2str(lxi(2,1,j)),',',num2str(lxi(3,1,j))... %base vec
            ,')--(',num2str(lxi(1,1,j)+lxi(1,2,j)),',',num2str(lxi(2,1,j)+lxi(2,2,j)),',',num2str(lxi(3,1,j)+lxi(3,2,j))...
            ,'),LineStyle',num2str(s+t+i),',Arrow3(',num2str(vars.ArrowHead),'));\n']);
         if ~isempty(vars.File)
            fprintf(fID,lineS);
         end
        lfileStr = [lfileStr,lineS];
    end
end
if ~isempty(vars.File)
    fclose(fID);
end
if nargout>0
    fileStr = lfileStr;
end
