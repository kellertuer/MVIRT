function fileStr = exportSphereSignals2Asymptote(varargin)
% fileStr = exportSphereSignals2Asymptote(f,colors)
% Export a onedimensional signal of spherical data to points on the sphere
%
% INPUT
%   f : data of 3xlxs, i.e. s sets of l data points on the sphere embedded in R3
%   colors : 3xs, i.e. columns of color values for the sets
% OUTPUT
%   fileStr : (if specified) string containing the file contents
%
% OPTIONAL PARAMETERS
%   'DotSize': ('1pt') dot size for the data points
%   'File' : ('') file name of the file to export to: If not specified,
%             only the string will be generated.
%   'ExportHeader' : (false) whether to export the Header of the .asy-File
%                    or just the drawing commands.
%   'OpacityVector' : ([0.5,1,...,1] Make some signals opaque. The standard
%                     sets the first signal to 1/2 the others to 1.
%
% ---
% Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-03-26 | 2015-04-08

% Lofgile
% 2015-04-08 Added the optional parameter OpacityVector
% 2015-03-27 Added the optional parameter ExportHeader

ip = inputParser;
addRequired(ip,'f');
addRequired(ip,'colors');
addParameter(ip,'File','');
addParameter(ip,'DotSize','1.5pt');
addParameter(ip,'OpacityVector',[]);
addParameter(ip,'ExportHeader', false);
parse(ip, varargin{:});
vars = ip.Results;
if (size(vars.f,1) ~= 3) || (length(size(vars.f)))>3
    error('onedimensional signal of data on S2 in R3 required)');
end
[~,l,s] = size(vars.f);
if any(size(vars.colors) ~= [3,s])
    error(['Wrong color dimensions ',num2str(size(vars.colors)),' (data requires',num2str([3,s]),')']);
end
if ~isempty(vars.File)
    fID = fopen(vars.File,'w');
end
if isempty(vars.OpacityVector)
    vars.OpacityVector = ones(1,s);
    vars.OpacityVector(1) = 0.5;
end
lfileStr = ''; %Header?
if vars.ExportHeader
    header = sprintf(['import settings;\nimport three;\nimport solids;',...
        'unitsize(4cm);\n\n',...
        'currentprojection=perspective( camera=(1,.4,.9), target = (0,0,0));',...
        'currentlight=nolight;\n\n',...
        'revolution S=sphere(O,1);\n',...
        'draw(surface(S),surfacepen=lightgrey+opacity(.6), meshpen=0.6*white+linewidth(.5pt));',...
    ]);
    lfileStr = [lfileStr,header];
    if ~isempty(vars.File)
        fprintf(fID,header);
    end
end
for i=1:s %linestyles
    line = sprintf(['\npen LineStyle',num2str(i),' = ',...
        'rgb(',num2str(vars.colors(1,i)),',',num2str(vars.colors(2,i)),',',num2str(vars.colors(3,i)),')+',...
        'linewidth(',num2str(vars.DotSize),')+opacity(',num2str(vars.OpacityVector(i)),');']);
    if ~isempty(vars.File)
        fprintf(fID,line);
    end
    lfileStr = [lfileStr,line]; %#ok<*AGROW>
end
for i=1:s %for each signal
    for j=1:l %for all dots each
        line = sprintf(['dot((',num2str(vars.f(1,j,i)),',',num2str(vars.f(2,j,i)),',',num2str(vars.f(3,j,i)),'),LineStyle',num2str(i),');\n']);
        if ~isempty(vars.File)
            fprintf(fID,line);
        end
        lfileStr = [lfileStr,line];
    end
end
if ~isempty(vars.File)
    fclose(fID);
end
if nargout>0
    fileStr = lfileStr;
end
