function [dataStr,gD] = exportEllipses2Tikz(varargin)
% dataStr = exportEllipses2TikZ(f) - generate a standalone tex File for a plot
%          of given 1D, 2D of ellipses, i.e. 2-by-2 symmetric positive
%          definite matrices.
%
% INPUT
%   f : given (1D or 2D) data of symmetric positive definite matrices, i.e.
%       an 2x2xImageSize aray
%
% OPTIONAL
%   'Colors'          : ([]) provide a color data image, i.e. an array
%       of the same size of the data ImageSizex3 providing color data for
%       each data item.
%   'ColorScheme': 'GA' for the geometric anisotropiy index, 'FA' tor the
%                   euclidean anisotropy index,
%                  'Direction' for the directional encoding
%   'File': (String) Write the generated Code directly into a file
%   'ExportHeader' : (false) whether to export the Header of the .tex-File
%                    or just the tikz environment.
%   'GridDistance : (1/3) a value or vector of distances between the center
%       points of the ellipsoids. Standard Value is 1. If it is set to [],
%       the maximal axis is set to sqrt(3)*a, where a is the maximal eigenvalue
%       of all matrices involved.
%   'GridScale'  : (1) Scale the Grid additionally in the directions, which
%                  can be used to scale the optimal grid from the previous
%                  option down after its computation
%
% OUTPUT
%   dataStr : the String containing the tex source code
% ---
% ManImRes | R. Bergmann ~ 2015-12-05
%
ip = inputParser;
addRequired(ip,'f');
addParameter(ip,'Colors',[]);
addParameter(ip,'GridDistance',[]);
addParameter(ip,'GridScale',1);
addParameter(ip,'ColorScheme', 'GA');
addParameter(ip,'ExportHeader',false);
addParameter(ip,'File', '');
parse(ip, varargin{:});
vars = ip.Results;
dataDim = size(vars.f);
manDim = dataDim(1:2);
assert( all(manDim==[2,2]), 'visualization only works for 2x2 matrices');
if (numel(vars.Colors)>0)
    assert( all(size(vars.Colors)==[dataDim(3:end),3]),'The color data has to be the same as the data dimensions with an additional last dimension with RGB channels');
end
d = manDim(1);
dataDim = dataDim(3:end);
if length(vars.GridDistance)==1
    gD = vars.GridDistance*ones(1,2);
elseif isempty(vars.GridDistance)
    % Compute maximal eigenvalue ... may take some time?
    allM = reshape(vars.f,2,2,[]);
    maxEV = 0;
    for i=1:size(allM,3)
        maxEV = max(maxEV,max(eig(allM(:,:,i))));
    end
    gD = sqrt(length(dataDim))*maxEV*ones(1,2);
elseif length(vars.GridDistance) > (length(dataDim))
    error(['gridDistance vector is too long (',num2str(length(GridDistance)),'), a number or a vector of length ',num2str(length(dataDim)),' was expected.']);
else
    gD = ones(1,2);
    gD(1:length(vars.GridDistance)) = vars.GridDistance;
end
gD = gD.*vars.GridScale;
numCol = 1000;
pt = zeros(1,3);
if length(dataDim)>2 %not 2D or 1D
    error('only signals or images of ellipses possible');
end
if length(dataDim)>1
    pt(2) = dataDim(2);
end
pt(1) = dataDim(1);
% Header
if ~isempty(vars.File) %Export given
    exportFile = fopen(vars.File,'w','n','UTF-8');
end
dataStr = '';
if vars.ExportHeader
    dataStr = ['\\documentclass{standalone}\n',...
        '\\usepackage{pgfplots}\n',...
        '\\usepgfplotslibrary{colormaps}\n',...
        '\\pgfplotsset{compat=1.12}\n',...
        '\\tikzset{ellC/.style={/utils/exec={\\pgfplotscolormapdefinemappedcolor{#1}},%%\n',...
        '\tdraw=mapped color!80!black,very thin, fill=mapped color!80!white}}\n',...
        '\\begin{document}%%\n'];
end
dataStr = [dataStr,'\t\\pgfmathsetmacro{\\xStep}{',num2str(gD(1)),'}%%\n',...
    '\t\\pgfmathsetmacro{\\yStep}{',num2str(gD(2)),'}%%\n',...
    '\t\\begin{tikzpicture}\n',...
    '\t\t\\begin{axis}[hide axis,\n',...
    '\t\t\tcolormap/hsv\n,',...
    '\t\t\txmin=-1*\\xStep, xmax=',num2str(dataDim(1)+1),'*\\xStep,\n',...
    '\t\t\tymin=-1*\\yStep, ymax=',num2str(dataDim(2)+1),'*\\yStep,\n',...
    '\t\t\taxis equal]\n'];
if ~isempty(vars.File) %Export given
    fprintf(exportFile,dataStr);
end
for i=1:prod(dataDim) %run through all matrices
    [x,y] = ind2sub(dataDim,i); %position on grid, works for 1 & 2-dim data
    if all(all(~isnan(vars.f(:,:,x,y))))
        [V,D] = eig(vars.f(:,:,x,y));
        lambda = diag(D);
        if numel(vars.Colors)>0
            switch length(size(dataDim)>1)
                case 1
                    color = vars.Colors(x,:);
                case 2
                    color = vars.Colors(x,y,:);
                case 3
                    color = vars.Colors(x,y,z,:);
            end
        else
            switch vars.ColorScheme
                case 'GA'
                    if all(all(vars.f(:,:,x,y)==0))
                        colorInd=0;
                    else
                        AR = sqrt(((d-1)/d)*sum( log(lambda).^2) - 2/d*sum(sum(tril(log(lambda)*log(lambda)',-1))));
                        ARnormed = AR/(1+AR); %See Also Moakher & Batchelor in [0,1]
                        colInd = round(ARnormed*(numCol-1))+1; %color index
                        if ~isreal(colInd) || isnan(colInd) %might happen due to rounding errors
                            lambda = zeros(1,3);
                            colorInd = 0;
                        else
                            colorInd = colInd;
                        end
                    end
                case 'FA'
                    if all(all(vars.f(:,:,x,y)==0))
                        colorInd=0;
                    else
                        AF = sqrt(((d-1)/d)*sum( lambda.^2) - 2/d*sum(sum(tril(lambda*lambda',-1))));
                        AFnormed = AF/norm(vars.f(:,:,x,y),'fro'); %See Also Moakher & Batchelor in [0,1]
                        colInd = round(AFnormed*(numCol-1))+1; %color index
                        if ~isreal(colInd) || isnan(colInd) %might happen due to rounding errors
                            lambda = zeros(1,3);
                            colorInd = 0;
                        else
                            colorInd = colInd;
                        end
                    end
                    %            case 'Direction'
                    %                color = V(:,2) / max(abs(V(:,2)));
                    %                color = abs(color);
                otherwise
                    error('Unknown ColorScheme; available Schemes are ''GA'', ''FA''.');
            end
        end
        if ~all(all(vars.f(:,:,x,y)==0)) % only export nonzero ones for computational efforts
            alpha = atan2(V(2,1),V(1,1));
            line = ['\t\t\t\\draw[ellC=',num2str(colorInd),', rotate around={',num2str(alpha*180/pi),':(',num2str(y-1),'*\\xStep,',num2str(dataDim(1)-x),'*\\yStep)}] ',...
                '(',num2str(y-1),'*\\xStep,',num2str(dataDim(1)-x),'*\\yStep) ellipse (',num2str(lambda(2)),' and ',num2str(lambda(1)),');\n'];
            if ~isempty(vars.File) %Export given
                fprintf(exportFile,line);
            end
            if nargout>0
                dataStr = [dataStr,line];
            end
        end
    end
    %    if round(100*(i-1)/prod(dataDim))~=round(100*(i)/prod(dataDim))
    %        disp([' ',num2str(round(100*(i)/prod(dataDim))),'%']);
    %    end
end
endLine = '\t\t\\end{axis}\n\t\\end{tikzpicture}\n';
if vars.ExportHeader
    endLine = [endLine,'\\end{document}'];
end
if ~isempty(vars.File) %Export given
    fprintf(exportFile,endLine);
end
if nargout>0
    dataStr = sprintf([dataStr,endLine]);
end
if ~isempty(vars.File) %Export given
    fclose(exportFile);
end
if nargout==0
    clear dataStr;
end
end