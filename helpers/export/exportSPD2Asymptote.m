function [dataStr,gD] = exportSPD2Asymptote(varargin)
% dataStr = exportSPD2Asymptote(f) - generate a File for a POVRay plot a given
%          1D, 2D or 3D signal of symmetric positive definite matrices.
%
% INPUT
%   f : given (1D or 2D) data of symmetric positive definite matrices, i.e.
%       an 3x3xImageSize-d aray
%
% OPTIONAL
%   'Colors'          : ([]) provide a color data image, i.e. an array
%       of the same size of the data ImageSizex3 providing color data for
%       each data item.
%   'ColorScheme'     : ('GA') values to be used for coloring:
%       'GA' [G]eometric [A]nisotropy, 'FA' [F]ractional [A]nisotropy,
%       'Direction' for the directional encoding. This values is ignored
%       if a color image (see previous Option) is provided
%   'EllipsoidPoints' : (6) number of points used to surf eny ellipse. The
%       standard value does not produce that nice results, but scales up to
%       at least 112^2; for smaller images use a larger value, e.g. 15.
%   'File'            : (String) Write the generated Code directly into a file
%   'ExportHeader'    : (false) whether to export the Header of the
%       .asy-File or just the drawing commands.
%   'GridDistance'    : a value or vector of distances between the center
%       points of the ellipsoids. Standard Value is 1. If it is set to [],
%       the maximal axis is set to sqrt(3)*a, where a is the maximal
%       eigenvalue of all matrices involved.
%   'GridScale'       : (1) Scale the Grid additionally in the directions,
%       which can be used to scale the optimal grid from the previous
%       option down after its computation
%
% OUTPUT
%   dataStr : the String containing the POV-Ray source code
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-02-24 | 2015-12-05

% Logfile
% 2015-03-27 Added the optional parameter ExportHeader
% 2015-04-01 Added GridScale
ip = inputParser;
addRequired(ip,'f');
addParameter(ip,'Colors',[]);
addParameter(ip,'GridDistance',[]);
addParameter(ip,'GridScale',1);
addParameter(ip,'ColorScheme', 'GA');
addParameter(ip,'ExportHeader',false);
addParameter(ip,'File', '');
addParameter(ip,'ColorMap',hsv);
parse(ip, varargin{:});
vars = ip.Results;
dataDim = size(vars.f);
manDim = dataDim(1:2);
assert( all(manDim==[3,3]), 'visualization only works for 3x3 matrices');
d = manDim(1);
if (numel(vars.Colors)>0)
    assert( all(size(vars.Colors)==[dataDim(3:end),3]),'The color data has to be the same as the data dimensions with an additional last dimension with RGB channels');
end
dataDim = dataDim(3:end);
if length(vars.GridDistance)==1
    gD = vars.GridDistance*ones(1,3);
elseif isempty(vars.GridDistance)
    % Compute maximal eigenvalue ... may take some time?
    allM = reshape(vars.f,3,3,[]);
    maxEV = 0;
    for i=1:size(allM,3)
        maxEV = max(maxEV,max(eig(allM(:,:,i))));
    end
    gD = sqrt(length(dataDim))*maxEV*ones(1,3);
else
    gD = ones(1,3);
    gD(1:length(vars.GridDistance)) = vars.GridDistance;
end
gD = gD.*vars.GridScale;
numCol = size(vars.ColorMap,1);%256;
CM = vars.ColorMap;
pt = zeros(1,3);
if length(dataDim)>2 %2D or 1D
    pt(3) = dataDim(3);
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
    dataStr = sprintf([...%'\\\\begin{asypicture}{name=',vars.File,'}\n',...
        '\timport three;\n\timport settings;\n\n',...
        '\tsurface ellipsoid(triple v1,triple v2,triple v3,real l1,real l2, real l3, triple pos=O) {\n',...
        '\t\ttransform3 T = identity(4);\n\t\tT[0][0] = l1*v1.x;\n\t\tT[1][0] = l1*v1.y;\n\t\tT[2][0] = l1*v1.z;\n',...
        '\t\tT[0][1] = l2*v2.x;\n\t\tT[1][1] = l2*v2.y;\n\t\tT[2][1] = l2*v2.z;\n',...
        '\t\tT[0][2] = l3*v3.x;\n\t\tT[1][2] = l3*v3.y;\n\t\tT[2][2] = l3*v3.z;\n',...
        '\t\tT[0][3] = pos.x;\n\t\tT[1][3] = pos.y;\n\t\tT[2][3] = pos.z;\n',...
        '\t\treturn T*unitsphere;\n\t}\n\n',...
        '\tsize(200);\n\n',...
        '\treal gDx=',num2str(gD(1)),';\n',...
        '\treal gDy=',num2str(gD(2)),';\n',...
        '\treal gDz=',num2str(gD(3)),';\n\n',...
        '\tcurrentprojection=perspective( camera=(gDx*',num2str((pt(1)-1)/2),',gDy*',num2str((pt(2)-1)/2),',',num2str(max(dataDim)),'),',...
        'up=Y, target=(gDx*',num2str((pt(1)-1)/2),',gDy*',num2str((pt(2)-1)/2),',gDz*',num2str((pt(3)-1)/2),'));\n']);
    %         '\tcurrentprojection=perspective( camera=(gDx*',num2str(pt(1)/2),',gDy*',num2str(pt(2)/2),',-',num2str(max(dataDim)),'-2*gDx),',...
    %         'up=Z, target=(gDx*',num2str(pt(1)/2),',gDy*',num2str(pt(2)/2),',gDz*',num2str(pt(3)/2),'));']);
end
if ~isempty(vars.File) %Export given
    fprintf(exportFile,dataStr);
end
for i=1:prod(dataDim) %run through all matrices
    [x,y,z] = ind2sub(dataDim,i); %position on grid, works for all 3 dimension (1,2,3 dim data)
    if all(all(~isnan(vars.f(:,:,x,y,z))))
        [V,D] = eig(vars.f(:,:,x,y,z));
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
                if all(all(vars.f(:,:,x,y,z)==0))
                    color=CM(1,:);
                else
                    AR = sqrt(((d-1)/d)*sum( log(lambda).^2) - 2/d*sum(sum(tril(log(lambda)*log(lambda)',-1))));
                    ARnormed = AR/(1+AR); %See Also Moakher & Batchelor in [0,1]
                    colInd = round(ARnormed*(numCol-1))+1; %color index
                    if ~isreal(colInd) || isnan(colInd) %might happen due to rounding errors
                        lambda = zeros(1,3);
                        color = CM(1,:);
                    else
                        color = CM(colInd,:);
                    end
                end
            case 'FA'
                if all(all(vars.f(:,:,x,y,z)==0))
                    color=CM(1,:);
                else
                    AF = sqrt(((d-1)/d)*sum( lambda.^2) - 2/d*sum(sum(tril(lambda*lambda',-1))));
                    AFnormed = AF/norm(vars.f(:,:,x,y,z),'fro'); %See Also Moakher & Batchelor in [0,1]
                    colInd = round(AFnormed*(numCol-1))+1; %color index
                    if ~isreal(colInd) || isnan(colInd) %might happen due to rounding errors
                        lambda = zeros(1,3);
                        color = CM(1,:);
                    else
                        color = CM(colInd,:);
                    end
                end
            case 'Direction'
                color = V(:,3) / max(abs(V(:,3)));
                color = abs(color);
            otherwise
                error('Unknown ColorScheme; available Schemes are ''GA'', ''FA'', and ''Direction''.');
        end
        end
        if ~all(all(vars.f(:,:,x,y,z)==0)) % only export nonzero ones for computational efforts
            line = ['\t\tdraw(\n\t\t\tellipsoid( ',...
                '(',num2str(V(1,1)),',',num2str(V(2,1)),',',num2str(V(3,1)),'), ',...
                '(',num2str(V(1,2)),',',num2str(V(2,2)),',',num2str(V(3,2)),'), ',...
                '(',num2str(V(1,3)),',',num2str(V(2,3)),',',num2str(V(3,3)),'), ',...
                num2str(lambda(1)),', ',num2str(lambda(2)),', ',num2str(lambda(3)),', ',...
                '(gDx*',num2str(x-1),',gDy*',num2str(y-1),',gDz*',num2str(z-1),')),',...
                'rgb(',num2str(color(1)),',',num2str(color(2)),',',num2str(color(3)),'));\n'];
            if ~isempty(vars.File) %Export given
                fprintf(exportFile,line);
            end
            if nargout>0
                dataStr = sprintf([dataStr,line]);
            end
        end
    end
%    if round(100*(i-1)/prod(dataDim))~=round(100*(i)/prod(dataDim))
%        disp([' ',num2str(round(100*(i)/prod(dataDim))),'%']);
%    end
end
if ~isempty(vars.File) %Export given
    fclose(exportFile);
end
if nargout==0
    clear dataStr;
end
end