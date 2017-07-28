function [fhB,gD] = plotSPD(varargin)
% figurehandle = plotSPD(f,gridDist) - plot a given 1D or 2D signal of
%          symmetric positive definite 3x3 matrices.
%
% INPUT
%   f : given (1D, 2D or 3D) data of symmetric positive definite matrices, i.e.
%       an 3x3xImageSize-d array
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
%   'GridDistance'    : a value or vector of distances between the center
%       points of the ellipsoids. Standard Value is 1. If it is set to [],
%       the maximal axis is set to sqrt(3)*a, where a is the maximal
%       eigenvalue of all matrices involved.
%   'GridScale'       : (1) Scale the Grid additionally in the directions,
%       which can be used to scale the optimal grid from the previous
%       option down after its computation
%
% OUTPUT
%   figurehandle : the handle of the created figure
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-01-19 | 2015-12-05
% see LICENSE.txt

% Logfile
% 2015-03-27 Added the optional parameter ExportHeader
% 2015-12-05 added the colors Option
ip = inputParser;
addRequired(ip,'f');
addParameter(ip,'Colors',[]);
addParameter(ip,'ColorScheme', 'GA');
addParameter(ip,'EllipsoidPoints',6);
addParameter(ip,'File', '');
addParameter(ip,'GridDistance',[]);
addParameter(ip,'GridScale',1);
parse(ip, varargin{:});
vars = ip.Results;
dataDim = size(vars.f);
manDim = dataDim(1:2);
assert( all(manDim==[3,3]), 'visualization only works for 3x3 matrices');
if (numel(vars.Colors)>0)
    assert( all(size(vars.Colors)==[dataDim(3:end),3]),'The color data has to be the same as the data dimensions with an additional last dimension with RGB channels');
end
d = manDim(1);
dataDim = dataDim(3:end);
dataDim = [dataDim,ones(1,3-length(dataDim))];
if length(vars.GridDistance)==1
    gD = vars.GridDistance*ones(1,3);
elseif isempty(vars.GridDistance)
    % Compute maximal eigenvalue ... may take some time?
    allM = reshape(vars.f,3,3,[]);
    maxEV = 0;
    for i=1:size(allM,3)
        if all(reshape(~isnan(allM(:,:,i)),9,[])) && all(reshape(~isinf(allM(:,:,i)),9,[]))
            maxEV = max(maxEV,max(eig(allM(:,:,i))));
        end
    end
    gD = sqrt(length(dataDim))*maxEV*ones(1,3)/2;
else
    gD = ones(1,3);
    gD(1:length(vars.GridDistance)) = vars.GridDistance;
end
gD = gD.*vars.GridScale;
fh = gcf; % create empty figure if none exists, save handle
clf(fh);
hold on;
axis image
light('Position',[0 0 1],'Style','infinite');
numCol = 256;
hsvCM = hsv(numCol);

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
                        color=hsvCM(1,:);
                    else
                        AR = sqrt(((d-1)/d)*sum( log(lambda).^2) - 2/d*sum(sum(tril(log(lambda)*log(lambda)',-1))));
                        ARnormed = AR/(1+AR); %See Also Moakher & Batchelor in [0,1]
                        colInd = round(ARnormed*(numCol-1))+1; %color index
                        if ~isreal(colInd) || isnan(colInd) %might happen due to rounding errors
                            lambda = zeros(1,3);
                            color = hsvCM(1,:);
                        else
                            color = hsvCM(colInd,:);
                        end
                    end
                case 'FA'
                    if all(all(vars.f(:,:,x,y,z)==0))
                        color=hsvCM(1,:);
                    else
                        AF = sqrt(((d-1)/d)*sum( lambda.^2) - 2/d*sum(sum(tril(lambda*lambda',-1))));
                        AFnormed = AF/norm(vars.f(:,:,x,y,z),'fro'); %See Also Moakher & Batchelor in [0,1]
                        colInd = round(AFnormed*(numCol-1))+1; %color index
                        if ~isreal(colInd) || isnan(colInd) %might happen due to rounding errors
                            lambda = zeros(1,3);
                            color = hsvCM(1,:);
                        else
                            color = hsvCM(colInd,:);
                        end
                    end
                case 'Direction'
                    color = V(:,3) / max(abs(V(:,3)));
                    color = abs(color);
                otherwise
                    error('Unknown ColorScheme; available Schemes are ''GA'', ''FA'', and ''Direction''.');
            end
        end
        if ~all(all(vars.f(:,:,x,y,z)==0))
            [X,Y,Z]=ellipsoid(0,0,0,lambda(1),lambda(2),lambda(3),vars.EllipsoidPoints);
            for i1=1:size(X,1)
                for i2=1:size(X,2)
                    %rotate each coordinate
                    a =V*([X(i1,i2);Y(i1,i2);Z(i1,i2)]);
                    X(i1,i2) = a(1);
                    Y(i1,i2) = a(2);
                    Z(i1,i2) = a(3);
                end
            end
            % Shift
            X=X+x*gD(1);
            Y=Y+y*gD(2); % invert y
            Z=Z+y*gD(3);
            surf(X,Y,Z,'FaceColor',color,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        end
    end
end
hold off
axis off
view([0 90]);
fhB = fh;
end
