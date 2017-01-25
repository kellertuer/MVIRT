function fileStr = exportSphere2Asymptote(varargin)
% exportSphere2Asy(f) Export sphere-valued data to Asymptote
%
% Exports the given array of sphere valued data f (3xnxmxl) to an asympote
% file and additionally returns the content.
%
% INPUT
%   f : sphere-valued data (1D,2D or 3D)
%
% OUTPUT
%  fileStr : (if specified) String of the file contents. If no return value
%             is specified, this will not be returned
%
% OPTIONAL PARAMETERS
%   'File' : ('') file name of the file to export to: If not specified,
%             only the string will be generated.
%   'ColorMap'  : (parula(256)) Colormap to encode the elevation with
%   'ScaleAxes' : (1/3[1,1,1]) usually the spheres would not overlap, hence
%                 scaling the distances down makes the plot look like a
%                 colored quiver plot.
%   'ExportHeader' : (false) whether to export the Header of the .asy-File
%                    or just the drawing commands. 
%
% ---
% Manifold-Valued Image Restoration Toolbox 1.0, R. Bergmann ~ 2015-03-10 | 2015-03-27

    % Lofgile
    % 2015-03-27 Added the optional parameter ExportHeader
    ip = inputParser;
    addRequired(ip,'f');
    addParameter(ip,'File','');
    addParameter(ip,'ColorMap', parula(256));
    addParameter(ip,'ScaleAxes', [1/3,1/3,1/3]);
    addParameter(ip,'ExportHeader', true);
    parse(ip, varargin{:});
    vars = ip.Results;
    fe = reshape(vars.f,3,[]);
    minB = [1,1,1];
    sizes = size(vars.f);
    maxB = [1,1,1]; maxB(1:length(sizes(2:end))) = sizes(2:end);
    [Xpl,Ypl,Zpl] = meshgrid(minB(1):maxB(1),minB(2):maxB(2),minB(3):maxB(3));
    Xe = reshape(Xpl,1,[]);
    Ye = reshape(Ypl,1,[]);
    Ze = reshape(Zpl,1,[]);
    % Asymptote Export
    scales = vars.ScaleAxes;
    if ~isempty(vars.File)
        fID = fopen(vars.File,'w');
    end
    pt(1) = size(vars.f,2);
    pt(2) = size(vars.f,3);
    
    lfileStr = ''; %Header?
    if vars.ExportHeader
        lineS = sprintf(['import three;\nimport settings;\n',... %Packages
            'size(7cm);DefaultHead.size=new real(pen p=currentpen) {return 1.8mm;};\n',...%Arrow specification
         '\tcurrentprojection=perspective( camera=(',num2str(scales(1)*pt(1)/2),',',num2str(scales(2)*pt(2)/2),',',num2str(scales(3)*max(pt)),'),',...
         'up=Y, target = (',num2str(maxB(1)*scales(1)/2),',',num2str(maxB(2)*scales(2)/2),',0));']);
%            'currentprojection=perspective( camera=(0,0,10), up=Z,',...
%            'target = (',num2str(maxB(1)*scales(1)/2),',',num2str(maxB(2)*scales(2)/2),',0));']);
        lfileStr = [lfileStr,lineS];
        if ~isempty(vars.File)
            fprintf(fID,lineS);
        end
    end
    for i=1:size(fe,2)
        % Color encode elevation angle
        cInd = min(max( round((asin(...
            min(max(fe(3,i),-1),1)...
        )+pi/2)/pi*size(vars.ColorMap,1)),0),size(vars.ColorMap,1)-1)+1;
        thisC = vars.ColorMap(cInd,:);    
        lineS = ['draw((',num2str(scales(1)*Xe(i)),',',num2str(scales(2)*Ye(i)),',',num2str(scales(3)*Ze(i))...
            ,')--(',num2str(scales(1)*Xe(i)+fe(1,i)),',',num2str(scales(2)*Ye(i)-fe(2,i)),',',num2str(scales(3)*Ze(i)+fe(3,i))...
            ,'),rgb(',num2str(thisC(1)),',',num2str(thisC(2)),',',num2str(thisC(3)),'),Arrow3);\n'];
        if ~isempty(vars.File)
            fprintf(fID,lineS);
        end
        lfileStr = [lfileStr,lineS];
    end
    if ~isempty(vars.File)
        fclose(fID);
    end
    if nargout>0
        fileStr = lfileStr;
    end
end

