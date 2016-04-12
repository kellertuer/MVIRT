function s = fig2TikZ(filename,t,f,ylim,yticks,yticklabels,xticks,xticklabels,scales)
% s = fig2TikZ(filename,t,f,ylim,ytsicks,ytickslabels,xticks,xticklabels,scales)
% Generate a TikZ-Figure from given data t (x-axis) and certain function
% values f. t automatically determines the x limits
% 
% INPUT
%   t      : Data points fot the x axis, in most cases a linspace vector
%   f      : one or more data terms
%   ylim   : limits for the y-axis [ymin,ymax]
%   yticks : Array containing strings (possibly incl math) of y-axis labels
%   xticks : Array containing strings (possibly incl math) of y-axis labels
%   scales : scale both axis (can also be 
%
% OUTPUT
%   s      : string containing the file content 
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-03-29
% see LICENSE.txt

numplots = size(f,1);
if (size(f,2) ~= length(t))
    error('t and f have to be same length');
end
yticknum = length(yticks);
xticknum = length(xticks);
tl = min(abs(max(ylim)-min(ylim)),abs(max(t)-min(t)));
xtickslength = 1/scales(2)*0.075*tl;
ytickslength = 1/scales(1)*0.075*tl;
% ytickslength = 1/scales(1)*0.01*(abs(max(t)-min(t))); %max(ylim)/max(t)

s = ['\begin{tikzpicture}[xscale=',num2str(scales(1)),', yscale=',num2str(scales(2)),']\n'];
% Axes
s = [s,'\t\draw[axis,->] (',num2str(min(t)),',',num2str(min(ylim)),') -- ('...
    ,num2str(max(t) + 1/2*1/xticknum*abs(max(t)-min(t))),',',num2str(min(ylim)),');\n'];
s = [s,'\t\draw[axis,->] (',num2str(min(t)),',',num2str(min(ylim)),') -- ('...
    ,num2str(min(t)),',',num2str(max(ylim) + 1/2*1/yticknum*abs(max(ylim)-min(ylim))),');\n'];
    % node[anchor=east] {$\\cdot 10^{',num2str(exp),'}$};\n'];


% ticks
for i=0:xticknum-1
    p = i/xticknum*(max(t)-min(t))+min(t);
 	s = [s, '\t\draw[ticks] (',num2str(xticks(i+1)),',',num2str(min(ylim)),'+',num2str(xtickslength),') -- (',...
		num2str(xticks(i+1)),',',num2str(min(ylim)),'-',num2str(xtickslength),') node[anchor=north] {$',...
		xticklabels{i+1},'$};\n'];
end
for i=0:yticknum-1
    p = i/yticknum*(max(ylim)-min(ylim))+min(ylim);
    s = [s, '\t\draw[ticks] (',num2str(min(t)),'+',num2str(ytickslength),',',num2str(yticks(i+1)),') -- (',...
        num2str(min(t)),'-',num2str(ytickslength),',',num2str(yticks(i+1)),') node[anchor=east] {$',...
        yticklabels{i+1},'$};\n'];
end

%grid
% for i=0:xticknum-1
%  	s = [s, '\t\\draw[grid] (',num2str(xticks(i+1)),',',num2str(min(ylim)),') -- (',...
% 		num2str(xticks(i+1)),',',num2str(max(ylim)),');\n'];
% end
% for i=0:yticknum-1
%  	s = [s, '\t\\draw[grid] (',num2str(min(t)),',',num2str(yticks(i+1)),') -- (',...
% 		num2str(max(t)),',',num2str(yticks(i+1)),');\n'];
% end

%separating lines
% for i=1:61
% 	if ismember(i,[2, 6, 14, 30,62])
% 	s = [s, '\t\\draw[dashed, gray,ultra thin] (',num2str(i-0.5),',',num2str(min(ylim)),') -- ('....
% 		num2str(i-0.5),',',num2str(max(ylim)),');\n'];
% 	end
% end

%s = [s,'\n'];
% plots
for i=1:numplots
    s = [s, '\draw[plot',num2str(i),']'];
    for j=1:length(t)
        if (j>1)
            s = [s,' --'];
        end
        s = [s, ' (',num2str(t(j),'%12.11f'),',',num2str(f(i,j),'%12.11f'),')'];
    end
    s = [s,';\n'];
end
s = [s '\end{tikzpicture}'];
fi = fopen(filename,'w'); 
fs = strrep(s,'\','\\');
fs = strrep(fs,'%','\%');
fs = strrep(fs,'\\t\','\t\');
fs = strrep(fs,'\\n\','\n\');
fprintf(fi,fs); 
fclose(fi);
end
