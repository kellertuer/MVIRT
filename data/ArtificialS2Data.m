function f = ArtificialS2Data(pts, sur, steps, stepSize)
% f = ArtificialS2Data(pts,sur,jumps,jumpheight)
% Generate artificial spherical Data of size pts based on sur-surroundings
% of the equator an evelation parametrization and several steps (jumps) of
% height stepheight
%
% INPUT
%   pts   : a number (for a squared data array) or a vector of the size of
%           the resulting S2 data field
%   sur   : the generation is based on a parametrization of the elevation
%           running from the south pole over north again to south and a
%           surrounding of the equator. This vector denotes the numeberof
%           surroundings around the equator (first entry) and the number
%           of north-south-north surroundingd (second entry, std: 1)   
%   steps : number(s) of equidispaced positions in the image where jumps,
%           if one number is given, it is used both for x and y.
%   stepSize : size of the steps in x and y direction, may be given as a
%           number, e.g. pi/4 or a vector giving different step sizes for x
%           and y 
% OUTPUT
%      f  : example data field of size [3,pts(1),pts(2)]
% ---
% ManImRes 1.0, R. Bergmann ~ 2015-03-25

pts = ones(1,2).*pts;
steps = ones(1,2).*steps;
stepSize = ones(1,2).*stepSize;
if (length(sur)==1)
    sur = [sur,1];
end
% for the y component: divide 0,2pi into pts steps
    x = linspace(0,sur(1)*2*pi,pts(1));
    y = linspace(0,sur(2)*2*pi,pts(2));
    [X,Y] = meshgrid(x,y);
    f = zeros(3,pts(1),pts(2));
    for i=1:pts(1)
        for j=1:pts(2)
            thisStepX = stepSize(1)*floor(i/pts(1)*steps(1));
            tX = X(i,j) + thisStepX;
            thisStepY = stepSize(2)*floor(j/pts(2)*steps(2));
            tY = Y(i,j) + thisStepY;
            f(:,i,j) = ... % rotation for X coordinate in x,y-plane
            [   cos(tX+tY), -sin(tX+tY), 0; ...
                sin(tX+tY), cos(tX+tY), 0;
                0         ,           0, 1]*... %rot for Y: in xz plane
            [   cos(tX-tY), 0, -sin(tX-tY); ...
                0         , 1, 0;...
                sin(tX-tY), 0, cos(tX-tY)]*[0;0;1]; %times north pole
        end
    end
end

