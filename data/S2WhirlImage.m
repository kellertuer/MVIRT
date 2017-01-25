function f = S2WhirlImage(pts)
% An S2-dataset of several whirls on a slowly chaning background
% 
% Input
%   pts: (64) a number or 2 intergers sepcifying the image size
%
% Output
% f : an [3,pts,pts] or [3,pts(1),pts(2)] size array, the S2 valued image
%
% ---
% Manifold-valued Image Restoration Toolbox 1.1
%  J. Persch, R. Bergmann ~ 2017-01-05
% see LICENSE.txt
if nargin<1
    pts=64;
end 
pts = ones(1,2).*pts;
sur = [.5,.5];
x = linspace(0,sur(1)*2*pi,pts(1));
y = linspace(0,sur(2)*2*pi,pts(2));
[X,Y] = meshgrid(x,y);
f = zeros(3,pts(1),pts(2));
for i=1:pts(1)
    for j=1:pts(2)
        f(:,i,j) = ... % rotation for X coordinate in x,y-plane
        [   cos(X(i,j)+Y(i,j)), -sin(X(i,j)+Y(i,j)), 0; ...
            sin(X(i,j)+Y(i,j)), cos(X(i,j)+Y(i,j)), 0;
            0         ,           0, 1]*... %rot for Y: in xz plane
        [   cos(X(i,j)-Y(i,j)), 0, -sin(X(i,j)-Y(i,j)); ...
            0         , 1, 0;...
            sin(X(i,j) - Y(i,j)), 0, cos(X(i,j)-Y(i,j))]*[0;0;1]; %times north pole
    end
end
f(2,:,:) = 1;

pSize1 = floor(5.*pts./64);
pos1 = floor([16,11].*pts./64);
%f(:,10+(1:11),5+(1:11)) = -S2WhirlPatch(11);
f(:,pos1(1)+(-pSize1(1):pSize1(1)),pos1(2)+(-pSize1(2):pSize1(2))) = -S2WhirlPatch(2*pSize1+1);
pos2 = floor([41,58].*pts./64);
%f(:,35+(1:11),52+(1:11)) = -S2WhirlPatch(11);
f(:,pos2(1)+(-pSize1(1):pSize1(1)),pos2(2)+(-pSize1(2):pSize1(2))) = -S2WhirlPatch(2*pSize1+1);
pos3 = floor([11,41].*pts./64);
%f(:,5+(1:11),35+(1:11)) = S2WhirlPatch(11);
f(:,pos3(1)+(-pSize1(1):pSize1(1)),pos3(2)+(-pSize1(2):pSize1(2))) = S2WhirlPatch(2*pSize1+1);

pSize2 = floor(4.*pts./64);
pos4 = floor([35,7].*pts./64);
%f(:,30+(1:9),2+(1:9)) = S2WhirlPatch(9);
f(:,pos4(1)+(-pSize2(1):pSize2(1)),pos4(2)+(-pSize2(2):pSize2(2))) = S2WhirlPatch(2*pSize2+1);

pos5 = floor([25,41].*pts./64);
%f(:,20+(1:9),36+(1:9)) = S2WhirlPatch(9);
f(:,pos5(1)+(-pSize2(1):pSize2(1)),pos5(2)+(-pSize2(2):pSize2(2))) = S2WhirlPatch(2*pSize2+1);

pos6 = floor([32,25].*pts./64);
%f(:,27+(1:9),20+(1:9)) = -S2WhirlPatch(9);
f(:,pos6(1)+(-pSize2(1):pSize2(1)),pos6(2)+(-pSize2(2):pSize2(2))) = -S2WhirlPatch(2*pSize2+1);

pos7 = floor([7,60].*pts./64);
%f(:,2+(1:9),55+(1:9)) = S2WhirlPatch(9);
f(:,pos7(1)+(-pSize2(1):pSize2(1)),pos7(2)+(-pSize2(2):pSize2(2))) = S2WhirlPatch(2*pSize2+1);

pSize3 = floor(7.*pts./64);
pos8 = floor([23,56].*pts./64);
%f(:,15+(1:15),48+(1:15)) = S2WhirlPatch(15);
f(:,pos8(1)+(-pSize3(1):pSize3(1)),pos8(2)+(-pSize3(2):pSize3(2))) = S2WhirlPatch(2*pSize3+1);

pos9 = floor([38,45].*pts./64);
%f(:,30+(1:15),37+(1:15)) = -S2WhirlPatch(15);
f(:,pos9(1)+(-pSize3(1):pSize3(1)),pos9(2)+(-pSize3(2):pSize3(2))) = -S2WhirlPatch(2*pSize3+1);

pos10 = floor([16,28].*pts./64);
%f(:,8+(1:15),20+(1:15)) = -S2WhirlPatch(15);
f(:,pos10(1)+(-pSize3(1):pSize3(1)),pos10(2)+(-pSize3(2):pSize3(2))) = -S2WhirlPatch(2*pSize3+1);

pSize4 = floor(10.*pts./64);
pos11 = floor([51,16].*pts./64);
%f(:,40+(1:21),5+(1:21)) = +S2WhirlPatch(21);
f(:,pos11(1)+(-pSize4(1):pSize4(1)),pos11(2)+(-pSize4(2):pSize4(2))) = S2WhirlPatch(2*pSize4+1);

pSize5 = floor(8.*pts./64);
pos12 = floor([55,42].*pts./64);
%f(:,46+(1:17),33+(1:17)) = -S2WhirlPatch(17);
f(:,pos12(1)+(-pSize5(1):pSize5(1)),pos12(2)+(-pSize5(2):pSize5(2))) = -S2WhirlPatch(2*pSize5+1);

% normalize
f = f./repmat(sqrt(sum(f.^2,1)),3,1,1);
end
