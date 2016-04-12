function [Zn,Z] = ArtificalSARImage(pts,sigma)
% ArtificialSARImage(pts) constructs the artificial SAR image of size
%    pts x pts pixel
%
%   INPUT
%      pts : pixel size of the resulting Image
%      sigma : standard deviation of the wrapped Gaussian noise imposed on
%           the image.
%
%   OUTPUT
%      Zn : the artificial SAR Image (i.e. noisy
%      Z  : the noise free wrapped original
%
% This image is used in
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-07-23 | 2016-02-24
% see LICENSE.txt

[X,Y] = meshgrid(0:1/(pts-1):1,0:1/(pts-1):1);
Xm = X -1/2; Ym = Y-1/2;
% Rotate by 35 Degree
alpha = 35/180*pi;
Xr = cos(alpha)*Xm - sin(alpha)*Ym;
Yr = cos(alpha)*Ym + sin(alpha)*Xm;
% Build an Ellipse
axis = [1/6,1/25];
Zellipse = ((1/axis(1)*Xr.^2 + 1/axis(2)*Yr.^2) <= 1);
Zlinellipse = Zellipse.*(10*pi*Yr);

% Steps
% Divide in the middle
alpha = 45/180*pi;
Xr = cos(alpha)*Xm - sin(alpha)*Ym;
Yr = cos(alpha)*Ym + sin(alpha)*Xm;
Zstep = 4*pi*(Yr<0)-2*pi;

m = [0.275,0.275]; rad = 0.18;

Zsphere = -4*pi/(rad^2)*(((Xm-m(1)).^2 + (Ym-m(2)).^2) <= rad^2 ).*(rad^2-(Xm-m(1)).^2 - (Ym-m(2)).^2);

% Divide in the middle
alpha = 60/180*pi;
Xr = cos(alpha)*Xm - sin(alpha)*Ym;
Yr = cos(alpha)*Ym + sin(alpha)*Xm;
l = 0.075;
Zsq = zeros(size(Xr));
for i=1:8
    m = ([cos(alpha),-sin(alpha);cos(alpha),sin(alpha)]*[-.675+i*0.15,-.25]')';
    Zsq = Zsq + 2*pi*(i/8)*((abs(Xr-m(1))+abs((Yr-m(2)))) <= l);
end
Z = Zsphere-Zstep.*(1-Zellipse) + Zlinellipse + Zsq;
Zn = symMod(Z + sigma*randn(size(Z)),2*pi); %noisy signal
Z = symMod(Z,2*pi);
