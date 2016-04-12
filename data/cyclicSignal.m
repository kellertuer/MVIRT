function [fn,fo, t] = cyclicSignal(n,sigma)
% cyclicSignal(n,sigma);
%  generate a S1-valued signal with n data points and a wrapped gaussian
%  noise with sigma standard deviation.
%
% INPUT
%     n     : length of the signal
%     sigma : amount of noise, standard deviation of the wrapped gaussian
%
% OUTPUT
%     fn : the generated (noisy) signal
%     fo : the generated signal without noise to compare with
%     t  : corresponding x-axis values for both signals on [0,1]
%
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-11-29 | 2016-02-24
% see LICENSE.txt

t = linspace(0,1,n); % on the interval [0,1]
fo = zeros(1,length(t));
%Define segments
fo(t<=1/4) = -24*pi*(t(t<=1/4)-1/4).^2 + 3*pi/4; %Parabel
fo((t>1/4)&(t<=3/8)) = 4*pi*t((t>1/4)&(t<=3/8)) - pi/4; %liear
fo((t>3/8)&(t<=1/2)) = -pi*t((t>3/8)&(t<=1/2)) -3/8*pi;
fo((t>1/2)&(t<=19/32)) = -7/8*pi;
fo((t>19/32)&(t<=11/16)) = -8/8*pi;
fo((t>11/16)&(t<=25/32)) = 7/8*pi;
fo((t>25/32)&(t<=7/8)) = 6/8*pi;
fo((t>7/8)&(t<=1)) = 3/2*pi/exp(8/7-1/((1-7/8))) * exp(8/7-1./(1-t((t>7/8)&(t<=1))) ) - 3/4*pi;

fn = fo + sigma*randn(size(fo)); %noisy signal
fn = symMod(fn,2*pi); % modulo 2pi -> wrapped Gaussian noise
fo = symMod(fo,2*pi); % modulo 2pi, original signal
