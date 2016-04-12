function [fn,fo,t] = CyclicPiecewiseSignal(n,sigma)
% [fn,fo,t] = CyclicPiecewiseSignal(n,sigma)
%    Create a signal on [0,1] that gets steeper and steeper towards 1.
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
% This image is used in
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann (2014). 
%       Second order differences of cyclic data and applications in variational denoising. 
%       SIAM Journal on Imaging Sciences. 7, (4), 2916?2953.
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2015-06-19 | 2016-02-24
% see LICENSE.txt

    t = linspace(0,1,n); % on the interval [0,1]
    a = 4;
    fo = zeros(1,length(t));
    fo( t<=1/6) = -pi/2+a*pi/8*t(t<=1/6);
    % first term shifts the function piece in order to be continuous
    fo( (1/6<t)&(t<=1/3) ) = max(fo(fo~=0)) - a*pi/4*1/6 + a*pi/4*t( (1/6<t)&(t<=1/3) );
    fo( (1/3<t)&(t<=1/2) ) = max(fo(fo~=0)) - a*pi/2*1/3 + a*pi/2*t( (1/3<t)&(t<=1/2) );
    fo( (1/2<t)&(t<=2/3) ) = max(fo(fo~=0)) - a*pi*1/2   + a*pi*t  ( (1/2<t)&(t<=2/3) );
    fo( (2/3<t)&(t<=5/6) ) = max(fo(fo~=0)) - a*2*pi*2/3 + a*2*pi*t( (2/3<t)&(t<=5/6) );
    fo( 5/6<t ) =            max(fo(fo~=0)) - a*4*pi*5/6 + a*4*pi*t( 5/6<t );
    fo = symMod(fo,2*pi);
     
    fn = fo + sigma*randn(size(fo)); %noisy signal
    fn = symMod(fn,2*pi); % modulo 2pi -> wrapped Gaussian noise
    fo = symMod(fo,2*pi); % modulo 2pi, original signal
end

