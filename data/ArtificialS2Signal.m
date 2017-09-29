function f = ArtificialS2Signal(pts,varargin)
% ArtificialS2Signal(pts) - Create a Lemniscate signal on the Sphere, i.e.
%      sample the Bernoullis Lemniscate (defined in the tangential space of
%      the north pole) at pts points
%
% INPUT:
%   pts : number of points to sample from
%
% OPTIONAL PARAMETERS
%   Interval : ([0,2*pi]) sample only a sub (or super) section of the
%              lemniscate defined on [0,2*pi] by specifying an
%              interval [a,b]
%   a        : (pi/(2*sqrt(2)) specify the size of the Lemniscate. Using
%              the standard, the lemniscate touches the equator twice,
%              cf. https://en.wikipedia.org/wiki/Lemniscate_of_Bernoulli
%  
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2016-10-07
% see LICENSE.txt

ip = inputParser;
addParameter(ip, 'Interval',[0,2*pi]);
addParameter(ip, 'a',pi/(2*sqrt(2)));
parse(ip, varargin{:});
par = ip.Results;
M = Sn(2);
%   Detailed explanation goes here
    %Create Data
    t = linspace(par.Interval(1), par.Interval(2),pts);
    a = par.a; % pi/(2*sqrt(2)); %Lemniscate of Bernoulli
    xc = a*sqrt(2)*(cos(t)./(sin(t).^2+1));
    yc = a*sqrt(2)*(cos(t).*sin(t)./(sin(t).^2+1));
    f = M.exp(repmat([0;0;1],[1,pts]),[xc;yc;zeros(size(xc))]); %length-keeping 
end

