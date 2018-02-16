function x = ArtificialS2SignalLemniscate(varargin)
% artificialS2SignalLemniscate(pts)
% create a curve in the north pole tangent space in the form of the
% Lemniscate of Benoulli sampled at pts point.
%
% INPUT
%   pts : number of points to sample on
%
% OPTIONAL
%   'a'  : (pi/2*sqrt(2)) parameter of the Lemniscate, the standard value
%           yields a curve raaching down to the equa
%   'range'
ip = inputParser;
addRequired(ip,'pts');
addParameter(ip,'a',pi/(2*sqrt(2)));
addParameter(ip,'range',[0,2*pi]);
parse(ip, varargin{:});
vars = ip.Results;
M = Sn(2);
pts = vars.pts;
t = linspace(vars.range(1),vars.range(2),pts);
a = vars.a; %Create a Lemniscate of Bernoulli
xc = a*sqrt(2)*(cos(t)./(sin(t).^2+1));
yc = a*sqrt(2)*(cos(t).*sin(t)./(sin(t).^2+1));
zc = ones(size(xc)); % i.e. in T_pS^2, where p is the north pole.
x = M.exp(repmat([0;0;1],[1,pts]),[xc;yc;zeros(size(xc))]); %get down to s2
end

