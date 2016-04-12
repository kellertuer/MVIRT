function p = symMod(v,T)
% symMod(v,T) compute the modulus with respect to [-T/2, T/2)
%
%   INPUT
%       v : value(s) of arbitrary data shape
%       T : period T
%
%   OUTPUT
%       p : the congruence class representant w.r.t. [-T/2, T/2)
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-02-12 | 2014-11-29
% see LICENSE.txt

%
% Change Log
% 2014-11-29 Adopted an optimization coined by J. Persch to spare half the
%            time by avoiding mod.

% p = mod(v+T/2,T)-T/2;
p =  v - T.*floor((v+0.5*T)./T);
end

