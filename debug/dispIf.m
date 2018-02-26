function dispIf(condition,str)
%  disIf(condition, str) display string if condition is true
% just a small helper for inline functions
%
% ---
% Manifold-valued Image Restoration Toolbox 1.2
% Ronny Bergmann | 2018-02-12
if condition
    disp(str)
end
end

