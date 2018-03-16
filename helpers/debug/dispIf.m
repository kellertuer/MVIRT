function dispIf(condition,str)
% dispIf(condition,str) display string if condition is true.
% This small helper function is necessary for inline functions since there
% is no inline if.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.2
% R. Bergmann | 2018-02-12
if condition
    disp(str)
end
end
