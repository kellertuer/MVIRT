function folder = getManImResPath()
% getManImResPath()
%    returns the base path of the MVIRT Toolbox and returns an error
%    message, if the Toolbox is not initialized yet.
%
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-11-29
% see LICENSE.txt

folder = fileparts(which('initMVIRT.m'));
assert(~isempty(folder),...
    'MVIRT not found in path, please run initMVIRT.m first.');
end

