function folder = getManImResPath()
% getManImResPath()
%    returns the base path of the ManImRes Toolbox and returns an error
%    message, if the Toolbox is not initialized yet.
%
% ---
% ManImRes 1.0, R. Bergmann ~ 2014-11-29

folder = fileparts(which('initManImRes.m'));
assert(~isempty(folder),...
    'ManImRes Toolbox not found in path, please run initManImRes.m first.');
end

