function act = isDebugActive()
% isDebugActive()
%    returns whether the debug system is activated or not
%
% ---
% LMDLO ~ R. Bergmann ? 2014-11-29

global LMDLO_active LMDLO_struct;
act = ~ (isempty(LMDLO_active)...
    || (~isfield(LMDLO_struct,'types'))...
    || (~isfield(LMDLO_struct,'levels'))  );
end

