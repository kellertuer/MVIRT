function deactivateDebug()
% deactivateDebug()
%   clear all local variables and hence deactivate Debug
%
% ---
% LMDLO ~ R. Bergmann ~ 2014-04-21

    debug('text'3,'Text','Clearing LDMDLO global variables, deactivating debug completely.');
    clearvars -global LDMDLO_active LMDLO_struct;
end