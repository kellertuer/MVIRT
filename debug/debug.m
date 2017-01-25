function debug(type, level, varargin)
% Debug(type,level)
%
% Provide Debug output for a certain Level of Debug (see @setDebugLevel)
% for different levels of debug available. Further information concerning
% the debug is given by optional parameters
% 
%   INPUT:
%        type  : type, the debug belongs to
%        level : number of the level, the debug belongs to
%
%   OPTIONAL
%        'StartTimer' : start global timer saved in the global struct timer
%        'StopTimer'  : stop global timer number
%        'GetTimer'   : get global timer
%        'TimeFormat' : number format for the timer (see fprintf, std: (%09.8f) 
%        'Text'       : Display a Text
%
%   OUTPUT:
%      None, despite Text output during the method
% ---
% LMDLO ~ R. Bergmann ~ 2014-04-21

    global LMDLO_active;
    if isempty(LMDLO_active) %Shortcut - shirter that isActive to spare time.
        return;
    end
    persistent timers;
    if ~isa(timers,'containers.Map') % not initialized
        timers = containers.Map();
        debug('text',3,'Text','Debug Timers initialized.');
    end
    % check input level struct
    if ~isnumeric(level)
        warning('Debug level not a number, ignoring this debug action');
        return;
    end
    % Parse parameters
    p = inputParser;
    addParameter(p, 'StartTimer','');
    addParameter(p, 'StopTimer','');
    addParameter(p, 'GetTimer','');
    addParameter(p, 'TimeFormat','%09.8f');
    addParameter(p, 'Text','');
    parse(p, varargin{:});
    debugparams = p.Results;
    %
    %
    %  In the following the first argument checks, whether this debug call
    %  should do something for a special case, while the second to last
    %  checks determine, whether the debug level contains the case from the
    %  first comparison
    %
    % Text outputs
    
    if ( (getDebugLevel('text')>=level) && (strcmp(type,'text')))
        disp(debugparams.Text);
    end
    % Timers
    if ( (getDebugLevel('time')>=level) && (strcmp(type,'time')) )
        if (numel(debugparams.StartTimer)>0)
            if (timers.isKey(debugparams.StartTimer))
                warning(['Key ''',debugparams.StartTimer,''' already in use for a running time measurement; restarting timer.']);
                timers.remove(debugparams.StartTimer);
            end
                timers = [timers; containers.Map({debugparams.StartTimer},{tic})];
                debug('text',3,'Text',['Timer started for ',debugparams.StartTimer]);
        end
        if (numel(debugparams.StopTimer)>0)
            if (timers.isKey(debugparams.StopTimer))
                temptime = toc(timers(debugparams.StopTimer));
                debug('text',3,'Text',['Timer stopped for ',debugparams.StopTimer,'.']);
                debug('text',2,'Text',['Computation for ',debugparams.StopTimer,' took ',sprintf(debugparams.TimeFormat,temptime),' sec.']);
                timers.remove(debugparams.StopTimer);
            else
                warning(['There was no timer running for ',debugparams.StopTimer]);
            end
        end
        if (numel(debugparams.GetTimer)>0)
            if (timers.isKey(debugparams.GetTimer))
                temptime = toc(timers(debugparams.GetTimer));
                debug('text',2,'Text',['Computation for ',debugparams.GetTimer,' runs already ',sprintf(debugparams.TimeFormat,temptime),' sec.']);
            else
                warning(['There was no timer running for ',debugparams.GetTimer]);
            end
        end
    end
end