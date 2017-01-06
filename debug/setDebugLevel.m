function setDebugLevel(varargin)
% setDebugLevel(type,level)
%   Set the global debug level of type to the value level, where
%   the type can be ommited, which yields setting types 'LevelMax' and 'LevelMin'
%      'LevelMax' is the upper bound, 'LevelMin' the lower bound for all types.
%   
%   INPUT
%       type  : (optional) a string of the type of debug the leve is set
%       level : level the type is set to
%
% ---
% LMDLO ~ R. Bergmann ~ 2014-04-21, last edit: 2014-04-22

    global LMDLO_struct;
    % if debug was not active before, calling setDebugLevel activates it.
    if ~isDebugActive() 
        initializeDebug();
    end
    if (nargin==0)
        setDebugLevel('LevelMax',Inf);
        setDebugLevel('LevelMin',0);
        return;
    elseif (nargin == 1)
        setDebugLevel('LevelMax',varargin{1});
        setDebugLevel('LevelMin',varargin{1});
        return;
    elseif (nargin>2)
        warning(['Wrong number of setDebugLevel (',num2str(nargin),...
            '). Ingoring function call.']);
        return;
    end
    type = varargin{1};
    level = varargin{2};
    if ~any(strcmp(type,LMDLO_struct.types))
        %type not yet included -> add types
        newind = length(LMDLO_struct.types)+1;
        LMDLO_struct.types{newind} = type;
        LMDLO_struct.levels(newind) = level;
        if strcmp(type,'LevelMin') && level > getDebugLevel('LevelMax')
            warning(['Trying to set ''LevelMin'' (',num2str(level),') > ''LevelMax'' (',num2str(getDebugLevel('LevelMax')),') is not possible, using the value of ''LevelMax''']);
            LMDLO_struct.levels(newind) = getDebugLevel('LevelMax');
        end
        debug('text',3,'Text',['Debug Level of ',type,' set to ',num2str(LMDLO_struct.levels(newind))]);
    else
        %update type
        thisind = find(strcmp(type,LMDLO_struct.types),1);
        LMDLO_struct.levels(thisind) = level;
        LMDLO_struct.levels(thisind) = level;
        if strcmp(type,'LevelMin') && level > getDebugLevel('LevelMax')
            warning(['Trying to set ''LevelMin'' (',num2str(level),') > ''LevelMax'' (',num2str(getDebugLevel('LevelMax')),') is not possible, using the value of ''LevelMax''']);
            LMDLO_struct.levels(thisind) = getDebugLevel('LevelMax');
        end
        debug('text',3,'Text',['Debug Level of ',type,' set to ',num2str(LMDLO_struct.levels(thisind))]);
    end
end