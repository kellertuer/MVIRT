function lvl = getDebugLevel( type )
%   getDebugLevel(type)
%       get Debug Level of a certain type, 0 if not set
%   INPUT
%       type    : type a of debug
%
%   OUTPUT
%       lvl     : level in {1,2,3} of that debug type, if set,
%                   this level is bounded by `all` (if set) from above
%                   this lebel is bounded by `any` (if set) from below
%                 0 is returned, if none of these three is set
% ---
% LMDLO, R. Bergmann, 2014-04-21
    global LMDLO_active LMDLO_struct;
    if (isempty(LMDLO_active)) %not initailized by any setDebugLevel
        warning('getDebugLevel() should not be called without initialization first');
        lvl = -1;
        return;
    end
    lvl=0;
    if any(strcmp(type,LMDLO_struct.types))
        %For this type exists an explicit value
        lvl = LMDLO_struct.levels(find(strcmp(type,LMDLO_struct.types),1));
    end
    lvlmax=Inf; lvlmin=0;
    if any(strcmp('LevelMax',LMDLO_struct.types)) %is all set? remember that level
        lvlmax = LMDLO_struct.levels(find(strcmp('LevelMax',LMDLO_struct.types),1));
    end
    if any(strcmp('LevelMin',LMDLO_struct.types)) %is any set? remember that level
        lvlmin = LMDLO_struct.levels(find(strcmp('LevelMin',LMDLO_struct.types),1));
    end
    lvl = min(max(lvlmin,lvl),lvlmax);
end