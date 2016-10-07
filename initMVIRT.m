function initMVIRT(varargin)
% initManImRes()
%   Initialize the Manifold Valued Data Processing Toolbox by adding all
%   needed folders relative to this file, where the neccessary files
%   reside and compile cpp-source files (optional).
%
% OPTIONAL PARAMETERS
%   'Make'   : (false) Set this flag to compile all mex-files
%   'UseMex' : (true) Set this flag to false to disable the usage of
%   C++Algorithms and use the (fallback) Matlab implementations
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-11-29 | 2015-04-12

ip = inputParser;
addParameter(ip,'Make',false);
addParameter(ip,'UseMex',true);
parse(ip, varargin{:});
vars = ip.Results;

%% Setup folders
folder = fileparts(which(mfilename));
addpath(...
    fullfile(folder, 'algorithms'),...
    genpath(fullfile(folder, 'data')),...
    fullfile(folder, 'debug'),...
    fullfile(folder, 'helpers'),...
    fullfile(folder, 'helpers/TikZ'),...
    fullfile(folder, 'helpers/export'),...
    fullfile(folder, 'manifolds'),...
    genpath(fullfile(folder,'examples')),...
    genpath(fullfile(folder,'mex'))...
);
%% Compile?
if vars.Make
    cd mex
    % The Sphere functions
    disp('--- Compiling Sphere functions ---');
    mex -I../include/eigen -Imanifolds SnDist.cpp manifolds/manifoldSn.cpp
    mex -I../include/eigen -Imanifolds SnExp.cpp manifolds/manifoldSn.cpp
    mex -I../include/eigen -Imanifolds SnLog.cpp manifolds/manifoldSn.cpp
    mex -I../include/eigen -Imanifolds SnMean.cpp manifolds/manifoldSn.cpp
    mex -I../include/eigen -Imanifolds SnParallelTransport.cpp manifolds/manifoldSn.cpp
    mex -I../include/eigen -Imanifolds SnGrad_X_D2.cpp manifolds/manifoldSn.cpp
    % Hyperbolic space
    disp('--- Compiling Hyperbolic Space functions ---');
    mex -I../include/eigen -Imanifolds HnDist.cpp manifolds/manifoldHn.cpp
    mex -I../include/eigen -Imanifolds HnDot.cpp manifolds/manifoldHn.cpp
    mex -I../include/eigen -Imanifolds HnExp.cpp manifolds/manifoldHn.cpp
    mex -I../include/eigen -Imanifolds HnLog.cpp manifolds/manifoldHn.cpp
    mex -I../include/eigen -Imanifolds HnMean.cpp manifolds/manifoldHn.cpp
    % SPD
    disp('--- Compiling SymPosDef functions ---');
    mex -I../include/eigen -Imanifolds SPDDist.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDDot.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDExp.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDGrad_X_D2.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDLog.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDMean.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDParallelTransport.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDMean.cpp manifolds/manifoldSPD.cpp
    cd ..
end
%% Disable Mex usage?
if ~vars.UseMex
    % where to save this?
end
%% Init Debug
setDebugLevel('LevelMin',0);
setDebugLevel('LevelMax',50);
setDebugLevel('time',3);
setDebugLevel('text',3);
% Add all necessary paths
debug('text',1,'Text','Initializing the Manifold-valued Image Restoration Toolbox.');
cd(getManImResPath());
%
