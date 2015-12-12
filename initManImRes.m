function initManImRes(varargin)
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
% Manifold Valued Image Restoration 1.0
% R. Bergmann ~ 2014-11-29 | 2015-04-09

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
    mex -I../include/eigen -Imanifolds SnGradX.cpp manifolds/manifoldSn.cpp
    % SPD
    disp('--- Compiling SymPosDef functions ---');
    mex -I../include/eigen -Imanifolds SPDDist.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDDot.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDExp.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDGradX.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDLog.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDParallelTransport.cpp manifolds/manifoldSPD.cpp
    mex -I../include/eigen -Imanifolds SPDMean.cpp manifolds/manifoldSPD.cpp
    cd ..
end% Init Debug
setDebugLevel(2);
% Add all necessary paths
debug('text',1,'Text','Initializing the MANifold Valued IMage REStoration Toolbox.');
cd(getManImResPath());
%
