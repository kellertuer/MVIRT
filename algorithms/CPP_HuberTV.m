function [x,recData] = CPP_HuberTV(varargin)
%   CPP_AddiviteTV12(M,f,alpha,tau,omega)
% Compute the CPP of the first and second order TV with additive coupling
% can also be called with a struct with the corresponding fields, i.e. as
%   CPP_AdditiveTV12(problem).
% 
% INPUT
%   M                 ; a manifold.
%   f                 : an image from M^{m_1,m_2,...,m_n} of manifold-valued data
%   alpha             : weight(s) for first order absolute difference terms
%                       either a number or one for each dimension.
%   tau, omega        : weight(s) for each weight for second order absolute
%                       difference terms either two nubers (straight &
%                       mixed differences) or one for each (n + n-1+...+1)
%   
%   stoppingCriterion : a functional @(x,xold,lambda,iter) as a stopping
%                       criterion based on x, the last iterate xold, lambda
%
%   OUTPUT:
%      x   : result of the cyclic proximal point algorithm
%  recData : recorded data array of the by the record-function
%
%   OPTIONAL PARAMETERS
%
%   lambdaIterate     : (@ (iter) lambda/iter) a functional @(iter) computing the iterates value of lambda
%   Debug             : ([]) a function @ (x,xold,iter) producing debug
%                          output
%   Record            : ([]) a function @(x,iter)) a function returning a column vector,
%                       if there's a second return value and this function
%                       is given, data is recorded in an array and returned
%
% 'UnknownMask'       : ([]) Specify a mask the same size as the data, i.e.,
%                        m_1xm_2x...xm_nw, which are labeled as unknown
%                        they are initialized in the cycles, when possible.
%
% 'FixedMask'         : ([]) Specify a binary mask (m_y1xm_2x...xm_n) for
%                       values that are fixed. They are fixed on a 'poor
%                       mans method', i.e. reset after each step
%
% R. Bergmann, M. Bacak, G. Steidl, A. Weinmann,
%     A second order non-smooth variational model for restoring
%     manifold-valued images,
%     SIAM Journal on Scientific Computing, 38, (1), A567?A597, 2016.
% ---
% MVIRT | R. Bergmann | 2018-01-22
if (length(varargin)==1 && isstruct(varargin{1})) % as struct
    % short hand call to surpass parser
    vars = varargin{1};
    assert(isfield(vars,{'M'}),'A manifold is not given');
    assert(isa(vars.M,'manifold'),'Input M is not a manifold');
    assert(isfield(vars,{'f'}),'no data x given');
    assert(isfield(vars,{'alpha'}),'a weight fot the TV term is missing.');
    assert(isfield(vars,{'tau'}),'the Huber parameter tau is missing.');
    assert(isfield(vars,{'omega'}),'the Huber parameter omega is missing.');
    assert(isfield(vars,{'stoppingCriterion'}),'a stopping criterion function is missing');
    if ~isfield(vars,{'Debug'})
        vars.Debug = [];
    end
    if ~isfield(vars,{'lambdaIterate'})
        vars.lambdaIterate = @(iter) pi/iter;
    end
    if ~isfield(vars,{'Record'})
        vars.Record = [];
    end
    if ~isfield(vars,{'UnknownMask'})
        vars.UnknownMask = [];
    end
    if ~isfield(vars,{'FixedMask'})
        vars.FixedMask = [];
    end
else % parse Input
    ip = inputParser;
    addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
    addRequired(ip,'f');
    addRequired(ip,'proximalMaps');
    addRequired(ip,'alpha');
    addRequired(ip,'tau');
    addRequired(ip,'omega');
    addRequired(ip,'stoppingCriterion');
    addOptional(ip,'Debug',[]);
    addOptional(ip,'lambdaIterate',@(iter) pi/iter);
    addOptional(ip,'Record',[]);
    addOptional(ip,'UnknownMask',[]);
    addOptional(ip,'FixedMask',[]);
    parse(ip, varargin{:});
    vars = ip.Results;
end
f = vars.f;
f(vars.M.allDims{:},vars.UnknownMask) = NaN;
proximalMaps = {...
    @(x,lambda) proxDistanceSquared(vars.M,x,f,lambda,vars.FixedMask),...
    @(x,lambda) proxTV(vars.M,x,lambda.*vars.alpha,...%vars.FixedMask,...
        'DifferenceProx',@(x1,x2,lambda) ...
            proxAbsoluteDifferenceHuber(vars.M,x1,x2,lambda,vars.tau,vars.omega))};
if isempty(vars.Record) || nargout < 2 % no record
    x = cyclicProximalPoint(vars.M,f,proximalMaps,vars.stoppingCriterion,...
            'lambdaIterate',vars.lambdaIterate,...
            'Debug',vars.Debug);
else
    [x,recData] = cyclicProximalPoint(vars.M,f,proximalMaps,vars.stoppingCriterion,...
            'lambdaIterate',vars.lambdaIterate,...
            'Debug',vars.Debug,...
            'Record',vars.Record);
end
end

