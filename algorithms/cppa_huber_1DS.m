function x = cppa_huber_1DS(f, alpha,lambda,tau,omega,varargin)
%cppa_tv12_1D(f,alpha,beta,lambda,q,f)
%
%   Compute the cyclic proximal point algorithm for first and
%   second order differences, weighted by alpha and beta
%
%   INPUT
%      f       : a vector from [-pi,pi)^d of phase valued data
%      alpha   : weight for the first order difference term
%      lambda  : weight for the proximal mappings
%      tau     : parameter for the huber functional to determine parabola
%      omega   : parameter for the huber functional to determine parabola
%                width
%
%   OUTPUT:
%      x       : result of the cyclic proximal point algorithm
%
%   OPTIONAL PARAMETERS
%      'MaxIterations' : (400) Maximal number of Iterations
%      'Epsilon'    : (10^{-6}) Lower bound for the max change in one cycle
%           
%           one of the parameters MaxIterations/Epsilon can be dactivated
%           by specifying them as Inf but not both.
%
%    SEE ALSO
%
%    R. Bergmann, F. Laus, G. Steidl, A. Weinmann: Second Order Differences
%        of S1-valued Data and Applications in Variational Image Denoising,
%        Preprint, 2014
%    R. Bergmann, A. Weinmann: Inpainting of Cyclic Data using First and
%        Second Order Differences, Preprint, 2014
% ---
% Manifold-valued Image Restoration Toolbox 1.0
% R. Bergmann ~ 2014-08-27 | 2016-02-22
% see LICENSE.txt

ip = inputParser;
addParameter(ip, 'MaxIterations',400);
addParameter(ip, 'Epsilon',10^(-6));
parse(ip, varargin{:});
par = ip.Results;
%
% Validate optional stuff
if sum(size(par.MaxIterations))>2
    error('The maximal number of iterations has to be a number.');
end
if floor(par.MaxIterations)~=par.MaxIterations
    error('The maximal number of iterations has to be an integer.');
end
maxiter = par.MaxIterations;
if sum(size(par.Epsilon))>2
    error('The upper bound for cyclic step change has to be a number');
end
epsilon = par.Epsilon;
% Error checks
if ( isinf(maxiter) && (isinf(epsilon)) )
    error('You can not set both the maximal Iterations and Epsilon to Inf');
end
%
n = size(f,2);
if (n==1) && (length(size(f))==2) % zero-dimensional-manifold i.e. values given as column -> prime
    f = f';
    n = length(f);
end
xold = zeros(size(f));
x = f;
i = 0;
lambdait=lambda;
%
debug('time',3,'StartTimer','proxtiming');
M = S1;
itD = Inf;
while ( (max(itD)>epsilon) && (i<maxiter) )
    xold = x;
    % First step : Distances, eventually masked
    x = M.proxDist(x,f,lambdait);
    % Second steps: Huber
    for j=0:1
        nh=floor((n-j)/2);
        x(j+1:2*nh+j) = reshape(proxS1Huber(reshape(x(j+1:2*nh+j),2,[]),alpha*lambdait,tau,omega),1,[]);
    end    
    i = i + 1;
    lambdait = lambda/(i);
    itD = M.dist(x,xold);
    if mod(i,getDebugLevel('IterationStep'))==0
        debug('text',3,'text',...
            ['i=',num2str(i),' where lastdiff is ',num2str(max(itD(:))...
            ),' and lambdait=',num2str(lambdait)]);
    end
end
debug('time',3,'StopTimer','proxtiming');
debug('text',2,'Text',...
    [num2str(i),' iterations, last difference:',...
    num2str(max(itD(:)))]);
end
