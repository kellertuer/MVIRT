function v = l2TVFunctional( varargin )
% TVFunctionalValue(M,f,alpha,beta) - compute the value of the TV-type
%   functional of first and second order for manifold valued data
%
% Input
%   M     : a manifold
%   f     : data (1D or 2D)
%   alpha : weight(s) for the TV term in the functional
%   beta  : weight(s) for the second order TV term in the functional 
%
% OUTPUT
%   v     : the value of the TV functional
%
% OPTIONAL PARAMETERS
%   'DomainMask'  : (true) only includes terms from the domain mask into the
%       functional in all terms
%   'UnknownMask' : (true) only includes the TV terms ot the indicated terms into
%       the functional but not the data terms.
% ---
% ManImRes, R. Bergmann ~ 2015-11-24

ip = inputParser;
addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}))
addRequired(ip,'f');
addRequired(ip,'x');
addRequired(ip,'alpha');
addRequired(ip,'beta');
addParameter(ip,'DomainMask',NaN);
addParameter(ip,'UnknownMask',NaN);
parse(ip, varargin{:});
vars = ip.Results;
dimen = size(vars.f);
manDims = length(vars.M.ItemSize);
imgDim = dimen((manDims+1):end);
manDim = dimen(1:manDims);
M = vars.M;
if isnan(vars.DomainMask)
    if isscalar(imgDim)
        dM = true(imgDim,1);
    else
        dM = true(imgDim);
    end
elseif all(size(vars.DomainMask)~=imgDim)
    error('Mask for Domain does not match image dimensions');
else
    dM = vars.DomainMask;
end
if isnan(vars.UnknownMask)
    if isscalar(imgDim)
        uM = false(imgDim,1);
    else
    uM = false(imgDim);
    end
elseif all(size(vars.UnknownMask)~=imgDim)
    error('Mask for Domain does not match image dimensions');
else
    uM = vars.UnknownMask;
end
dM = ~dM;
dists = M.dist(vars.f,vars.x);
dists(dM|uM) = 0;
v =1/2* sum(dists(:).^2);

%
% Signals
if length(imgDim)==1 % signal
    xm = reshape(vars.x,[prod(manDim),imgDim]);
    if vars.alpha>0
        firstOrder =  M.dist(reshape(xm(:,1:end-1),[manDim,imgDim-1]),...
            reshape(xm(:,2:end),[manDim,imgDim-1]));
        firstOrder(dM(1:end-1)|dM(2:end))=0;
        v = v + sum(vars.alpha*firstOrder(:));
    end
    if vars.beta>0
        secondOrder =  M.dist(reshape(xm(:,2:end-1),[manDim,imgDim-2]),...
            M.midPoint(reshape(xm(:,1:end-2),[manDim,imgDim-2]),reshape(xm(:,3:end),[manDim,imgDim-2])));
        secondOrder(dM(2:end-1)|dM(1:end-2)|dM(3:end)) = 0;
        v = v + sum(vars.beta*secondOrder(:));
    end
end
%
% Images
if length(imgDim)==2
    if (length(vars.alpha)==1) % use Axis TV
        alpha = [1,1,0,0]*vars.alpha;
    elseif (length(vars.alpha)<5) % set the first ones
        alpha = zeros(1,4);
        alpha(1:length(vars.alpha)) = vars.alpha;
    else
        error('the vector alpha is too long');
    end
    %
    % beta (Check and form to a vector for future simplicity)
    if (length(vars.beta)==1) % number
        beta = [1,1,1]*vars.beta;
    elseif (length(vars.beta)==2) %
        beta = [vars.beta(1),vars.beta(1),vars.beta(2)];
    elseif (length(vars.beta)==3)
        beta = vars.beta;
    else
        error('beta is a vector longer than 3 and hence too long.');
    end
    xm = reshape(vars.x,[prod(manDim),imgDim]);
    if alpha(1)>0
        firstOrder =  M.dist(reshape(xm(:,1:end-1,:),[manDim,imgDim-[1,0]]),...
            reshape(xm(:,2:end,:),[manDim,imgDim-[1,0]]));
        firstOrder(dM(1:end-1,:)|dM(2:end,:)) =0;
        v = v + alpha(1)*sum(firstOrder(:));
    end
    if alpha(2)>0
        firstOrder =  M.dist(reshape(xm(:,:,1:end-1),[manDim,imgDim-[0,1]]),...
            reshape(xm(:,:,2:end),[manDim,imgDim-[0,1]]));
        firstOrder(dM(:,1:end-1)|dM(:,2:end)) =0;
        v = v + alpha(2)*sum(firstOrder(:));
    end
    if alpha(3)>0
        firstOrder =  M.dist(reshape(xm(:,1:end-1,1:end-1),[manDim,imgDim-[1,1]]),...
            reshape(xm(:,2:end,2:end),[manDim,imgDim-[1,1]]));
        firstOrder(dM(1:end-1,1:end-1)|dM(2:end,2:end)) =0;
        v = v + alpha(3)*sum(firstOrder(:));
    end
    if alpha(4)>0
        firstOrder =  M.dist(reshape(xm(:,1:end-1,2:end),[manDim,imgDim-[1,1]]),...
            reshape(xm(:,2:end,1:end-1),[manDim,imgDim-[1,1]]));
        firstOrder(dM(1:end-1,2:end)|dM(2:end,1:end-1)) =0;
        v = v + alpha(4)*sum(firstOrder(:));
    end
    if beta(1)>0
        secondOrder =  M.dist(reshape(xm(:,2:end-1,:),[manDim,imgDim-[2,0]]),...
            M.midPoint(reshape(xm(:,1:end-2,:),[manDim,imgDim-[2,0]]),...
            reshape(xm(:,3:end,:),[manDim,imgDim-[2,0]])));
        secondOrder(dM(2:end-1,:)|dM(1:end-2,:)|dM(3:end,:)) = 0;
        v = v + beta(1)*sum(secondOrder(:));
    end
    if beta(2)>0
        secondOrder =  M.dist(reshape(xm(:,:,2:end-1),[manDim,imgDim-[0,2]]),...
            M.midPoint(reshape(xm(:,:,1:end-2),[manDim,imgDim-[0,2]]),...
            reshape(xm(:,:,3:end),[manDim,imgDim-[0,2]])));
        secondOrder(dM(:,2:end-1)|dM(:,1:end-2)|dM(:,3:end)) = 0;
        v = v + beta(2)*sum(secondOrder(:));
    end
    if beta(3)>0
        secondOrder =  M.dist(M.midPoint(reshape(xm(:,1:end-1,1:end-1),[manDim,imgDim-[1,1]]),...
            reshape(xm(:,2:end,1:end-1),[manDim,imgDim-[1,1]])),...
            M.midPoint(reshape(xm(:,1:end-1,2:end),[manDim,imgDim-[1,1]]),...
            reshape(xm(:,2:end,2:end),[manDim,imgDim-[1,1]])));
        secondOrder(dM(1:end-1,1:end-1)|dM(2:end,1:end-1)|dM(1:end-1,1:end-1)|dM(2:end,2:end)) = 0;
        v = v + beta(3)*sum(secondOrder(:));
    end
end
end