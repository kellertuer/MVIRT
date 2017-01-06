function [u_final,u_oracle] = NL_MMSE_2D(varargin)
% NL_MMSE_2D(M,f,sigma)
% Denoises a manifold valued image with the non-local mmmse estimator
% INPUT:
% M:                                manifold of the image
% f:                                given image
% sigma:                            standard deviation of the noise
% OPTIONAL:
% _1,_2 for first and second step only _1 uses the parameters for 2 steps
% patch_size x patch_size:          patch size, must be odd
% window_size x window_size:        size of window around each pixel in which we search similar patches, must be odd
% K:                                number of similar patches that are kept
% oracle:                           states if an orcalce image should be
%                                   used
% processed                         [true] switches the accelerartion
% hom_test:                         [true] indicates if homogeneous area
% L_1:                              changes the similarity measure to L1
% simple:                           [False] states if only means should be
%                                   considered
% OUTPUT:
% u_final:                          denoised image
% u_oracle:                         oracle image
%
%   For further details see
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%     Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
% ---
% Manifold-valued Image Restoration Toolbox 1.2
%  J. Persch  ~2017-01-06 | R. Bergmann 2016-01-07
% see LICENSE.txt
if (length(varargin)==1 && isstruct(varargin{1}))
    % short hand call to surpass parser
    vars = varargin{1};
else
    ip = inputParser;
    
    addRequired(ip,'M', @(x) validateattributes(x,{'manifold'},{}));
    addRequired(ip,'f');
    addRequired(ip,'sigma');
    
    addParameter(ip,'patch_size_1', 3);
    addParameter(ip,'patch_size_2', []);
    addParameter(ip,'window_size_1',[]);
    addParameter(ip,'window_size_2',[]);
    addParameter(ip,'K_1',[]);
    addParameter(ip,'K_2',[]);
    addParameter(ip,'gamma',1.2);
    addParameter(ip,'oracle',true);
    addParameter(ip,'L1',false);
    addParameter(ip,'processed',true);    
    addParameter(ip,'hom_test',true);
    addParameter(ip,'simple',false);
    parse(ip, varargin{:});
    vars = ip.Results;
end

assert(all(isfield(vars,{'M','f','sigma'})),'There are required input fields missing'); % at least required
M = vars.M;
u = vars.f;
sigma = vars.sigma;
% parameter initialization
dimen = size(u);
% N = m*n;
m = dimen(end-1);
n = dimen(end);
manDim = M.ItemSize;
if isfield(vars,'processed')
    processed_switch = vars.processed;
else
    processed_switch = true;
end
if isfield(vars,'patch_size_1')
    patch_size_1 = vars.patch_size_1;
else
    patch_size_1=3;
end
if mod(patch_size_1,2)==0
    patch_size_1 = patch_size_1+1;
    warning('Even patch_size_1 got incremented by 1');
end
kappa_1 = (patch_size_1 - 1)/2;
if ~isfield(vars,'patch_size_2') || isempty(vars.patch_size_2)
    patch_size_2 = patch_size_1;
else
    patch_size_2 = vars.patch_size_2;
    if mod(patch_size_2,2)==0
        patch_size_2 = patch_size_2+1;
        warning('Even patch_size_2 got incremented by 1');
    end
end
kappa_2 = (patch_size_2 - 1)/2;

if ~isfield(vars,'window_size_1') ||isempty(vars.window_size_1)
    window_size_1 = ceil(sqrt(M.Dimension))*patch_size_1*4+1;
else
    window_size_1 = vars.window_size_1;
    if mod(window_size_1,2)==0
        window_size_1 = window_size_1+1;
        warning('Even window_size_1 got incremented by 1');
    end
end
if ~isfield(vars,'window_size_2') || isempty(vars.window_size_2)
    window_size_2 = window_size_1;
else
    window_size_2 = vars.window_size_2;
    if mod(window_size_2,2)==0
        window_size_2 = window_size_2+1;
        warning('Even window_size_2 got incremented by 1');
    end
end

if ~isfield(vars,'K_1') || isempty(vars.K_1)
    K_1 = M.Dimension*patch_size_1.^2*3;
else
    K_1 = vars.K_1;
end
if ~isfield(vars,'K_2') ||  isempty(vars.K_2)
    K_2 = K_1;
else
    K_2 = vars.K_2;
end
if K_1 > ((window_size_1-1)/2-(patch_size_1-1)/2+1)^2
    K_1 = ((window_size_1-1)/2-(patch_size_1-1)/2+1)^2;
    warning(['K_1 too large, set it to',num2str(K_1)]);
end
if K_2 > ((window_size_2-1)/2-(patch_size_2-1)/2+1)^2
    K_2 = ((window_size_2-1)/2-(patch_size_2-1)/2+1)^2;
    warning(['K_2 too large, set it to',num2str(K_2)]);
end
if isfield(vars,'gamma')
    gamma = vars.gamma;
else
    gamma = 1.2;
end
if isfield(vars,'hom_test')
    hom_test = vars.hom_test;
else
    hom_test = true;
end
if isfield(vars,'L1')
    L1 = vars.L1;
else
    L1 = false;
end
if isfield(vars,'oracle')
    oracle = vars.oracle;
else
    oracle = true;
end
if isfield(vars,'simple')
    simple = vars.simple;
else
    simple = false;
end
%Initialization
u_oracle = cell(m,n);
processed = zeros(m,n);
processed2 = zeros(m,n);
runs = 0;
number_mean = 0;
step = 1;
debug('time',3,'StartTimer','NL_MMSE');
for pos_x = ( kappa_1+1) : m-(kappa_1)
    for pos_y = (kappa_1+1) : n-(kappa_1)
        if processed(pos_x,pos_y)==0
            runs = runs +1;
            % select K most similar patches from the neighborhood
            [u_selected,u_selected_ref,center_x, center_y] = get_patches(M,u,u,pos_x,pos_y,patch_size_1, window_size_1,K_1,L1);           
            % compute MMSE estimate
            [u_estimate,mu,S] = compute_MMSE(M,u_selected,u_selected_ref,sigma,step,simple);
            % homogeneous area test (After MMSE because we want to check condition of S as well)
            if hom_test && (M.var(u_selected)< gamma*sum(sigma)^2 || abs(cond(S))>10^6)
                u_estimate = repmat(M.mean(reshape(mu,[M.ItemSize,1,patch_size_1^2])),[ones(1,length(M.ItemSize)),patch_size_1,patch_size_1,K_1]);
                number_mean = number_mean+1;            
            end           
            u_estimate = reshape(mat2cell(reshape(u_estimate,prod(manDim),[]),prod(manDim),ones(1,patch_size_1^2*K_1)),[patch_size_1,patch_size_1,K_1]);
            % aggregation
            for i = 1:K_1
                u_oracle((center_x(i)-kappa_1):(center_x(i)+kappa_1),(center_y(i)-kappa_1):(center_y(i)+kappa_1)) = ...
                    cellfun(@(x,y) {[x,y]},u_oracle((center_x(i)-kappa_1):(center_x(i)+kappa_1),(center_y(i)-kappa_1):(center_y(i)+kappa_1)),...
                    u_estimate(:,:,i));
                if processed_switch
                    processed(center_x(i),center_y(i)) = 1;
                end
            end
        end
    end
    
end
if number_mean/runs  > 0.5
    warning('More than 50% positive homogeneous area tests in the first step eventually choose smaller gamma');
end
runs = 0;
number_mean = 0;
u_oracle = cell2mat(cellfun(@(x) {M.mean(reshape(x,[manDim,1,size(x,2)]))},u_oracle));
% Reshapes work only for vector or matrix manifolds
if isscalar(manDim)
    u_oracle = reshape(u_oracle,[manDim,m,n]);
else
    u_oracle = reshape(permute(reshape(permute(reshape(u_oracle,manDim(1)*m,manDim(2),[]),[2,1,3]),manDim(2),manDim(1),[]),[2,1,3]),manDim(1),manDim(2),m,n);
end
u_final = u_oracle;
debug('time',3,'GetTimer','NL_MMSE');
% eventually, perform a second step using the first step denoise image
if oracle
    debug('text',3,'text',...
            'Oracle image is calculated starting final step.');
    step = 2;
    u_final = cell(m,n);   
    
    for pos_x = ( kappa_2+1) : m-(kappa_2)
        for pos_y = (kappa_2+1) : n-(kappa_2)
            if processed2(pos_x,pos_y)==0
                runs = runs +1;
                % select K most similar patches from the neighborhood
                [u_selected,u_selected_ref,center_x, center_y] = get_patches(M,u,u_oracle,pos_x,pos_y,patch_size_2, window_size_2,K_2,L1);
                
                
                % compute Bayes estimate
                [u_estimate,mu,S] = compute_MMSE(M,u_selected,u_selected_ref,sigma,step,simple);
                %u_estimate = repmat(mu,1,1,1,K_2);
                % homogeneous area test
                if hom_test && (M.var(u_selected)< gamma*sigma^2|| abs(cond(S))>10^6)
                    u_estimate = repmat(M.mean(reshape(mu,[M.ItemSize,1,patch_size_2^2])),[ones(1,length(M.ItemSize)),patch_size_2,patch_size_2,K_2]);
                    number_mean = number_mean+1;
                end
               u_estimate = reshape(mat2cell(reshape(u_estimate,prod(manDim),[]),prod(manDim),ones(1,patch_size_2^2*K_2)),[patch_size_2,patch_size_2,K_2]);
                % aggregation
                for i = 1:K_2
                    u_final((center_x(i)-kappa_2):(center_x(i)+kappa_2),(center_y(i)-kappa_2):(center_y(i)+kappa_2)) = ...
                        cellfun(@(x,y) {[x,y]},u_final((center_x(i)-kappa_2):(center_x(i)+kappa_2),(center_y(i)-kappa_2):(center_y(i)+kappa_2)),...
                        u_estimate(:,:,i));
                    if processed_switch
                        processed2(center_x(i),center_y(i)) = 1;
                    end
                end
            end
        end        
    end    
    u_final = cell2mat(cellfun(@(x) {M.mean(reshape(x,[manDim,1,size(x,2)]))},u_final));
    % Reshapes work only for vector or matrix manifolds
    if isscalar(manDim)
        u_final = reshape(u_final,[manDim,m,n]);
    else
        u_final = reshape(permute(reshape(permute(reshape(u_final,manDim(1)*m,manDim(2),[]),[2,1,3]),manDim(2),manDim(1),[]),[2,1,3]),manDim(1),manDim(2),m,n);
    end
end
if number_mean/runs  > 0.5
    warning('More than 50% positive homogeneous area tests in the first step eventually choose smaller gamma');
end
debug('time',3,'StopTimer','NL_MMSE');
end




