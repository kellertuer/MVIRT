function [u_selected,u_selected_ref,center_x, center_y] = get_patches(M,u,u_oracle,pos_x,pos_y,patch_size, window_size,K,L1)
% get_patches(M,u,u_oracle,pos_x,pos_y,patch_size, window_size,K,L1)
% Finds K similar patches of size patch_size to the patch centered at
% pos_x,pos_y. Similarity is measured in u_oracle and only patches inside a
% window_size window are considered.

% INPUT:
% M                                 manifold on which the patch lives
% u                                 noisy image
% u_oracle                          oracle image
% pos_x,pos_y                       position of the center of the patch
% patch_size x patch_size           patch size, must be odd
% window_size x window_size         size of window around each pixel in which we search similar patches, must be odd
% K                                 number of similar patches
% L1                                States if L1-distance is used to select similar
%                                   patches
% OUTPUT:
% u_selected                        selected patches in noisy image
% u_selected_ref                    selected patches in the oracle image
% [center_x, center_y]                center of selected patches
%
%   For further details see
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%     Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
% ---
% Manifold-valued Image Restoration Toolbox 1.2
% J. Persch  ~2017-07-05 | R. Bergmann 2016-01-07
% see LICENSE.txt

dimen = size(u);
m = dimen(length(M.ItemSize)+1);
n =  dimen(length(M.ItemSize)+2);
d = prod(M.ItemSize);
u = reshape(u,[],m,n);
u_oracle = reshape(u_oracle,[],m,n);

kappa = (patch_size - 1)/2;
nu = (window_size - 1)/2;
delta = nu-kappa;

% reference patch
u_ref = u_oracle(:,(pos_x-kappa):(pos_x+kappa),(pos_y-kappa):(pos_y+kappa));

% indices of the patches
x = max(pos_x-delta,kappa+1):min(pos_x+delta,m-kappa);
y = max(pos_y-delta,kappa+1):min(pos_y+delta,n-kappa);

X = kron(x,ones(1,patch_size))+kron(ones(1,length(x)),-kappa:kappa);
Y = kron(y,ones(1,patch_size))+kron(ones(1,length(y)),-kappa:kappa);

% indices of centers
[C_x,C_y] = meshgrid(x,y);
C_x = C_x';
C_y = C_y';

centers_x = C_x(:);
centers_y = C_y(:);

% all patches of the neighborhood
u_patches = u(:,X,Y);
u_patches_ref = u_oracle(:,X,Y,:);

% reorder patches
data = permute(reshape(permute(reshape(u_patches,[d,length(x)*patch_size,patch_size,length(y)]),[1,3,2,4]),[d,patch_size,patch_size,length(x)*length(y)]),[1,3,2,4]);
data_ref = permute(reshape(permute(reshape(u_patches_ref,[d,length(x)*patch_size,patch_size,length(y)]),[1,3,2,4]),[d,patch_size,patch_size,length(x)*length(y)]),[1,3,2,4]);


% distances = sum(sum((data_ref - repmat(u_ref,[1,1,length(x)*length(y)])).^2,2),1);

if L1 
    distances = sum(sum(abs(M.dist(reshape(data_ref,[M.ItemSize,patch_size,patch_size,length(x)*length(y)]),...
                                  reshape(repmat(u_ref,[1,1,1,length(x)*length(y)]),[M.ItemSize,patch_size,patch_size,length(x)*length(y)]))),2),1);
else
    
distances = sum(sum(M.dist(reshape(data_ref,[M.ItemSize,patch_size,patch_size,length(x)*length(y)]),...
                                  reshape(repmat(u_ref,[1,1,1,length(x)*length(y)]),[M.ItemSize,patch_size,patch_size,length(x)*length(y)])).^2,2),1);
end

 [~,idx] = sort(distances(:));

u_selected = reshape(data(:,:,:,idx(1:K)),[M.ItemSize,patch_size,patch_size,K]);
u_selected_ref = reshape(data_ref(:,:,:,idx(1:K)),[M.ItemSize,patch_size,patch_size,K]);

center_x = centers_x(idx(1:K));
center_y = centers_y(idx(1:K));
end

