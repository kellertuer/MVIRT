function Res = NL_Mean(M, F, patch_size, window_size, K, sigma_patch, sigma_weights,center_element)
% NL_Mean(M, F, patch_size, window_size, k, sigma_patch, sigma_weights)
% Calculates a denoised image Res by an non local mean, where the similar
% pixels have similiar neighborhoods
%
% INPUT:
% M             manifold
% F             given image ManDim x ImgDim
% patch_size    size of patches we compare, must be odd
% window_size   size of box around each pixel in which we search similar
%               patches, must be odd
% K             number of neighbors for each pixels.
% sigma_patch   standard deviation of guassian for patch comparison (center more important than outer pixels)
% sigma_weights standard deviation of guassian for computation of actual
%               weights
% OPTIONAL
% center_element ['m'] Char indicating which weight to use for the center
%                   patch before averaging, 'm' maximum over neighbors
%                                           'z' dist=0 i.e. w = 1
% OUTPUT:
% Res           Restored manifold valued image
%
% See also
% non_local_weights_fast
%
% MVRIT 1.2 J. Persch 05.07.2016
if nargin < 8
    center_element = 'm';
end
assert(mod(patch_size,2) == 1,'Patch_size must be odd!');
assert(mod(window_size,2) == 1,'Window_size must be odd!');
dimen = size(F);
imgDim = dimen(length(M.ItemSize)+1:end);
N = prod(imgDim);
W = non_local_weights(M, F, patch_size, window_size, K, sigma_patch, sigma_weights,center_element);
disp('Weights have been calculated.');
F = reshape(F,prod(M.ItemSize),[]);
weights = zeros(N,max(sum(W~=0)));
FF = zeros([size(F),size(weights,2)]);
% Extract the non zero weights
[neighbors,ind_i,w] = find(W);
% get the indeces of each pixel
[~,ind_i,~] = unique(ind_i);
% add the last element for the loop
ind_i = [ind_i;length(w)];
% Capture the weights and neighbors of each pixel
for i = 1:N
    weights(i,1:(ind_i(i+1)-ind_i(i))) = w(ind_i(i):(ind_i(i+1)-1));
    FF(:,i,1:(ind_i(i+1)-ind_i(i))) = F(:,neighbors(ind_i(i):(ind_i(i+1)-1)));
end
% Calculate the mean
Res = reshape(M.mean(reshape(FF,[M.ItemSize,N,size(weights,2)]),'Weights',weights),dimen);
toc
end
