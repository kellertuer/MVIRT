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
%   For futher details see
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%     Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
% ---
% Manifold-valued Image Restoration Toolbox 1.1
%  J. Persch  ~ 2017-01-06 | R. Bergmann 2016-01-07
% see LICENSE.txt

if nargin < 8
    center_element = 'm';
end
assert(mod(patch_size,2) == 1,'Patch_size must be odd!');
assert(mod(window_size,2) == 1,'Window_size must be odd!');
dimen = size(F);
imgDim = dimen(length(M.ItemSize)+1:end);
N = prod(imgDim);
debug('time',3,'StartTimer','NL_Mean');
W = non_local_weights(M, F, patch_size, window_size, K, sigma_patch, sigma_weights,center_element);
debug('text',3,'text',...
            'Weights have been calculated.');
F = reshape(F,prod(M.ItemSize),[]);
weights = zeros(N,max(sum(W~=0)));
FF = zeros([size(F),size(weights,2)]);
% Capture the weights and neighbors
for i = 1:N
    w = W(i,W(i,:)~=0);
    weights(i,1:length(w)) = w;
    FF(:,i,1:length(w)) = F(:,W(i,:)~=0);
end
% Calculate the mean
Res = reshape(M.mean(reshape(FF,[M.ItemSize,N,size(weights,2)]),'Weights',weights),dimen);
debug('time',3,'StopTimer','NL_Mean');
end
