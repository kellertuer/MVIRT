function W = non_local_weights(M, F, patch_size, window_size, k, sigma_patch, sigma_weights,center_element)
% non_local_weights(M, F, patch_size, window_size, K, sigma_patch, sigma_weights)
% computes a sparse matrix W where W_{i,j} specifies the similarity between
% pixel i and j of a vectorized version of the image F, where at least K
% similiar pixels are chosen, maybe more due to symmetries.
%
% INPUT:
% M                 Manifold
% F                 given image ManDim x ImgDim
% patch_size  (odd) specifies the neighborhood in which the pixels should
%                   be similar
% window_size (odd) Size of search window        
% K                 minimal number of neighbors for each pixels 
% sigma_patch       standard deviation of guassian for patch comparison
% sigma_weights     standard deviation of guassian for actual weights
% OPTIONAL
% center_element ['m'] Char indicating which weight to use for the center
%                   patch before averaging, 'm' maximum over neighbors
%                                           'z' dist=0 i.e. w = 1
%
%   For further details see
%     F. Laus, M. Nikolova, J. Persch, G. Steidl: A Nonlocal Denoising
%     Algorithm for Manifold-Valued Images Using Second Order Statistics.
%     (ArXiv Preprint 1607.08481)
% ---
% Manifold-valued Image Restoration Toolbox 1.0
%  J. Persch  ~2017-07-05 | R. Bergmann 2016-01-07
% see LICENSE.txt
if nargin < 8
    center_element = 'm';
end    
dimen = size(F);
imgDim = dimen(length(M.ItemSize)+1:end);
N = prod(imgDim);
m = imgDim(1);
n = imgDim(2);
p = (patch_size - 1)/2;
w = (window_size - 1)/2;
F = reshape(F,[prod(M.ItemSize),imgDim]);
F_mirror = F(:,[1+p:-1:2,1:end,(end-1):-1:(end-p)],...
                    [1+p:-1:2,1:end,(end-1):-1:(end-p)]);
ExtImgDim = imgDim + 2*p*ones(size(imgDim));
%Gaussian for computing the patch differences
g = exp(-(-p:p).^2/(2*sigma_patch^2));
g = g/sum(g);
Gx = toeplitz([g(1),zeros(1,m-1)],[g,zeros(1,m-1)]);
Gy = toeplitz([g(1),zeros(1,n-1)],[g,zeros(1,n-1)]);


d = inf*ones(1,(2*w+1)*w +w+1);  %non-zero diagonals of W
d(1) = 0;
B = zeros(N,(2*w+1)*w +w+1);     %elements of the non-zero elements of W
B(:,1) = inf;

l = 2;
%compute all weights and store them in W
for ii = -w:w
    for jj = 0:w
       if ii+jj*m>0

         K = zeros(m,n);
         x_1 = max(1,1-ii):min(ExtImgDim(1),m+2*p-ii);
         x_2 = max(1,1+ii):min(ExtImgDim(1),m+2*p+ii);
         y_1 = max(1,1-jj):min(ExtImgDim(2),n+2*p-jj);
         y_2 = max(1,1+jj):min(ExtImgDim(2),n+2*p+jj);
         Q = M.dist(reshape(F_mirror(:,x_1,y_1),[M.ItemSize,length(x_1),length(y_1)]),...
               reshape(F_mirror(:,x_2,y_2),[M.ItemSize,length(x_2),length(y_2)]));
 
         K(max(1,1+ii):min(end,end+ii),max(1,1+jj):min(end,end+jj)) = Gx(1:end-abs(ii),1:end-abs(ii))*Q*Gy(1:end-abs(jj),1:end-abs(jj))';
     
         d(l) = ii+jj*m;
         B(:,l) = K(:);
                
         l = l+1;
       end
    end
end

W = spdiags(B,d,N,N);
W = W + W';

W = -W/(2*sigma_weights^2);
W = spfun(@exp,W);

[weights,ind_i] = sort(W,'descend');
weights = reshape(weights(1:k,:),[],1);
ind_i = reshape(ind_i(1:k,:),[],1);
ind_j = kron((1:N)',ones(k,1));
W = sparse(ind_i,ind_j,weights,N,N);
if center_element == 'm'
    W = W+spdiags(max(W)',0,N,N);
else
    W = W+speye(N);
end
    
W = max(W,W');
W = spdiags(1./sum(W,2),0,N,N)*W;
end
