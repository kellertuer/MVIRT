function patch = whirlpatch(ps)
% S2WhirlPatch(pts)
%   constructs a small whirl of size pts x pts points of S2-data
%
% Input
%   pts: a number specifying the image size (pts x pts)
%
% Output
% f : an [3,pts,pts] size array, the S2 valued image
%
% ---
% Manifold-valued Image Restoration Toolbox 1.1
%  J. Persch, R. Bergmann ~ 2017-01-05
% see LICENSE.txt
patch = zeros(3,ps(1),ps(2));
[x,y] = meshgrid(-floor(ps(1)/2):floor(ps(2)/2),-floor(ps(1)/2):floor(ps(2)/2));
alpha = atan2(y,x);
beta = sqrt(x.^2+y.^2);
beta = beta/max(beta(:))*pi/3;
patch(1,:,:) = sin(alpha).*sin(pi/2-beta);
patch(2,:,:) = cos(alpha).*sin(pi/2-beta);
patch(3,:,:) = cos(pi/2-beta);
patch(:,floor((ps(1)+1)/2),floor((ps(2)+1)/2)) = -[0;0;1];
end