function f = ArtificialSPD2Image(pts,s)
% ArtificialSPDImage(pts)
% Create an Image of SPD-valued points of size pts x pts.
% Opitonal: s [4] says how aniisotrope the ellipses get
%
% ---
% ManImRes 1.0, J. Persch ~ 2016-06-30
if nargin<2
    s = 8;
end
[X,Y] =  meshgrid(-floor((pts-1)/2):ceil((pts-1)/2),-floor((pts-1)/2):ceil((pts-1)/2));

angles = atan2(X,Y);
radii = sqrt(Y.^2+X.^2);

U = zeros(2,2,pts,pts);
U(1,1,:,:) = cos(angles);
U(2,2,:,:) = cos(angles);
U(1,2,:,:) = sin(angles);
U(2,1,:,:) = -sin(angles);
circ = zeros(2,2,pts,pts);
for i = 1:pts
    for j = 1:pts
        circ(:,:,i,j) = U(:,:,i,j)'*[1,0;0,1+s*radii(i,j)/pts]*U(:,:,i,j);   
    end
end
f = circ;
f(:,:,1:floor(pts/2),1:floor((pts+1)/2)) = circ(:,:,1:floor(pts/2),ceil((pts+1)/2):end);
f(:,:,1:floor(pts/2),ceil((pts+1)/2):end) = circ(:,:,1:floor(pts/2),1:floor((pts+1)/2));