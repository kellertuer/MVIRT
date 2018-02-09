% ====
% Illustration of the gradient of the second order difference from the
% mid point model for three points on the sphere
%
% ---
% Manifold-valued Image Restoration Toolbox 1.2 | R. Bergmann | 2018-02-09
start = pwd;
if ~isempty(fileparts(which(mfilename)))
    cd(fileparts(which(mfilename)));
end
run('../initMVIRT.m')

% Lets take three points on S2 excluding strange special cases (hence
% moving y a little
M = Sn(2);
x = [1;0;0];
z = [0;1;0];
y = M.geopoint([0;0;1], M.geopoint(x,z,0.5),0.05);
c = M.midPoint(x,z);
% the distance from y to c is the second order difference in the mid point
% model
M.dist(y,c)

% we take the inner differential (for x and z)
inner = -M.log(c,y)./M.dist(c,y);
% and evaluate gradX and gradZ using adjoint jacobian, gradY is just a log
gradX = M.AdjDxGeo(x,z,0.5,inner);
gradY = - M.log(y,c)./M.dist(c,y);
gradZ = M.AdjDyGeo(x,z,0.5,inner);

%collect everything to plot
pts = cat(2,x,y,z);
vecs = -cat(2,gradX,gradY,gradZ);

%first the vecors
quiver3(pts(1,:),pts(2,:),pts(3,:),vecs(1,:),vecs(2,:),vecs(3,:));
hold on
%Then the 4 points
plotS2(cat(2,x,y,z,c),'o','MarkerSize',3);
%xz for help where the mid points at
plotS2(M.geodesic(x,z,'t',0:1/100:1),'-');
%Dashed line indicating the length computed above
plotS2(M.geodesic(y,c,'t',0:1/100:1),'--b');
hold off

%we do a step in gradient direction, just with step size 1
xN = M.exp(x,-gradX);
yN = M.exp(y,-gradY);
zN = M.exp(z,-gradZ);

% and its reduced
M.dist(yN,M.midPoint(xN,zN))

hold on
% as one can see, printing them again into the sphere.
plotS2(cat(2,xN,yN,zN,M.midPoint(xN,zN)),'o','MarkerSize',3);
plotS2(M.geodesic(xN,zN,'t',0:1/100:1),'-');
plotS2(M.geodesic(yN,M.midPoint(xN,zN),'t',0:1/100:1),'--k');
hold off


