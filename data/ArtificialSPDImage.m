function f = ArtificialSPDImage(pts, jumpSize )
% ArtificialSPDImage(pts,jumpSize)
% Create an Image of SPD-valued points of size pts x pts having jumps in
% their eigenvalues along (across) the diagonal and the vertical and
% horizontal middle line.
%
% ---
% MVIRT 1.0, R. Bergmann ~ 2015-04-15

    t = linspace(0,1-1/pts,pts);
    t1 = abs(2*pi*t-pi);
    t2 = pi*t;
    t3 = linspace(0,3*(1-1/pts),2*pts);
    f = zeros(3,3,pts,pts);
    for x=1:pts
        for y=1:pts
            B = [1,0,0; 0,cos(t2(x)), -sin(t2(x));0,sin(t2(x)),cos(t2(x))];
            A = [cos(t1(y)),-sin(t1(y)),0;sin(t1(y)),cos(t1(y)),0;0,0,1];
            C = [cos(t1(mod(y-x,pts)+1)),0,-sin(t1(mod(y-x,pts)+1));0,1,0;sin(t1(mod(y-x,pts)+1)),0,cos(t1(mod(y-x,pts)+1))];
            f(:,:,x,y) = ...
                A*B*C*diag([...
                1+jumpSize/2*((y+x)>pts),...%
                1+t3(x+y)-jumpSize*(y>pts/2),...%
                4-t3(x+y)+jumpSize*(x>pts/2)
            ])*C'*B'*A';
         end
    end
end

