/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Small helper part for all SPD C++-functions used in the wrappers
 *
 * ---
 * Manifold Valued Image Restoration 1.0
 * R. Bergmann ~ 2015-04-22
 */
#include "manifoldSn.h"
#include "mex.h"

VectorXd mSnExp(VectorXd X, VectorXd V) {
    return mSnExp(X,V,1.0);
}
VectorXd mSnExp(VectorXd X, VectorXd V,double t) {
    int ItemSize = X.size();
    int i;
    double n = V.norm();
    if (n < std::numeric_limits<double>::denorm_min())
        return X;
    else {
        VectorXd Y = X * cos(t*n) + V * ((sin(t*n)/n));
        // Needed?
        Y = Y * (1/Y.norm());
        return Y;
    }
}
VectorXd mSnLog(VectorXd X, VectorXd Y) {
    double scp = max(min(X.dot(Y),1.0),-1.0);
    VectorXd V = Y - scp*X;
    double nv = V.norm();
    if (nv>0)
        return V*acos(scp)/nv;
    else
        return VectorXd::Zero(X.size());
}
double mSnDist(VectorXd X, VectorXd Y)
{
    return acos(max(min(X.dot(Y),1.0),-1.0));
}
VectorXd mSnGradX(VectorXd X, VectorXd Y, VectorXd Z) {
    int i=0;
    // 1. Compute mid point
    VectorXd M = mSnExp(X,mSnLog(X,Z),0.5);
    VectorXd R = mSnLog(M,Y);
    if (R.norm()==0) 
        return VectorXd::Zero(X.size());
    R = R/R.norm();
    // Compute ONBs of tangential Planes at at X and M
    //                V = this.TpMONB(m(:,proxpts),x(:,proxpts,3));
    //                W = this.TpMONB(x(:,proxpts,1), x(:,proxpts,3));
    VectorXd G = VectorXd::Zero(X.size());
    if (X.size()!=3)
        mexErrMsgTxt("Not yet implemented in C for Sn, n>2.");
    // On S2 we have 2 vectors
    Vector3d V = (Vector3d) mSnLog(M,Z);
    if (V.norm()>0)
        V = V/V.norm();
    Vector3d W = mSnLog(X,Z);
    if (W.norm()>0)
        W = W/W.norm();
    G = G+V.dot(R)*0.5*W; //alpha*0.5*W
    
    // These would have to be performed on a basis of the nullspace of [X V V V V ... V]
    V = V.cross((Vector3d)M); //No normlization needed, because both are units
    W = W.cross((Vector3d)X); //Same    
    G = G + V.dot(R)*1/(2*cos(mSnDist(X,Z)/2))*W; //alpha*1/(2cos(dist(X,Z)/2))(W without 1/2    
    return (VectorXd) G;
}