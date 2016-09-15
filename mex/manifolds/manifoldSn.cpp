/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Small helper part for all SPD C++-functions used in the wrappers
 *
 * ---
 * ManImRes, R. Bergmann ~ 2015-04-22
 */
#include "manifoldSn.h"
#include "mex.h"
#include <iostream>

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
VectorXd mSnParallelTransport(VectorXd X, VectorXd Y, VectorXd V) {
    VectorXd dir = mSnLog(X,Y);
    double norm_dir = dir.norm();
    dir = dir/norm_dir;
    return V - dir.dot(V)*(dir+mSnLog(Y,X)/norm_dir);
}
MatrixXd mSnTxM_ONB(VectorXd X, VectorXd Y) {
// internal function using X as base point and Log_xY as first direction in TxM to build the ONB
    int i;
    MatrixXd LGS = MatrixXd(2,X.size());
    LGS.row(0) = X;
    LGS.row(1) = mSnLog(X,Y);
    LGS.row(1) = 1/LGS.row(1).norm()*LGS.row(1);
    MatrixXd TxMONB = MatrixXd(X.size(),X.size()-1);
    if (LGS.row(1).norm()==0)
        return MatrixXd::Zero(X.size(),X.size()-1);
    TxMONB.col(0) = LGS.row(1);
    //The kernel should be X.size()-2-dimensional
    MatrixXd ker = LGS.fullPivLu().kernel();
    // mexPrintf("%d %d \n",ker.cols(),ker.rows());
    // Gram-Schmidt to form an ONB
    for (i=1; i< (X.size()-1);i++){
        TxMONB.col(i) = ker.col(i-1);
        for (int j=0;j<i;j++){
            TxMONB.col(i) -= TxMONB.col(j).dot(ker.col(i-1))*TxMONB.col(j);
        }
        TxMONB.col(i) = 1/TxMONB.col(i).norm()*TxMONB.col(i);
    }
    return TxMONB;
}
VectorXd mSnGrad_X_D2(VectorXd X, VectorXd Y, VectorXd Z) {
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
    if (X.size()==3) {
        // On S2 we have 2 vectors
        // (1) X->Z gets transported to M->Z
        Vector3d V = (Vector3d) mSnLog(M,Z);
        if (V.norm()>0)
            V = V/V.norm();
        Vector3d W = mSnLog(X,Z);
        if (W.norm()>0)
            W = W/W.norm();
        G = G+V.dot(R)*0.5*W; //alpha*0.5*W
        // (2) The orthogonal one
        // These would have to be performed on a basis of the nullspace of [X V V V V ... V]
        V = V.cross((Vector3d)M); //No normlization needed, because both are units
        W = W.cross((Vector3d)X); //Same    
        G = G + V.dot(R)*1/(2*cos(mSnDist(X,Z)/2))*W; //alpha*1/(2cos(dist(X,Z)/2))(W without 1/2
     } else {//Sn, n>2
        double alpha = 0;
        MatrixXd V = mSnTxM_ONB(M,Z);
        MatrixXd W = V;
        for (i=0; i<W.cols(); i++) {
            W.col(i) = mSnParallelTransport(M,X,V.col(i));
        }
        // MatrixXd W = mSnTxM_ONB(X,Z);
        for (i=0; i<V.cols(); i++) {
            //compute alpha
            VectorXd vT = V.col(i);
            alpha = vT.dot(R);
            if (i==0)
                alpha = alpha/2; //towards y
            else
                alpha = alpha/(2*cos(mSnDist(X,Z)/2));
            VectorXd wT = W.col(i);
            G += alpha*wT;
            //mexPrintf("%d:%f ",i,alpha);
        }
     }
    return (VectorXd) G;
}

VectorXd mSnMean(VectorXd *F, double *W, double E, double I,size_t L){
    VectorXd X,V,Xold;
    X = F[0];
    int ItemSize = X.size();
    int iter =  0;
    double step_size = 1.;
    double change = 1+E;    
    mwSize i,j,k;    
    while( ( iter < I) && (change > E)){
        Xold = X;  
        V = VectorXd::Zero(ItemSize);
        // Sum the Tangentvectors up according to their weights
        for (i = 0; i < L; i++)
        {  
            V = V+W[i]*mSnLog(X,F[i]); 
        }
        // go in the right direction
        X = mSnExp(X,V,step_size);        
        change = mSnDist(X,Xold);
        iter++;        
       // mexPrintf("\n%d, %3.16f",iter,change);
    }        
    return X;
}
