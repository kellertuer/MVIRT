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
    VectorXd T = X-Y;
    VectorXd V = Y - scp*X;
    double nv = V.norm();
    if (T.norm() > std::numeric_limits<double>::denorm_min() && nv> std::numeric_limits<double>::denorm_min())
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
    if (norm_dir==0)
        return V; //X=Y -> no transport necessary
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