/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Small helper part for all Hn C++-functions used in the wrappers
 *
 * ---
 * ManImRes, R. Bergmann ~ 2015-04-22
 */
#include "manifoldHn.h"
#include "mex.h"

VectorXd mHnExp(VectorXd X, VectorXd V) {
    return mHnExp(X,V,1.0);
}
VectorXd mHnExp(VectorXd X, VectorXd V,double t) {
    int ItemSize = X.size();
    int i;
    double n = sqrt(max(mHnDot(V,V),0.0));
    if (n < std::numeric_limits<double>::denorm_min())
        return X;
    else {
        VectorXd Y = X * cosh(t*n) + V * ((sinh(t*n)/n));
        Y = Y * (1/sqrt(-mHnDot(Y,Y)));
        return Y;
    }
}
VectorXd mHnLog(VectorXd X, VectorXd Y) {
    double scp = mHnDot(X,Y);
    VectorXd V = Y + scp*X;
    double nv = max(sqrt(scp*scp-1.0),0.0);
    if (nv>0)
        return V*acosh(max(-scp,1.0))/nv;
    else
        return VectorXd::Zero(X.size());
}
double mHnDist(VectorXd X, VectorXd Y)
{
    return acosh(max(-mHnDot(X,Y),1.0));
}
double mHnDot(VectorXd V, VectorXd W)
{
    int i=0;
    double d=0;
    int n = V.size()-1;
    for (i=0; i<n; i++)
        d = d+V(i)*W(i);
    d = d-V(n)*W(n);
    return d;
}
VectorXd mHnMean(VectorXd *F, double *W, double E, double I,size_t L){
    VectorXd X,V,Xold;
    X = F[0];
    int ItemSize = X.size();
    int iter =  0;
    double step_size = 1.;
    double max_dist = 0.;
    double change = 1+E;    
    mwSize i,j,k;
    for(i =1;i<L;i++){
        max_dist = max(mHnDist(F[0],F[i]),max_dist);
    }
    if(max_dist > 2){
        step_size = 1/max_dist;
        E = E*step_size;
        I = I/step_size;
    }
    while( ( iter < I) && (change > E)){
        Xold = X;  
        V = VectorXd::Zero(ItemSize);
        // Sum the Tangentvectors up according to their weights
        for (i = 0; i < L; i++)
        {  
            V = V+W[i]*mHnLog(X,F[i]); 
        }
        // go in the right direction
        X = mHnExp(X,V,step_size);        
        change = mHnDist(X,Xold);
        iter++;        
       // mexPrintf("\n%d, %3.16f",iter,change);
    }        
    return X;
}