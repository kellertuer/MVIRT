/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Small helper part for all SPD C++-functions used in the wrappers
 *
 * ---
 * ManImRes, R. Bergmann ~ 2015-04-12
 */
#include "manifoldSPD.h"
#include "mex.h"
#include "iostream"
using namespace std;
MatrixXd mSPDExp(MatrixXd X, MatrixXd V) {
    return mSPDExp(X,V,1.0);
}
MatrixXd mSPDExp(MatrixXd X, MatrixXd V,double t) {
    int ItemSize = X.rows();
    int i;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return MatrixXd::Zero(1,1);
    MatrixXd Y(ItemSize,ItemSize), U_i(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize), lY(ItemSize,ItemSize);
    if ( (X.isZero()) || (V.isZero()) || (t==0) )
        Y = X;
    else {
        svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
        S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues().cwiseSqrt());  
        U = svd.matrixU();
        Xsqrt = U*S*U.transpose();
        for (i=0; i<ItemSize; i++)
            S(i,i) = 1/S(i,i);
        lY = U*S*U.transpose()*(t*V)*U*S*U.transpose(); //A in the matlab code: sX\Vt/sX
        eig = EigenSolver<MatrixXd>(0.5*(lY+lY.transpose()), true);
        S = DiagonalMatrix<double,Dynamic,Dynamic>(eig.eigenvalues().real());
        U = eig.eigenvectors().real();
        for (i=0; i<ItemSize; i++)
            S(i,i) = exp(S(i,i));
        U_i = U*U.transpose()-MatrixXd::Identity(ItemSize, ItemSize);
        if(U_i.isZero()){
            U_i = U.transpose();
        }
        else{
            U_i = U.inverse();
        }
        Y = Xsqrt*U*S*U_i*Xsqrt;
    }
    return  0.5*(Y+Y.transpose());
}
MatrixXd mSPDLog(MatrixXd X, MatrixXd Y) {
    int ItemSize = X.rows();
    int i;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return MatrixXd::Zero(1,1);
    MatrixXd V(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize), lV(ItemSize,ItemSize);
    if ( (X.isZero()) || (Y.isZero()) )
        return MatrixXd::Zero(ItemSize,ItemSize);
    svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues().cwiseSqrt());
    U = svd.matrixU();
    Xsqrt = U*S*U.transpose();
    // % V = Xsqrt * [ logm(U*diag(1./sqrt(diag(D)))*U'*Y*U*diag(1./sqrt(diag(D)))*U') * Xsqrt
    for (i=0; i<ItemSize; i++)
        S(i,i) = 1/S(i,i);
    lV = U*S*U.transpose()*Y*U*S*U.transpose();
    // Is the following stable enough or do we need another svd?
    svd.compute(0.5*(lV+lV.transpose()), ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
            S(i,i) = log(S(i,i));
    U = svd.matrixU();
    V = Xsqrt*U*S*U.transpose()*Xsqrt;
    return 0.5*(V+V.transpose());
}
double mSPDDist(MatrixXd X, MatrixXd Y) {
    int ItemSize = X.rows();
    int i;
    GeneralizedEigenSolver<MatrixXd> ges;
    ArrayXXd S;
    if (ItemSize != X.cols())
        return 0;
    if ( (X.isZero()) || (Y.isZero()) )
        return 0;
    if  ((X.array()-Y.array()).isZero() )
        return 0;    
    ges.compute(X,Y); 
    S = ges.eigenvalues().transpose().cwiseAbs();
    S = log(S);
    return sqrt((S*S).sum()) ;
}
double mSPDDot(MatrixXd X, MatrixXd V, MatrixXd W) {
    int ItemSize = X.rows();
    int i;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return 0;
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize);
    if ( X.isZero() )
        return 0;
    svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
            S(i,i) = 1/S(i,i);
    U = svd.matrixU();
    S = V*U*S*U.transpose()*W*U*S*U.transpose(); //V*Xinv*W*Xinv
    return S.trace();
}
MatrixXd mSPDParallelTransport(MatrixXd X, MatrixXd Y, MatrixXd V) {
    return mSPDParallelTransport(X,Y,V,1.0);
}
MatrixXd mSPDParallelTransport(MatrixXd X, MatrixXd Y, MatrixXd V, double t) {
    int ItemSize = X.cols();
    int i;
    if (V.isZero())
        return V;
    MatrixXd S(ItemSize,ItemSize),U(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize);
    MatrixXd lX(ItemSize,ItemSize),lY(ItemSize,ItemSize),lV(ItemSize,ItemSize), lW(ItemSize,ItemSize);
    JacobiSVD<MatrixXd> svd(X, ComputeFullU);
    // Create Xsqrt
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues().cwiseSqrt());
    U = svd.matrixU();
    //for (i=0; i<ItemSize; i++)
    //S(i,i) = sqrt(S(i,i));
    Xsqrt = U*S*U.transpose();
    for (i=0; i<ItemSize; i++)
        S(i,i) = 1/S(i,i);
    lX = U*S*U.transpose(); // Save sqrtX^-1 in lX
    lV = lX*V*lX; // this is A from the matlab code
    lY = lX*Y*lX; //B in the Matlab code
    JacobiSVD<MatrixXd> svd2(0.5*(lY+lY.transpose()), ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd2.singularValues());
    for (i=0; i<ItemSize; i++)
        S(i,i) = log(S(i,i));
    U = svd2.matrixU();
    for (i=0; i<ItemSize; i++)
        S(i,i) = exp(0.5*S(i,i));
    lW = Xsqrt*(U*S*U.transpose())*0.5*(lV+lV.transpose())*(U*S*U.transpose())*Xsqrt;
    return lW;
}
// OldMean for better understanding of the one below
MatrixXd OldmSPDMean(MatrixXd *F, double *W, double E, double I,size_t L){
    MatrixXd X,V,Xold;
    X = F[0];
    int ItemSize = X.cols();
    int iter =  0;
    double step_size = 1.;
    double change = 1+E;
    double max_dist = 0.;
    mwSize i,j,k;
    MatrixXd *V_;
    V_ = new MatrixXd[L];
    for(i =1;i<L;i++){
        max_dist = max(mSPDDist(F[0],F[i]),max_dist);
    }
    if(max_dist > 2){
        step_size = 1/max_dist;
        E = E*step_size;
        I = I/step_size;
    }
    while( ( iter < I) && (change > E)){
        Xold = X;  
        V = MatrixXd::Zero(ItemSize,ItemSize);
        // Sum the Tangentvectors up according to their weights
        for (i = 0; i < L; i++)
            V = V+W[i]*mSPDLog(X,F[i]); 
        // go in the right direction
        X = mSPDExp(X,V,step_size);        
        change = mSPDDist(X,Xold);
        iter++;        
    }        
    //mexPrintf("\n%d, %f",iter,change);
    delete[] V_;
    V_ = NULL;
    return X;
}
// Not good to read but 25-30% better performance as old one, does nor use 
// Log/Exp from the manifold, as the last multiplication and the several 
// calculations of Xsqrt/Xinvsqrt can be saved. OldMean still included
MatrixXd mSPDMean(MatrixXd *F, double *W, double E, double I,size_t L){
    MatrixXd X,V,Xold;
    X = F[0];
    int ItemSize = X.cols();
    int iter =  0;        
    double step_size = 1.;
    double change = 1+E;
    double max_dist = 0.;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize),U_i(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize),Xinvsqrt(ItemSize,ItemSize), lV(ItemSize,ItemSize);
    mwSize i,j,k;
    MatrixXd *V_;
    V_ = new MatrixXd[L];
    
    for(i =1;i<L;i++){
        max_dist = max(mSPDDist(F[0],F[i]),max_dist);
    }
    if(max_dist > 2){
        step_size = 1/max_dist;
        E = E*step_size;
        I = I/step_size;
    }
    
    while( ( iter < I) && (change > E)){
        Xold = X;  
        V = MatrixXd::Zero(ItemSize,ItemSize);
        // Calculations of Xsqrt, Xinvsqrt
        svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
        S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues().cwiseSqrt());        
        U = svd.matrixU();
        Xsqrt = U*S*U.transpose();
        for (i=0; i<ItemSize; i++)
            S(i,i) = 1/S(i,i);
        Xinvsqrt =  U*S*U.transpose();
        // Sum the Tangentvectors up according to their weights
        // (no multiplication with Xsqrt)
        for (i = 0; i < L; i++){
            lV = Xinvsqrt*F[i]*Xinvsqrt;    
            svd.compute(0.5*(lV+lV.transpose()), ComputeFullU);
            S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues().unaryExpr<double(*)(double)>(&std::log));
            U = svd.matrixU();
            lV = U*S*U.transpose();
            V = V+W[i]*lV;             
        }
        // go in the right direction by calcculating Xsqrt*exp(V)*Xsqrt
        eig = EigenSolver<MatrixXd>(0.5*(V+V.transpose()), true);
        S = DiagonalMatrix<double,Dynamic,Dynamic>(eig.eigenvalues().real().unaryExpr<double(*)(double)>(&std::exp));
        U = eig.eigenvectors().real();        
        U_i = U*U.transpose()-MatrixXd::Identity(ItemSize, ItemSize);
        if(U_i.isZero()){
            U_i = U.transpose();
        }
        else{
            U_i = U.inverse();
        }
        V = Xsqrt*U*S*U_i*Xsqrt;
        
        X = 0.5*(V+V.transpose());        
        change = mSPDDist(X,Xold);
        iter++;        
    }        
    //mexPrintf("\n%d, %f",iter,change);
    delete[] V_;
    V_ = NULL;
    return X;     
}
MatrixXd mSPDGeo(MatrixXd X, MatrixXd Y,double t) {
    int ItemSize = X.rows();
    int i;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return MatrixXd::Zero(1,1);
    MatrixXd V(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), lV(ItemSize,ItemSize),U(ItemSize,ItemSize), U_i(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize);
    if ( (X.isZero()) || (Y.isZero()) )
        return MatrixXd::Zero(ItemSize,ItemSize);
    svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues().cwiseSqrt());
    U = svd.matrixU();
    Xsqrt = U*S*U.transpose();
    // % V = Xsqrt * [ logm(U*diag(1./sqrt(diag(D)))*U'*Y*U*diag(1./sqrt(diag(D)))*U') * Xsqrt
    for (i=0; i<ItemSize; i++)
        S(i,i) = 1/S(i,i);
    lV = U*S*U.transpose()*Y*U*S*U.transpose();
    // Is the following stable enough or do we need another svd?
    svd.compute(0.5*(lV+lV.transpose()), ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
            S(i,i) = pow(S(i,i),t);
    U = svd.matrixU();
    Y = Xsqrt*U*S*U.transpose()*Xsqrt;
    
    return  0.5*(Y+Y.transpose());
}
