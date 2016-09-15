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
    MatrixXd Y(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize), lY(ItemSize,ItemSize);
    if ( (X.isZero()) || (V.isZero()) || (t==0) )
        Y = X;
    else {
        svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
        S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
        for (i=0; i<ItemSize; i++)
            S(i,i) = sqrt(S(i,i));
        U = svd.matrixU();
        Xsqrt = U*S*U.transpose();
        for (i=0; i<ItemSize; i++)
            S(i,i) = 1/S(i,i);
        lY = U*S*U.transpose()*(t*V)*U*S*U.transpose(); //A in the matlab code: sX\Vt/sX
        eig = EigenSolver<MatrixXd>(0.5*(lY+lY.transpose()), true);
        S = eig.pseudoEigenvalueMatrix();
        U = eig.pseudoEigenvectors();
        for (i=0; i<ItemSize; i++)
            S(i,i) = exp(S(i,i));
        Y = Xsqrt*U*S*U.transpose()*Xsqrt;
    }
    return Y;
}
MatrixXd mSPDExpAtId(MatrixXd X, MatrixXd V,double t) {
    int ItemSize = X.rows();
    int i;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return MatrixXd::Zero(1,1);
    MatrixXd Y(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize), lY(ItemSize,ItemSize);
    if ( (X.isZero()) || (V.isZero()) || (t==0) )
        Y = X;
    else {
        svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
        S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
        for (i=0; i<ItemSize; i++)
            S(i,i) = sqrt(S(i,i));
        U = svd.matrixU();
        Xsqrt = U*S*U.transpose();
        for (i=0; i<ItemSize; i++)
            S(i,i) = 1/S(i,i);
        lY = t*V; 
        eig = EigenSolver<MatrixXd>(0.5*(lY+lY.transpose()), true);
        S = eig.pseudoEigenvalueMatrix();
        U = eig.pseudoEigenvectors();
        for (i=0; i<ItemSize; i++)
            S(i,i) = exp(S(i,i));
        Y = Xsqrt*U*S*U.transpose()*Xsqrt;
    }
    return Y;
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
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
        S(i,i) = sqrt(S(i,i));
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
    return Xsqrt*U*S*U.transpose()*Xsqrt;
}
MatrixXd mSPDLogAtId(MatrixXd X, MatrixXd Y) {
    int ItemSize = X.rows();
    int i;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return MatrixXd::Zero(1,1);
    MatrixXd V(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize), XsqrtInv(ItemSize,ItemSize), lV(ItemSize,ItemSize);
    if ( (X.isZero()) || (Y.isZero()) )
        return MatrixXd::Zero(ItemSize,ItemSize);
    svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
        S(i,i) = 1/sqrt(S(i,i));
    U = svd.matrixU();
    XsqrtInv = U*S*U.transpose();
    // % V = Xsqrt * [ logm(U*diag(1./sqrt(diag(D)))*U'*Y*U*diag(1./sqrt(diag(D)))*U') * Xsqrt
    for (i=0; i<ItemSize; i++)
        S(i,i) = 1/S(i,i);
    lV = XsqrtInv*Y*XsqrtInv;
    // Is the following stable enough or do we need another svd?
    svd.compute(0.5*(lV+lV.transpose()), ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
            S(i,i) = log(S(i,i));
    U = svd.matrixU();
    return XsqrtInv*U*S*U.transpose()*XsqrtInv;
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
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    U = svd.matrixU();
    for (i=0; i<ItemSize; i++)
        S(i,i) = sqrt(S(i,i));
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
    lY = U*S*U.transpose(); //temp for 3rd svd
    EigenSolver<MatrixXd> eig(0.5*0.5*(lY+lY.transpose()), true);// 0.5 from formula and 0.5 for stability
    S = eig.pseudoEigenvalueMatrix();
    U = eig.pseudoEigenvectors();
    for (i=0; i<ItemSize; i++)
        S(i,i) = exp(S(i,i));
    lW = Xsqrt*(U*S*U.transpose())*0.5*(lV+lV.transpose())*(U*S*U.transpose())*Xsqrt;
    return lW;
}
MatrixXd mSPDGrad_X_D2(MatrixXd X, MatrixXd Y, MatrixXd Z) {
    int ItemSize = X.rows();
    if (ItemSize!=X.cols())
        return MatrixXd::Zero(1,1);
    int i=0,j=0,s=0;
    int n = ItemSize*(ItemSize+1)/2;
    double r=0, alpha=0,lambda=0;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return MatrixXd::Zero(1,1);
    MatrixXd G(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize);
    MatrixXd V(ItemSize,ItemSize), M(ItemSize,ItemSize);
    MatrixXd Wx(ItemSize,ItemSize), Wm(ItemSize,ItemSize), R(ItemSize,ItemSize);
    if ( (X.isZero() || Y.isZero() || Z.isZero() ) )
        return MatrixXd::Zero(ItemSize,ItemSize);
    M = mSPDExp(X,mSPDLog(X,Z)/2.0); //midpoint
    //Compute one Jacobian Eigenframe in X and transport it to M
    //...and its lambdas
    V = mSPDLog(X,Z); //Starting point X->Z
    svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
        S(i,i) = sqrt(S(i,i)); //X^{-1/2}
    U = svd.matrixU();
    Xsqrt = U*S*U.transpose();
    for (i=0; i<ItemSize; i++)
        S(i,i) = 1/S(i,i);
    //called VE in matlab (loclEigenFrame)
    // X^{-1/2}*V*X^{-1/2}
    V = U*S*U.transpose()*V*U*S*U.transpose();
    eig = EigenSolver<MatrixXd> (0.5*(V+V.transpose()),true);
    S = eig.pseudoEigenvalueMatrix();
    U = eig.pseudoEigenvectors();
    // Compute all basis elements
    R = mSPDLog(M,Y);
    r = sqrt(mSPDDot(M,R,R));
    // Init G;
    G = MatrixXd::Zero(ItemSize,ItemSize);
    for (i=0; i<ItemSize; i++) {
        for (j=i; j<ItemSize; j++) {
            //Compute sth eigenvalue based on i&j
            lambda = abs((S(i,i) - S(j,j)))/2*mSPDDist(X,Z);
            //Compute corresponding basis element (matrix Wx)
            if (i==j)
                Wx = 0.5*( U.col(i)*U.col(j).transpose() + U.col(j)*U.col(i).transpose() );
            else
                Wx = sqrt(0.5)*( U.col(i)*U.col(j).transpose() + U.col(j)*U.col(i).transpose() );
            Wx = Xsqrt*Wx*Xsqrt;
            // Basis Vector in X, transport it to M
            Wm = mSPDParallelTransport(X, M, Wx);
            if ( r==0 )
                alpha = 0;
            else
                alpha = mSPDDot(M,R,Wm)/r;
             if ( lambda>0 ) // eigenvalue > 0
                G = G + ( (sinh(lambda/2)/sinh(lambda)*alpha)*Wx);
             else // eigenvalue 0
                G = G + (0.5*alpha*Wx);
        }
    } //end running through all basis matrices W
    return G;
}

MatrixXd mSPDGrad_X_D2_Sq(MatrixXd X, MatrixXd Y, MatrixXd Z) {
    int ItemSize = X.rows();
    if (ItemSize!=X.cols())
        return MatrixXd::Zero(1,1);
    int i=0,j=0,s=0;
    int n = ItemSize*(ItemSize+1)/2;
    double r=0, alpha=0,lambda=0;
    JacobiSVD<MatrixXd> svd;
    EigenSolver<MatrixXd> eig;
    if (ItemSize != X.cols())
        return MatrixXd::Zero(1,1);
    MatrixXd G(ItemSize,ItemSize);
    MatrixXd S(ItemSize,ItemSize), U(ItemSize,ItemSize), Xsqrt(ItemSize,ItemSize);
    MatrixXd V(ItemSize,ItemSize), M(ItemSize,ItemSize);
    MatrixXd Wx(ItemSize,ItemSize), Wm(ItemSize,ItemSize), R(ItemSize,ItemSize);
    if ( (X.isZero() || Y.isZero() || Z.isZero() ) )
        return MatrixXd::Zero(ItemSize,ItemSize);
    M = mSPDExp(X,mSPDLog(X,Z)/2.0); //midpoint
    //Compute one Jacobian Eigenframe in X and transport it to M
    //...and its lambdas
    V = mSPDLog(X,Z); //Starting point X->Z
    svd = JacobiSVD<MatrixXd>(X, ComputeFullU);
    S = DiagonalMatrix<double,Dynamic,Dynamic>(svd.singularValues());
    for (i=0; i<ItemSize; i++)
        S(i,i) = sqrt(S(i,i)); //X^{-1/2}
    U = svd.matrixU();
    Xsqrt = U*S*U.transpose();
    for (i=0; i<ItemSize; i++)
        S(i,i) = 1/S(i,i);
    //called VE in matlab (loclEigenFrame)
    // X^{-1/2}*V*X^{-1/2}
    V = U*S*U.transpose()*V*U*S*U.transpose();
    eig = EigenSolver<MatrixXd> (0.5*(V+V.transpose()),true);
    S = eig.pseudoEigenvalueMatrix();
    U = eig.pseudoEigenvectors();
    // Compute all basis elements
    R = mSPDLog(M,Y);
    r = sqrt(mSPDDot(M,R,R));
    // Init G;
    G = MatrixXd::Zero(ItemSize,ItemSize);
    for (i=0; i<ItemSize; i++) {
        for (j=i; j<ItemSize; j++) {
            //Compute sth eigenvalue based on i&j
            lambda = abs((S(i,i) - S(j,j)));
            //Compute corresponding basis element (matrix Wx)
            if (i==j)
                Wx = 0.5*( U.col(i)*U.col(j).transpose() + U.col(j)*U.col(i).transpose() );
            else
                Wx = sqrt(0.5)*( U.col(i)*U.col(j).transpose() + U.col(j)*U.col(i).transpose() );
            Wx = Xsqrt*Wx*Xsqrt;
            // Basis Vector in X, transport it to M
            Wm = mSPDParallelTransport(X, M, Wx);
            if ( r==0 )
                alpha = 0;
            else
                alpha = 2*mSPDDot(M,R,Wm);
             if ( lambda>0 ) // eigenvalue > 0
                G = G + ( (sinh(lambda/2)/sinh(lambda)*alpha)*Wx);
             else // eigenvalue 0
                G = G + (0.5*alpha*Wx);
            
        }
    } //end running through all basis matrices W
    return G;
}
MatrixXd mSPDMean(MatrixXd *F, double *W, double E, double I,size_t L){
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
