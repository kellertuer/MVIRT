#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <limits>
#include "manifolds/manifoldSn.h"

using namespace std;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])   {
    /* Input */
    const mxArray *IXA, *IYA, *IZA; //Just names
    double *IX, *IY, *IZ, *OG;
    /* variables */
    if (nrhs != 3)
        mexErrMsgTxt("Wrong number of inputs; three arrays of vectors required.");
    IXA = prhs[0];
    IYA = prhs[1];
    IZA = prhs[2];
    IX = mxGetPr(IXA); // Base Points
    IY = mxGetPr(IYA); // Base Points
    IZ = mxGetPr(IZA); // Base Points
    size_t Xs = mxGetNumberOfDimensions(IXA);
    const mwSize *Xn = mxGetDimensions(IXA);
    size_t Ys = mxGetNumberOfDimensions(IYA);
    const mwSize *Yn = mxGetDimensions(IYA);
    size_t Zs = mxGetNumberOfDimensions(IZA);
    const mwSize *Zn = mxGetDimensions(IZA);
    mwSize numel=1,i=0,j=0,ItemSize = Xn[0];
    if (Xs>1) {
        for (i=1; i<Xs; i++) {
            numel *= Xn[i];
        }
    }    
    /* get dimensions of the input vectors */
    /*
     * Checks
     */
    if ( (Xs!=Ys) || (Xs!=Zs) )
        mexErrMsgTxt("Vector arrays have to be of same dimension.");
    for (i=0; i<Xs; i++) {
        if ( (Xn[i]!=Yn[i]) || (Xn[i]!=Zn[i]) )
            mexErrMsgTxt("Vector arrays are not of same size.");
    }
    /* Create Output */
        plhs[0] = mxCreateNumericArray(Xs,Xn, mxDOUBLE_CLASS, mxREAL);
        OG = mxGetPr(plhs[0]);
    //Compute result
    double nv=0;
    VectorXd X = VectorXd::Zero(ItemSize), Y = VectorXd::Zero(ItemSize), Z = VectorXd::Zero(ItemSize), G = VectorXd::Zero(ItemSize);
    for (i=0; i<numel; i++) {
        for (j=0; j<ItemSize; j++) {
            X(j) = IX[j+ItemSize*i];
            Y(j) = IY[j+ItemSize*i];
            Z(j) = IZ[j+ItemSize*i];
        }
        G = mSnGrad_X_D2(X,Y,Z);
        for (j=0; j<ItemSize; j++) {
            OG[j + ItemSize*i] = G(j);
        }
    }
}