/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the gradient in X of the second order differences of SPD matrices
 * 
 * ---
 * ManImRes, R. Bergmann ~ 2015-04-12
 */
#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <algorithm>    // std::max
#include <Eigen/Dense>
#include "manifolds/manifoldSPD.h"
using namespace std;
using namespace Eigen;
        
size_t ItemSize;  /* size of matrix */
size_t N;         /* size of matrix */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])   {
    /* variables (IX, IY input) */
    const mxArray *IXA, *IYA, *IZA; //Just names
    double *OG, *IX, *IY, *IZ;
    /*For-Loops*/
    mwSize i,j,k;
    /* Input Names */
    if (nrhs != 3)
        mexErrMsgTxt("Wrong number of inputs; three arrays of matrices required.");
    IXA = prhs[0];
    IYA = prhs[1];
    IZA = prhs[2];
    IX = mxGetPr(IXA); // Base Points
    IY = mxGetPr(IYA); // Second Points
    IZ = mxGetPr(IZA); // Second Points
    size_t Xs = mxGetNumberOfDimensions(IXA);
    const mwSize *Xn = mxGetDimensions(IXA);
    size_t Ys = mxGetNumberOfDimensions(IYA);
    const mwSize *Yn = mxGetDimensions(IYA);
    size_t Zs = mxGetNumberOfDimensions(IZA);
    const mwSize *Zn = mxGetDimensions(IZA);
    mwSize numel=1;
    if (Xs<2)
        mexErrMsgTxt("No Matrix X given, only an array.");
    else if (Xs>2) {
        for (i=2; i<Xs; i++) {
            numel *= Xn[i];
        }
    }    
    /* get dimensions of the input matrix */
    ItemSize = Xn[0];
    /*
     * Checks
     */
    if (ItemSize != Xn[1])
        mexErrMsgTxt("Matrices have to be suqare.");
    if ( (Xs!=Ys) || (Xs!=Zs) )
        mexErrMsgTxt("Matrix arrays have to be of same dimension.");
    for (i=0; i<Xs; i++) {
        if ( (Yn[i]!=Xn[i]) || (Zn[i]!=Xn[i]) )
            mexErrMsgTxt("Matrices of Matrix arrays are not of same size.");
    }
    /* Create Output */
    plhs[0] = mxCreateNumericArray(Xs, Xn, mxDOUBLE_CLASS, mxREAL);
    OG = mxGetPr(plhs[0]);
    // Initialize result: V matrix array
    MatrixXd lX(ItemSize,ItemSize), lY(ItemSize,ItemSize), lZ(ItemSize,ItemSize),lG(ItemSize,ItemSize);
    for (i=0; i<numel; i++) {//for all matrices
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract local copies of submatrices
                lX(j,k) = IX[j + ItemSize*k + ItemSize*ItemSize*i];
                lY(j,k) = IY[j + ItemSize*k + ItemSize*ItemSize*i];
                lZ(j,k) = IZ[j + ItemSize*k + ItemSize*ItemSize*i];
            }
        }
        lG = mSPDGrad_X_D2(lX,lY,lZ);
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract copy back to result V
                OG[j + ItemSize*k + ItemSize*ItemSize*i] = lG(j,k);
            }
        }

    }
    // Eigen clears itself, so there's just dims left
}