/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the dot product in tangential planes of SPD matrices
 * 
 * ---
 * Manifold-valued Image Restoration Toolbox 1.0
 * R. Bergmann ~ 2015-04-12
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
    const mxArray *IXA, *IVA, *IWA; //Just names
    mxArray *OdA;
    double *Od, *IX, *IV, *IW;
    /*For-Loops*/
    mwSize i,j,k;
    /* Input Names */
    if (nrhs != 3)
        mexErrMsgTxt("Wrong number of inputs; three arrays of matrices required.");
    IXA = prhs[0];
    IVA = prhs[1];
    IWA = prhs[2];
    IX = mxGetPr(IXA); // Base Points
    IV = mxGetPr(IVA); // Second Points
    IW = mxGetPr(IWA); // Second Points
    size_t Xs = mxGetNumberOfDimensions(IXA);
    const mwSize *Xn = mxGetDimensions(IXA);
    size_t Vs = mxGetNumberOfDimensions(IVA);
    const mwSize *Vn = mxGetDimensions(IVA);
    size_t Ws = mxGetNumberOfDimensions(IWA);
    const mwSize *Wn = mxGetDimensions(IWA);
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
        mexErrMsgTxt("Matrices have to be square.");
    if ( (Xs!=Vs) || (Xs!=Ws) )
        mexErrMsgTxt("Matrix arrays have to be of same dimension.");
    for (i=0; i<Xs; i++) {
        if ( (Vn[i]!=Xn[i]) || (Wn[i]!=Xn[i]) )
            mexErrMsgTxt("Matrices of Matrix arrays are not of same size.");
    }
    /* Create Output */
    if (Xs==2) {
        mwSize dims[1];
        dims[0]=1;
        plhs[0] = mxCreateNumericArray(1,dims, mxDOUBLE_CLASS, mxREAL);
        Od = mxGetPr(plhs[0]);
    }
    else if (Xs>2) {
         mwSize* dims = (mwSize*)mxCalloc(Xs-2,sizeof(mwSize));
       /* mwSize dims[Xs-2];*/
        for (i=2; i<Xs; i++) {
            dims[i-2] = Xn[i];
        }
        plhs[0] = mxCreateNumericArray(Xs-2,dims, mxDOUBLE_CLASS, mxREAL);
        Od = mxGetPr(plhs[0]);
        mxFree(dims);
    }
    // Initialize result: V matrix array
    MatrixXd lX(ItemSize,ItemSize), lV(ItemSize,ItemSize), lW(ItemSize,ItemSize);
    for (i=0; i<numel; i++) {//for all matrices
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract local copies of submatrices
                lX(j,k) = IX[j + ItemSize*k + ItemSize*ItemSize*i];
                lV(j,k) = IV[j + ItemSize*k + ItemSize*ItemSize*i];
                lW(j,k) = IW[j + ItemSize*k + ItemSize*ItemSize*i];
            }
        }
        Od[i] = mSPDDot(lX,lV,lW);
    }
    // Eigen clears itself, so there's just dims left
}
