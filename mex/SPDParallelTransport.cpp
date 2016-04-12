/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the C++-Function to compute the parallel transport on 
 * symmetric posivive definte matrices.
 *
 * ---
 * Manifold-valued Image Restoration Toolbox 1.0
 * R. Bergmann ~ 2015-04-12
 */
#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "manifolds/manifoldSPD.h"

using namespace std;
using namespace Eigen;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])   {
    /* variables (IX, IY input) */
    const mxArray *IXA, *IYA, *IVA, *ItA; //Just names
    double *OW, *IX, *IY, *IV, *It;
    double t;
    bool isTdouble=false;
    /*For-Loops*/
    mwSize i,j,k, ItemSize, numel;
    if (nlhs > 1)
        mexErrMsgTxt("Too many output parameters");
    if (nrhs < 3)
        mexErrMsgTxt("Not enough input parameters, at least three (arrays of) matrices needed");
    else if (nrhs < 4) {
        t = 1.0;
        isTdouble = true;
    }
    else if (nrhs == 4) {
        if ( mxIsDouble(prhs[3]) ) {
            t = *mxGetPr(prhs[3]);
            isTdouble=true;
        } else {
            ItA = prhs[3];
            It = mxGetPr(ItA);
            isTdouble = false;
        }
    } else
        mexErrMsgTxt("Too many input parameters");
    /* Input Names */
    IXA = prhs[0];
    IYA = prhs[1];
    IVA = prhs[2];
    IX = mxGetPr(IXA); // Base Points
    IY = mxGetPr(IYA); // Second Points
    IV = mxGetPr(IVA); // Directions
    size_t Xs = mxGetNumberOfDimensions(IXA);
    const mwSize *Xn = mxGetDimensions(IXA);
    size_t Ys = mxGetNumberOfDimensions(IYA);
    const mwSize *Yn = mxGetDimensions(IYA);
    size_t Vs = mxGetNumberOfDimensions(IVA);
    const mwSize *Vn = mxGetDimensions(IVA);
    size_t Ts;
    const mwSize *Tn;
    if (!isTdouble) {
        Ts = mxGetNumberOfDimensions(ItA);
        Tn = mxGetDimensions(ItA);
    }
    numel=1;
    if (Xs<2)
        mexErrMsgTxt("No Matrix X given, only an array.");
    else if ( (Xs==2) && (!isTdouble) )
        mexErrMsgTxt("The coefficient t is of wrong dimensions");        
    else if (Xs>2) {
        for (i=2; i<Xs; i++) {
            numel *= Xn[i];
            if ( (!isTdouble) && (Xn[i] != Tn[i-2]))
                mexErrMsgTxt("The coefficient t is of wrong dimensions");        
        }
    }    
    /* get dimensions of the input matrix */
    ItemSize = Xn[0];
    /*
     * Checks
     */
    if (ItemSize != Xn[1])
        mexErrMsgTxt("Matrices have to be suqare.");
    if ( (Xs!=Vs) || (Xs!=Ys) )
        mexErrMsgTxt("Matrix arrays have to be of same dimension.");
    for (i=0; i<Xs; i++) {
        if ( (Vn[i]!=Xn[i]) || (Yn[i]!=Xn[i]) )
            mexErrMsgTxt("All 3 Matrices of Matrix arrays have to be of same size");
    }
    /* Output Name */
    plhs[0] = mxCreateNumericArray(Xs, Xn, mxDOUBLE_CLASS, mxREAL);
    OW = mxGetPr(plhs[0]);

    // Initialize result: V matrix array
    MatrixXd lW(ItemSize,ItemSize), lX(ItemSize,ItemSize), lV(ItemSize,ItemSize), lY(ItemSize,ItemSize);
    for (i=0; i<numel; i++) {//for all matrices
        if (!isTdouble)
            t = It[i];
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract local copies of submatrices
                lX(j,k) = IX[j + ItemSize*k + ItemSize*ItemSize*i];
                lV(j,k) = IV[j + ItemSize*k + ItemSize*ItemSize*i];
                lY(j,k) = IY[j + ItemSize*k + ItemSize*ItemSize*i];
            }
        }
        lW = mSPDParallelTransport(lX,lY,lV,t);
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract copy back to result V
                OW[j + ItemSize*k + ItemSize*ItemSize*i] = lW(j,k);
            }
        }
    }
    //Nothing to clear, Eigen clears itself, all others are input/output
}