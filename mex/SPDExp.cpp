/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the exponential function of SPD matrices
 * 
 * ---
 * Manifold-valued Image Restoration Toolbox 1.0
 * R. Bergmann ~ 2015-04-12
 */
#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <Eigen/Dense>
#include "manifolds/manifoldSPD.h"

using namespace std;
using namespace Eigen;

size_t ItemSize;  /* size of matrix */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])   {
    /* variables (IX, IY input) */
    const mxArray *IXA, *IVA, *ItA; //Just names
    mxArray *OYA;
    double *OY, *IX, *IV, *It;
    double t=0;
    bool isTdouble=false;
    /*For-Loops*/
    mwSize i,j,k;
    if (nlhs > 1)
        mexErrMsgTxt("Too many output parameters");
    if (nrhs < 2)
        mexErrMsgTxt("Not enough input parameters, at least two (arrays of) matrices needed");
    else if (nrhs < 3) {
        t = 1.0;
        isTdouble = true;
    }
    else if (nrhs == 3) {
        if ( mxIsDouble(prhs[2]) ) {
            t = *mxGetPr(prhs[2]);
            isTdouble=true;
        } else {
            ItA = prhs[2];
            It = mxGetPr(ItA);
            isTdouble = false;
        }
    } else
        mexErrMsgTxt("Too many input parameters");
    /* Input Names */
    IXA = prhs[0];
    IVA = prhs[1];
    IX = mxGetPr(IXA); // Base Points
    IV = mxGetPr(IVA); // Second Points
    size_t Xs = mxGetNumberOfDimensions(IXA);
    const mwSize *Xn = mxGetDimensions(IXA);
    size_t Vs = mxGetNumberOfDimensions(IVA);
    const mwSize *Vn = mxGetDimensions(IVA);
    size_t Ts=0;
    const mwSize *Tn;
    if (!isTdouble) {
        Ts = mxGetNumberOfDimensions(ItA);
        Tn = mxGetDimensions(ItA);
    }
    mwSize numel=1;
    if (Xs<2)
        mexErrMsgTxt("No Matrix X given, only an array.");
    else if ( (Xs==2) && (!isTdouble) )
        mexErrMsgTxt("The coefficient t is of wrong dimension.");        
    else if( (Xs>2) && (!isTdouble) ) {
        if (Vs!=Xs)
            mexErrMsgTxt("The input parameter V is of wrong size.");
        if (Ts!=(Xs-2))
            mexErrMsgTxt("The input parameter t is of wrong size.");
        for (i=2; i<Xs; i++) {
            numel *= Xn[i];
            if ( Xn[i] != Tn[i-2] )
                mexErrMsgTxt("The coefficient t is of wrong dimensions");        
        }
    }
    else if (isTdouble) {
        for (i=2; i<Xs; i++) 
            numel *= Xn[i];             
    }
    /* Output Name */
    plhs[0] = mxCreateNumericArray(Xs, Xn, mxDOUBLE_CLASS, mxREAL);
    OY = mxGetPr(plhs[0]);
    /* get dimensions of the input matrix */
    ItemSize = Xn[0];
    /*
     * Checks
     */
    if (ItemSize != Xn[1])
        mexErrMsgTxt("Matrices have to be square.");
    if (Xs!=Vs)
        mexErrMsgTxt("Matrix arrays have to be of same dimension.");
    for (i=0; i<Xs; i++) {
        if (Vn[i]!=Xn[i])
            mexErrMsgTxt("Matrices of Matrix arrays are not of same size.");
    }
    // Initialize result: V matrix array
    MatrixXd lX(ItemSize,ItemSize), lY(ItemSize,ItemSize), lV(ItemSize,ItemSize);
    for (i=0; i<numel; i++) {//for all matrices
        if (!isTdouble)
            t = It[i];
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract local copies of submatrices
                lX(j,k) = IX[j + ItemSize*k + ItemSize*ItemSize*i];
                lV(j,k) = IV[j + ItemSize*k + ItemSize*ItemSize*i];
            }
        }
        lY = mSPDExp(lX,lV,t);
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract copy back to result V
                OY[j + ItemSize*k + ItemSize*ItemSize*i] = lY(j,k);
            }
        }
    }
    //Nothing to clear, Eigen clears itself, all others are input/output
}
