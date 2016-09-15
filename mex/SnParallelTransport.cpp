/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the C++-Function to compute the parallel transport on 
 * symmetric posivive definte matrices.
 *
 * ---
 * ManImRes, R. Bergmann ~ 2015-04-12
 */
#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "manifolds/manifoldSn.h"

using namespace std;
using namespace Eigen;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])   {
    /* variables (IX, IY input) */
    const mxArray *IXA, *IYA, *IVA; //Just names
    double *OW, *IX, *IY, *IV;
    /*For-Loops*/
    mwSize i,j,k, ItemSize, numel;
    if (nlhs > 1)
        mexErrMsgTxt("Too many output parameters");
    if (nrhs < 3)
        mexErrMsgTxt("Not enough input parameters, at least three (arrays of) matrices needed");
    else if (nrhs > 3) 
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
    numel=1;
    if (Xs<1)
        mexErrMsgTxt("No Vector X given, only a number.");
    else if (Xs>1) {
        for (i=1; i<Xs; i++) {
            numel *= Xn[i];             
        }
    }    
    /* get dimensions of the input matrix */
    ItemSize = Xn[0];
    /*
     * Checks
     */
    if ( (Xs!=Vs) || (Xs!=Ys) )
        mexErrMsgTxt("Vectors have to be of same dimension.");
    for (i=0; i<Xs; i++) {
        if ( (Vn[i]!=Xn[i]) || (Yn[i]!=Xn[i]) )
            mexErrMsgTxt("All 3 Vector arrays have to be of same size");
    }
    /* Output Name */
    plhs[0] = mxCreateNumericArray(Xs, Xn, mxDOUBLE_CLASS, mxREAL);
    OW = mxGetPr(plhs[0]);

    // Initialize result: V matrix array
    VectorXd lW(ItemSize), lX(ItemSize), lV(ItemSize), lY(ItemSize);
    for (i=0; i<numel; i++) {//for all matrices       
        for (j=0; j<ItemSize; j++) {
                //Extract local copies of submatrices
                lX(j) = IX[j + ItemSize*i];
                lV(j) = IV[j + ItemSize*i];
                lY(j) = IY[j + ItemSize*i];
        }
        lW = mSnParallelTransport(lX,lY,lV);
        for (j=0; j<ItemSize; j++) {
                //Extract copy back to result V
                OW[j + ItemSize*i] = lW(j);
        }
    }
    //Nothing to clear, Eigen clears itself, all others are input/output
}
