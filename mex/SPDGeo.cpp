/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the geodesic point function of SPD matrices
 * 
 * ---
 * Manifold-valued Image Restoration Toolbox 1.0
 * J. Persch ~ 2017-03-31
 */
#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <Eigen/Dense>
#include "manifolds/manifoldSPD.h"

using namespace std;
using namespace Eigen;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])   {
    /* variables (IX, IY input) */
    const mxArray *IXA, *IYA, *ItA; //Just names
    mxArray *OVA;
    double *OV, *IX, *IY, *It;
    double t=0;
    bool isTdouble=false;
    /*For-Loops*/
    mwSize i,j,k,N,ItemSize;
    /* Input Names */
    if (nrhs != 3)
        mexErrMsgTxt("Wrong number of inputs; two arrays of matrices required and a array of scalars.");
    if ( mxIsScalar(prhs[2]) ) {
        t = *mxGetPr(prhs[2]);
        isTdouble=true;
    } else {
        ItA = prhs[2];
        It = mxGetPr(ItA);
        isTdouble = false;
    }
    IXA = prhs[0];
    IYA = prhs[1];
    IX = mxGetPr(IXA); // Base Points
    IY = mxGetPr(IYA); // Second Points
    size_t Xs = mxGetNumberOfDimensions(IXA);
    const mwSize *Xn = mxGetDimensions(IXA);
    size_t Ys = mxGetNumberOfDimensions(IYA);
    const mwSize *Yn = mxGetDimensions(IYA);
    size_t Ts=0;
    const mwSize *Tn;
    if (!isTdouble) {
        Ts = mxGetNumberOfDimensions(ItA);
        Tn = mxGetDimensions(ItA);
        if (Ts == 2 & Tn[1] == 1) 
            Ts -= 1;            
    }
    mwSize numel=1;
    if (Xs<2)
        mexErrMsgTxt("No Matrix X given, only an vector.");
    else if ( (Xs==2) && (!isTdouble) )
        mexErrMsgTxt("The coefficient t is of wrong dimension."); 
    else if( (Xs>2) && (!isTdouble) ) {
        if (Ys!=Xs)
            mexErrMsgTxt("The input parameter Y is of wrong size.");
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
    OV = mxGetPr(plhs[0]);
    
    /* get dimensions of the input matrix */
    ItemSize = Xn[0];
    /*
     * Checks
     */
    if (ItemSize != Xn[1])
        mexErrMsgTxt("Matrices have to be suqare.");
    if (Xs!=Ys)
        mexErrMsgTxt("Matrix arrays have to be of same dimension.");
    for (i=0; i<Xs; i++) {
        if (Yn[i]!=Xn[i])
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
                lY(j,k) = IY[j + ItemSize*k + ItemSize*ItemSize*i];
            }
        }
        // Would it be better to generate the complete array and pass that? We have to for-loop here anyways...
        lV = mSPDGeo(lX,lY,t);
        for (j=0; j<ItemSize; j++) {
            for (k=0; k<ItemSize; k++) {
                //Extract copy back to result V
                OV[j + ItemSize*k + ItemSize*ItemSize*i] = lV(j,k);
            }
        }
    }
    // Eigen clears itself, all others are in/out or not pointers
    // -> nothing to clear
}