/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the exponential function of Sn
 * 
 * ---
 * Manifold-valued Image Restoration Toolbox 1.0
 * R. Bergmann ~ 2015-05-22
 */
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
    const mxArray *IPA, *IVA, *ItA; //Just names
    double *IP, *IV, *It, *OQ;
    double t=0;
    bool isTdouble=false;
    /* variables */
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
    IPA = prhs[0];
    IVA = prhs[1];
    IP = mxGetPr(IPA); // Base Points
    IV = mxGetPr(IVA); // Second Points
    size_t Ps = mxGetNumberOfDimensions(IPA);
    const mwSize *Pn = mxGetDimensions(IPA);
    size_t Vs = mxGetNumberOfDimensions(IVA);
    const mwSize *Vn = mxGetDimensions(IVA);
    size_t ts = 0;
    const mwSize *tn;
    if (!isTdouble) {
        ts = mxGetNumberOfDimensions(ItA);
        tn = mxGetDimensions(ItA);
    }
    mwSize numel=1,i=0,j=0,ItemSize = Pn[0];
    if ((Ps==1) && (!isTdouble))
        mexErrMsgIdAndTxt("SnExp:WrongSizeOfT","The coefficient t is of wrong dimension."); 
    else if (Ps>1) {
        if ((!isTdouble) && (ts!=(Ps-1)))
            mexErrMsgIdAndTxt("SnExp:WrongSizeOfT","The input t is of wrong size.");
        for (i=1; i<Ps; i++) {
            numel *= Pn[i];
            if ( (!isTdouble) && (Pn[i] != tn[i-1]))
                mexErrMsgTxt("The coefficient t is of wrong dimensions");
        }
    }
    /* get dimensions of the input vectors */
    /*
     * Checks
     */
    if (Ps!=Vs)
        mexErrMsgTxt("Vector P and V arrays have to be of same dimension.");
    for (i=0; i<Ps; i++) {
        if (Pn[i]!=Vn[i])
            mexErrMsgTxt("Vector arrays are not of same size.");
    }
    /* Create Output */
    plhs[0] = mxCreateNumericArray(Ps,Pn, mxDOUBLE_CLASS, mxREAL);
    OQ = mxGetPr(plhs[0]);
    //Compute result
    double nv=0;
    VectorXd X = VectorXd::Zero(ItemSize), V = VectorXd::Zero(ItemSize), Y = VectorXd::Zero(ItemSize);
    for (i=0; i<numel; i++) {
        if (!isTdouble)
            t = It[i];
        for (j=0; j<ItemSize; j++) {
            X(j) = IP[j+ItemSize*i];
            V(j) = IV[j+ItemSize*i];
        }
        Y = mSnExp(X,V,t);
        for (j=0; j<ItemSize; j++) {
            OQ[j + ItemSize*i] = Y(j);
        }
    }
}