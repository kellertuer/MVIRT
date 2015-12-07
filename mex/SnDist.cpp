#include "mex.h"
#include "math.h"
#include <stdlib.h>
#include <algorithm>
#include "manifolds/manifoldSn.h"

using namespace std;

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])   {
    /* Input */
    const mxArray *IPA, *IQA; //Just names
    double *IP, *IQ, *Od;
    /* variables */
    if (nrhs != 2)
        mexErrMsgTxt("Wrong number of inputs; two arrays of vectors required.");
    IPA = prhs[0];
    IQA = prhs[1];
    IP = mxGetPr(IPA); // Base Points
    IQ = mxGetPr(IQA); // Second Points
    size_t Ps = mxGetNumberOfDimensions(IPA);
    const mwSize *Pn = mxGetDimensions(IPA);
    size_t Qs = mxGetNumberOfDimensions(IQA);
    const mwSize *Qn = mxGetDimensions(IQA);
    mwSize numel=1,i=0,j=0,ItemSize = Pn[0];
    if (Ps>1) {
        for (i=1; i<Ps; i++) {
            numel *= Pn[i];
        }
    }    
    /* get dimensions of the input vectors */
    /*
     * Checks
     */
    if (Ps!=Qs)
        mexErrMsgTxt("Vector arrays have to be of same dimension.");
    for (i=0; i<Ps; i++) {
        if (Pn[i]!=Qn[i])
            mexErrMsgTxt("Vector arrays are not of same size.");
    }
    /* Create Output */
    if (Ps==1) {
        mwSize dims[1];
        dims[0]=1;
        plhs[0] = mxCreateNumericArray(1,dims, mxDOUBLE_CLASS, mxREAL);
        Od = mxGetPr(plhs[0]);
    }
    else if (Ps>1) {
        mwSize* dims = (mwSize*)mxCalloc(Ps-1,sizeof(mwSize));
       /* mwSize dims[Ps-1];*/
        for (i=1; i<Ps; i++) {
            dims[i-1] = Pn[i];
        }
        plhs[0] = mxCreateNumericArray(Ps-1,dims, mxDOUBLE_CLASS, mxREAL);
        Od = mxGetPr(plhs[0]);
        mxFree(dims);
    }
    VectorXd X=VectorXd::Zero(ItemSize), Y=VectorXd::Zero(ItemSize);
    for (i=0; i<numel; i++) {
        for (j=0; j<ItemSize; j++) {
            X(j) = IP[j + ItemSize*i];
            Y(j) = IQ[j + ItemSize*i];
        }       
        Od[i] = mSnDist(X,Y);
    }
}