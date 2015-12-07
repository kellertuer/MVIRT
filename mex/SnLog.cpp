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
    const mxArray *IPA, *IQA; //Just names
    double *IP, *IQ, *OV;
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
        plhs[0] = mxCreateNumericArray(Ps,Pn, mxDOUBLE_CLASS, mxREAL);
        OV = mxGetPr(plhs[0]);
    //Compute result
    VectorXd X=VectorXd::Zero(ItemSize), Y=VectorXd::Zero(ItemSize), V = VectorXd::Zero(ItemSize);
    double scp,normpr=0;
    for (i=0; i<numel; i++) {
        scp=0;
        for (j=0; j<ItemSize; j++) {
            X(j) = IP[j + ItemSize*i];
            Y(j) = IQ[j + ItemSize*i];
        }
        V = mSnLog(X,Y);
        for (j=0; j<ItemSize; j++)
            OV[j+ItemSize*i] = V(j);
    }
}