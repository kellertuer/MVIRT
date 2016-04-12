/*
 * Manifold Valued Data Processing - Toolbox
 * 
 * Wrapper for the gradient in X of the second order differences of SPD matrices
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
    const mxArray *PF, *PW, *PE, *PI; //Just names
    double *OX, *F, *W;
    double E = 0;
    double I = 0;
    /*For-Loops*/
    mwSize i,j,k,l;
    /* Input Names */
    if (nrhs != 4)
        mexErrMsgTxt("Wrong number of inputs; four arguments are required.");
    if (nlhs > 1)
        mexErrMsgTxt("Too many output parameters.");
    PF = prhs[0];
    PW = prhs[1];   
    PE = prhs[2];
    PI = prhs[3];
    F = mxGetPr(PF); // Base Points
    W = mxGetPr(PW); // Weights

    
    size_t Fs = mxGetNumberOfDimensions(PF);
    const mwSize *Fn = mxGetDimensions(PF);
    size_t Xs = Fs-1;
    mwSize Xn[3];
    size_t Ws = mxGetNumberOfDimensions(PW);
    const mwSize *Wn = mxGetDimensions(PW);
    if ( mxIsDouble(PE) )
        E = *mxGetPr(PE);            
    else
        mexErrMsgTxt("The stopping critireon should be double");
    if (mxIsDouble(PI)){
        I = *mxGetPr(PI);         
        if (I<=0.0)
            mexErrMsgTxt("MaxIteration should be positive");
    }
    else
        mexErrMsgTxt("MaxIteration should be a positive number");   
    
    if (Fs < 2 && Ws<2)
        mexErrMsgTxt("No support for empty arrays as values or weights");
    ItemSize = Fn[0];
    if (ItemSize != Fn[1])
        mexErrMsgTxt("Only squared matrices can be handled");
    if (Fs > 2){
        if (Fs > 4)
            mexErrMsgTxt("Input formate Matrix x (n means) x (m variable per mean)");
        for (i=2;i<Fs;i++){
            if (Fn[i] != Wn[i-2])
                mexErrMsgTxt("Sizes of matrix array and weights do not fit");
        }            
    }
    // Is work to do?
    if (Fs <= 3){
        // No mean just copy the initial Values back
        
        plhs[0] = mxCreateNumericArray(Fs, Fn, mxDOUBLE_CLASS, mxREAL);
        OX = mxGetPr(plhs[0]);
        if (Fs == 2){
            for (j=0; j<ItemSize; j++) {
                for (k=0; k<ItemSize; k++) {
                    //Extract copy back to result V
                    OX[j + ItemSize*k] = F[j + ItemSize*k];
                }
            }
        }
        else{
            for ( i=0;i<Fn[2];i++){
                for (j=0; j<ItemSize; j++) {
                    for (k=0; k<ItemSize; k++) {
                        //Extract copy back to result V
                        OX[j + ItemSize*k + ItemSize*ItemSize*i] = F[j + ItemSize*k + ItemSize*ItemSize*i];
                    }
                }
            }        
        }
    }else{ 
        size_t Xs = Fs-1;
        mwSize Xn[3];
        for (i = 0;i<3;i++)
            Xn[i] = Fn[i];
        plhs[0] = mxCreateNumericArray(Xs, Xn, mxDOUBLE_CLASS, mxREAL);
        OX = mxGetPr(plhs[0]);    
        
        MatrixXd lX(ItemSize,ItemSize);       
        MatrixXd *lF;        
        lF = new MatrixXd[Fn[3]];
        double *lW;
        lW = new double[Fn[3]];
        for (i=0;i<Fn[3];i++){
            lF[i] = Eigen::MatrixXd::Zero(ItemSize,ItemSize);
        }
        // Calculate the means
        for (i=0;i<Fn[2];i++){
            for (j=0;j<Fn[3];j++){
                for (l = 0;l< ItemSize;l++){
                    for (k = 0;k<ItemSize;k++){
                        lF[j](l,k) = F[l+k*ItemSize+i*ItemSize*ItemSize+j*ItemSize*ItemSize*Fn[2]];                       
                    }                    
                }
                lW[j] = W[i+j*Fn[2]];
            }                      
           
            lX = mSPDMean(lF,lW,E,I,Fn[3]);
            for (j=0; j<ItemSize; j++) {
                for (k=0; k<ItemSize; k++) {
                    //Extract copy back to result X
                    OX[j + ItemSize*k + ItemSize*ItemSize*i] = lX(j,k);  
                }
            }
        }
        delete[] lF;
        delete[] lW;
        lF = NULL;
        
    } 
}
     