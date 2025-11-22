/* 
mexDisAggSol.cpp
[x0] = mexDisAggSol(x1,K0s,Ks)
*/
// #include <string> 
#include <iostream>
#include "mex.h"

using namespace std;

/* gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{    
    /* [A,c,Ks] = mexAggregateSDPcone(A0,c0,K0s,Ks) */
    /* ************************************************** */
    
        /* Check for proper number of input and output arguments */    
    if (nrhs != 3) {
        mexErrMsgTxt("3 input arguments required.");
    } 
    if(nlhs != 1){
        mexErrMsgTxt("1 output argument required");
    } 

    /* ---> */ 
       
    /* reading the dense column vector x1 */
    const mxArray* x1_ptr = prhs[0];
    if (mxIsSparse(x1_ptr) == 1) {
        mexErrMsgTxt("x1 must be a dense column vector.");
    }
    mwSize  x1colSize = mxGetN(x1_ptr);
    if (x1colSize > 1) {
        mexErrMsgTxt("x1 must be a column vector.");
    }     
    // mwSize  x1Size = mxGetM(x1_ptr);
    double* x1pr = mxGetPr(x1_ptr);
//    std::cout << "x1pr = ";    
//    for (unsigned int j=0; j < x1Size; j++) 
//        std::cout << x1pr[j] << " "; 
//    std::cout << std::endl;
    
    /* reading the dense mex vector K0s */ 
    const mxArray* K0s_ptr = prhs[1];
    if (mxIsSparse(K0s_ptr) == 1) {
        mexErrMsgTxt("K0s must be a dense column vector.");
    }
    mwSize  K0colSize = mxGetN(K0s_ptr);
    if (K0colSize > 1) {
        mexErrMsgTxt("K0s must be a column vector.");
    }     
    mwSize  K0sSize = mxGetM(K0s_ptr);
    double* K0spr = mxGetPr(K0s_ptr);
//    std::cout << "K0spr = ";
//    for (unsigned int k0=0; k0 < K0sSize; k0++) 
//        std::cout << K0spr[k0] << " "; 
//    std::cout << std::endl;
    
    /* reading the dense mex vector K1s */ 
    const mxArray* K1s_ptr = prhs[2];
    if (mxIsSparse(K1s_ptr) == 1) {
        mexErrMsgTxt("K1s must be a dense column vector.");
    }
    mwSize  K1colSize = mxGetN(K1s_ptr);
    if (K1colSize > 1) {
        mexErrMsgTxt("K1s must be a column vector.");
    }     
    // mwSize  K1sSize = mxGetM(K1s_ptr);
    double* K1spr = mxGetPr(K1s_ptr);
//    std::cout << "K1spr = ";
//    for (unsigned int k1=0; k1 < K1sSize; k1++) 
//        std::cout << K1spr[k1] << " "; 
//    std::cout << std::endl;
    
    /* Constructing the dense column vector x0 */
    
    /*   the column size of x0 be constructed */
    mwSize x0Size = 0;
    for (unsigned int k0 = 0; k0 < K0sSize; k0++)
        x0Size += (mwSize)(K0spr[k0]*K0spr[k0]);
//    std::cout << "x0Size = " << x0Size << std::endl;

    plhs[0] = mxCreateDoubleMatrix(x0Size,1,mxREAL);
    double*  x0pr = mxGetPr(plhs[0]);
    
    /* Construction of x0 */
    
    int k1 = -1; 
    int topNDblks = 0; 
    int bottomNDblks = 0; 
    unsigned int c0colPt = 0;
    unsigned int c1colPt = 0; 
    for (unsigned int k0 = 0; k0 < K0sSize; k0++) {
        bottomNDblks -= (int) K0spr[k0];
        if (bottomNDblks < 0) {
            k1++; 
            topNDblks = 0;
            bottomNDblks = (int) (K1spr[k1] - K0spr[k0]); 
        } 
        else {
            topNDblks += (int) K0spr[k0-1]; 
        }             
//        std::cout << " k0 = " << k0 << " topNDblks = " << topNDblks << " bottomNDblks =  " << bottomNDblks << std::endl;
        for (unsigned int j0 = 0; j0 < K0spr[k0]; j0++) {
            c1colPt += topNDblks;             
            for (unsigned int i0 = 0; i0 < K0spr[k0]; i0++) {
                x0pr[c0colPt+i0] = x1pr[c1colPt+i0]; 
            }
            c0colPt += (int) K0spr[k0]; 
            c1colPt += (int) (K0spr[k0]+bottomNDblks); 
        }
    }
    
//    for (unsigned int j=0; j < x0Size; j++) 
//        std::cout << x0pr[j] << " "; 
//    std::cout << std::endl;

    return;
}
