#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mex.h"
#include "matrix.h"
#include "tmwtypes.h"

#include "GMPtypes.h"
#include "GMPmat.h"
#include "dMat.h"
#include "mexInterface.h"
#include "translation.h"
#include "projection.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if ( nrhs != 3 ) {
        mexErrMsgIdAndTxt("MATLAB:mexTest:rhs","Input has to be A, B, dim!");
    } 

    struct dMat *A, *B, *retMat;
    struct GMPmat *calcMat;
    int dim;

    A = readMXArray(prhs[0]);
    B = readMXArray(prhs[1]);
    dim = (int) mxGetScalar(prhs[2]);

    A = dMat_negate(A);
    
    retMat = dMat_horCon(B, A);

    calcMat = dMat2GMPmat(retMat);
    dMat_destroy(retMat);
    calcMat = H2V(calcMat);
    retMat = GMPmat2dMat(calcMat);
    GMPmat_destroy(calcMat);

    plhs[0] = writeMXArray(retMat);
    dMat_destroy(retMat);   
}