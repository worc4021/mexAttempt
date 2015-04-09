#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include "mex.h"
#include "matrix.h"
#include "dMat.h"
#include "mexInterface.h"

#define OUTPUT_PRECISION mxDOUBLE_CLASS


struct dMat *readMXArray(const mxArray *pm){
    struct dMat *retVal;
    retVal = malloc (sizeof (*retVal));
    assert (retVal != NULL);

    double *ptr = mxGetPr(pm);

    retVal->m = mxGetM(pm);
    retVal->n = mxGetN(pm);
    retVal->data = malloc (retVal->m * retVal->n *sizeof(*retVal->data));
    assert (retVal->data != NULL);
    for (int i = 0; i != dMat_Rows(retVal); ++i)
    {
        for (int j = 0; j != dMat_Cols(retVal); ++j)
        {
            dMat_setValue(retVal, i, j, ptr[i + j*dMat_Rows(retVal)]);
        }
    }
    return retVal;
}


mxArray *writeMXArray(const struct dMat *A){
    mxArray *retVal;
    double *ptr;
    retVal = mxCreateNumericMatrix(dMat_Rows(A), dMat_Cols(A), OUTPUT_PRECISION, mxREAL);
    ptr = mxGetPr(retVal);
    for (int i = 0; i != dMat_Rows(A); ++i)
    {
        for (int j = 0; j != dMat_Cols(A); ++j)
        {
            ptr[i + j*dMat_Rows(A)] = dMat_getValue(A, i, j);
        }
    }
    return retVal;
}