#include "GMPtypes.h"

struct GMPmat *dMat2GMPmat(struct dMat *A)
{
	struct GMPmat *retVal;
	retVal = GMPmat_create(dMat_Rows(A),dMat_Cols(A),0);
	for (int i = 0; i < dMat_Rows(A); ++i)
	{
		for (int j = 0; j < dMat_Cols(A); ++j)
		{
			mpq_init(retVal->data[i*dMat_Cols(A)+j]);
			mpq_set_d(retVal->data[i*dMat_Cols(A)+j], dMat_getValue(A, i, j));
		}
	}
	return retVal;
}

struct dMat *GMPmat2dMat(struct GMPmat *A)
{
	struct dMat *retVal;
	retVal = dMat_create(GMPmat_Rows(A),GMPmat_Cols(A),0);
	for (int i = 0; i < GMPmat_Rows(A); ++i)
	{
		for (int j = 0; j < GMPmat_Cols(A); ++j)
		{
			dMat_setValue(retVal, i, j, mpq_get_d(A->data[i*GMPmat_Cols(A)+j]));
		}
	}
	return retVal;
}