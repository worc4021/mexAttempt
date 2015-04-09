#pragma once
#include "GMPtypes.h"


void dMat_print(struct dMat *A);
struct dMat *dMat_create (size_t m, size_t n, int init);
size_t dMat_Rows (const struct dMat *A);
size_t dMat_Cols (const struct dMat *A);
void dMat_setValue (const struct dMat *A, size_t r, size_t c, double val);
double dMat_getValue (const struct dMat *A, size_t r, size_t c);
static void __dMat_validate_indices (const struct dMat *A, size_t r, size_t c);
#ifdef NOTMATLAB
	struct dMat *dMatFromFile(int *dim);
	void toFile(struct dMat *A);
#endif /* NOTMATLAB */
void dMat_destroy (struct dMat *A);
struct dMat *dMat_horCon(struct dMat *A, struct dMat *B);
struct dMat *dMat_negate(struct dMat *A);
