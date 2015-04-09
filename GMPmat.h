#pragma once
#include "GMPtypes.h"

void GMPmat_getRow(mpz_t *ropN, mpz_t *ropD, struct GMPmat *A, size_t r);
void GMPmat_destroy (struct GMPmat *A);
void GMPmat_print(const struct GMPmat *A);
struct GMPmat *GMPmat_create (size_t m, size_t n, int init);
size_t GMPmat_Rows (const struct GMPmat *A);
size_t GMPmat_Cols (const struct GMPmat *A);
void GMPmat_setValue (const struct GMPmat *A, size_t r, size_t c, double val);
void GMPmat_getValue (mpq_t rop, const struct GMPmat *A, size_t r, size_t c);
static void __GMPmat_validate_indices (const struct GMPmat *A, size_t r, size_t c);
void GMPmat_printRow(const struct GMPmat *A, size_t r);
struct GMPmat *GMPmat_dropCols(struct GMPmat *A, size_t d);
void GMPmal_everyNrows(struct GMPmat *A, size_t N, char *type);
void GMPmat_invertSignForFacetEnumeration(struct GMPmat *A);