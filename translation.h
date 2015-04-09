#pragma once
#include "GMPtypes.h"
#include <time.h>

struct dMat *GMPmat2dMat(struct GMPmat *A);
struct GMPmat *dMat2GMPmat(struct dMat *A);
void mpz_to_mpq(mpq_t *rop, mpz_t *op, size_t m);
struct GMPmat *GMPmat_appendRow(struct GMPmat *A, mpq_t *row);
void mpz_norm(mpz_t norm, mpz_t *row, size_t m);
void mpq_row_init(mpq_t *row, size_t m);
void mpz_row_init(mpz_t *row, size_t m);
void mpq_row_clean(mpq_t *row, size_t m);
void mpz_row_clean(mpz_t *row, size_t m);
void mpz_row_print(mpz_t *row, size_t n);
void mpq_row_print(mpq_t *row, size_t n);
void mpq_print(mpq_t op);
void mpz_print(mpz_t op);
void mpz_print_product(mpz_t numA, mpz_t denA, mpz_t numB, mpz_t denB);
mpq_t *mpq_row_extract(const struct GMPmat *A, size_t r);
void timeIt(char *proName);
#ifdef NOTMATLAB
    void clockout(clock_t in, char *type);
#endif /* NOTMATLAB */
void lprint(long *array, size_t l);