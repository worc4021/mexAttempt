#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <gmp.h>
#include "tmwtypes.h"
#include "GMPmat.h"

static void __GMPmat_validate_indices (const struct GMPmat *A, size_t r, size_t c)
     {
       assert (A != NULL);
       assert (r < A->m);
       assert (c < A->n);
     }

void GMPmat_getValue (mpq_t rop, const struct GMPmat *A, size_t r, size_t c)
{
       __GMPmat_validate_indices (A, r, c);
       mpq_set (rop, A->data[r*(A->n) + c]);
}

void GMPmat_setValue (const struct GMPmat *A, size_t r, size_t c, double val)
{
       __GMPmat_validate_indices (A, r, c);
       mpq_set_d(A->data[r*(A->n) + c], val);
}

size_t GMPmat_Cols (const struct GMPmat *A)
     {
       assert (A != NULL);

       return (A->n);
     }


size_t GMPmat_Rows (const struct GMPmat *A)
     {
       assert (A != NULL);

       return (A->m);
     }


struct GMPmat *GMPmat_create (size_t m, size_t n, int init)
     {
       struct GMPmat *A;
       size_t         i, j;

       A = calloc (1, sizeof (*A));
       assert (A != NULL);

       A->m    = m;
       A->n    = n;
       A->data = calloc (m*n, sizeof(*A->data));
       assert (A->data != NULL );

       if (init != 0)
       {
         for (i=0; i!=m; ++i){
           for (j=0; j!=n; ++j){
            mpq_init (A->data[i*(A->n) + j]);
            
          }
        }
       }

       return (A);
}

void GMPmat_print(const struct GMPmat *A)
{
    const char *val;

    val = getenv ("SILENCE");
    if (val != NULL)
        return;

    assert( A!= NULL );
    fprintf(stdout, "\n");
    for (size_t i = 0; i < GMPmat_Rows(A); ++i)
    {
        for (size_t j = 0; j < GMPmat_Cols(A); ++j)
        {
            mpq_out_str(stdout, 10, A->data[i*(A->n) + j]);
            fprintf(stdout, " ");
        }
        fprintf(stdout, "\n" );
    }
    fprintf(stdout, "\n" );
}

void GMPmat_destroy (struct GMPmat *A)
     {
       assert (A != NULL);
       for (size_t i = 0; i < GMPmat_Rows(A); ++i)
       {
        for (size_t j = 0; j < GMPmat_Cols(A); ++j)
        {
            mpq_clear(A->data[i*(A->n) + j]);
        }
       }
       if (A->data != NULL) free (A->data);
       A->data = NULL;
       A->n   = 0;
       A->m   = 0;
       free (A);
     }

void GMPmat_getRow(mpz_t *ropN, mpz_t *ropD, struct GMPmat *A, size_t r)
{
    assert( r < A->m );
    assert( ropN != NULL && ropD != NULL );
    for (size_t i = 0; i < GMPmat_Cols(A); ++i)
    {
        mpz_set(ropN[i], mpq_numref(A->data[r*(A->n) + i]));
        mpz_set(ropD[i], mpq_denref(A->data[r*(A->n) + i]));
    }
}

void GMPmat_printRow(const struct GMPmat *A, size_t r)
{
  assert ( r < A->m );
  size_t i, n;
  n = GMPmat_Cols(A);
  mpq_t curVal;
  mpq_init(curVal);
  fprintf(stdout, "\n" );
  for ( i = 0; i < n ; ++i){
    GMPmat_getValue(curVal, A, r, i);
    mpq_out_str(stdout, 10, curVal);
    fprintf(stdout, " ");
  }
  fprintf(stdout, "\n" );
  mpq_init(curVal);
}

struct GMPmat *GMPmat_dropCols(struct GMPmat *A, size_t d)
{
  assert ( A != NULL );
  assert ( d < GMPmat_Cols(A) );
  size_t m, n, i, j;
  m = GMPmat_Rows(A);
  n = GMPmat_Cols(A) - d;

  struct GMPmat *retVal;
  retVal = GMPmat_create( m, n, 0 );

  for (i=0; i!=m; ++i){
    for (j=0; j!=n; ++j){
      mpq_init ( retVal->data[i*n + j] );
      mpq_set ( retVal->data[i*n + j], A->data[i*(A->n) + j] );
    }
  }

  GMPmat_destroy(A);
  return retVal;
}

void GMPmal_everyNrows(struct GMPmat *A, size_t N, char *type)
{
  assert ( A != NULL );
  if ( GMPmat_Rows(A) % N == 0 )
  {
    fprintf(stdout, "So far %zu %s found.\n", (GMPmat_Rows(A)/N)*N, type);
  }
}

void GMPmat_invertSignForFacetEnumeration(struct GMPmat *A)
{
  assert ( A != NULL );
  size_t i, j, m, n;
  m = GMPmat_Rows(A);
  n = GMPmat_Cols(A);
  mpq_t holder;
  mpq_init(holder);
  for (i = 0; i < m; i++)
  {
    for (j = 1; j < n; j++)
    {
      GMPmat_getValue(holder, A, i, j);
      mpq_neg(A->data[i*(A->n) + j], holder);
    }
  }
  mpq_clear(holder);
}