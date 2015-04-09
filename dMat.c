#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <gmp.h>
#include "tmwtypes.h"
#include "dMat.h"


#ifdef NOTMATLAB 
  extern char *path;
#endif /* NOTMATLAB */

struct dMat *dMat_negate(struct dMat *A)
{
    assert( A != NULL );
    size_t i, j, m, n;
    double curVal;
    m = dMat_Rows(A);
    n = dMat_Cols(A);
    struct dMat *retVal;
    retVal = dMat_create(m, n, 0);
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            curVal = dMat_getValue(A, i, j);
            dMat_setValue(retVal, i, j, -1*curVal);
        }
    }
    dMat_destroy(A);
    return retVal;
}

struct dMat *dMat_horCon(struct dMat *A, struct dMat *B)
{
    assert( A != NULL && B != NULL );
    size_t mA, nA, mB, nB, i;
    mA = dMat_Rows(A);
    nA = dMat_Cols(A);
    mB = dMat_Rows(B);
    nB = dMat_Cols(B);
    assert ( mA == nA );

    struct dMat *retVal;
    retVal = dMat_create( mA, nA + nB, 0);
    
    for (i = 0; i < mA; ++i)
    {
        memcpy(retVal->data + (nA + nB)*i, A->data + nA*i, nA*sizeof(*A->data));
        memcpy(retVal->data + (nA + nB)*i + nA, B->data + nB*i, nB*sizeof(*B->data));
    }
    
    dMat_destroy(A);
    dMat_destroy(B);

    return retVal;
}

#ifdef NOTMATLAB
  struct dMat *dMatFromFile(int *dim)
  {
    FILE *fp;
    char helper[100];
    strcpy(helper, path);
    strcat(helper, "MATLABtoLOWLEVEL");
    fp = fopen(helper ,"r");
    assert( fp != NULL );
    uint32_T out, m, n;
    assert( fread(&m ,1 , sizeof(uint32_T),fp) == sizeof(uint32_T) );
    assert( fread(&n ,1 , sizeof(uint32_T),fp) == sizeof(uint32_T) );
    struct dMat *retVal;
    retVal = dMat_create((size_t)m, (size_t)n, 0);
    for (size_t i = 0; i != dMat_Cols(retVal); ++i)
    {
      for (size_t j = 0; j != dMat_Rows(retVal); ++j)
      {
        assert( fread(retVal->data + i + j*n,1 , sizeof(double),fp) == sizeof(double) );
      }
    }
    if( fread(&out ,1 , sizeof(uint32_T),fp) == sizeof(uint32_T) ){
      *dim = (int)out;
    } else {
      *dim = 0;
    }
    fclose(fp);
    return retVal;
  }

  void toFile(struct dMat *A)
  {
    FILE *fp;
    char helper[100];
    strcpy(helper, path);
    strcat(helper, "LOWLEVELtoMATLAB");
    fp = fopen(helper, "w");
    assert( fp != NULL );
    uint32_T m = (uint32_T)dMat_Rows(A);
    uint32_T n = (uint32_T)dMat_Cols(A);
    assert( fwrite(&m , 1, sizeof(uint32_T), fp) == sizeof(uint32_T) );
    assert( fwrite(&n , 1, sizeof(uint32_T), fp) == sizeof(uint32_T) );
    for (int i = 0; i != dMat_Cols(A); ++i)
    {
      for (int j = 0; j != dMat_Rows(A); ++j)
      {
        assert( fwrite(A->data + i + j*n,1 , sizeof(double), fp) == sizeof(double) );
      }
    }
    fclose(fp);
  }
#endif /* NOTMATLAB */

static void __dMat_validate_indices (const struct dMat *A, size_t r, size_t c)
     {
       assert (A != NULL);
       assert (r < A->m);
       assert (c < A->n);
     }

double dMat_getValue (const struct dMat *A, size_t r, size_t c)
{
       __dMat_validate_indices (A, r, c);

       return (A->data[r*(A->n) + c]);
}

void dMat_setValue (const struct dMat *A, size_t r, size_t c, double val)
{
       __dMat_validate_indices (A, r, c);

       A->data[r*(A->n) + c] = val;
}

size_t dMat_Cols (const struct dMat *A)
     {
       assert (A != NULL);

       return (A->n);
     }


size_t dMat_Rows (const struct dMat *A)
     {
       assert (A != NULL);

       return (A->m);
     }


struct dMat *dMat_create (size_t m, size_t n, int init)
     {
       struct dMat *A;
       size_t         i, j;

       A = malloc (sizeof (*A));
       assert (A != NULL);

       A->m    = m;
       A->n    = n;
       A->data = malloc (m*n*sizeof(*A->data));
       assert (A->data != NULL);

       if (init != 0)
       {
         for (i=0; i!=m; ++i){
           for (j=0; j!=n; ++j){
            dMat_setValue(A, i, j, 0);
          }
        }
       }

       return (A);
}

void dMat_print(struct dMat *A)
{
    const char *val;

    val = getenv ("SILENCE");
    if (val != NULL)
        return;


    assert (A != NULL);  
    size_t m, n;
      m = dMat_Rows(A);
      n = dMat_Cols(A);
      printf("\n");
      for (size_t i = 0; i < m; ++i)        
      {
        for (size_t j = 0; j < n; ++j)
        {
          printf("%9.2e, ", dMat_getValue(A, i, j));
        }
        printf("\n");
      }
      printf("\n");
}

void dMat_destroy (struct dMat *A)
     {
       free (A->data);
       A->data = NULL;
       A->n   = 0;
       A->m   = 0;
       free (A);
     }
