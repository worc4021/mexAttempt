#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include "GMPmat.h"
#include "dMat.h"
#include "translation.h"

#define d_time(t) ((double)(t.tv_sec)+(double)(t.tv_usec)/1000000)

void mpz_row_clean(mpz_t *row, size_t m)
{
    assert(row != NULL);
    for (size_t i = 0; i < m; ++i)
    {
        mpz_clear(row[i]);
    }
    free(row);
}

void mpq_row_clean(mpq_t *row, size_t m)
{
    assert(row != NULL);
    for (size_t i = 0; i < m; ++i)
    {
        mpq_clear(row[i]);
    }
    free(row);
}

void mpz_row_init(mpz_t *row, size_t m)
{
    assert(row != NULL);
    for (size_t i = 0; i < m; ++i)
    {
        mpz_init(row[i]);
    }
}

void mpq_row_init(mpq_t *row, size_t m)
{
    assert(row != NULL);
    for (size_t i = 0; i < m; ++i)
    {
        mpq_init(row[i]);
    }
}

void mpz_norm(mpz_t norm, mpz_t *row, size_t m)
{
    assert( row != NULL );
    mpz_t help;
    mpz_init(help);
    for (size_t i = 0; i < m; ++i)
    {
        mpz_addmul(help, row[i], row[i]);
    }
    mpz_sqrt(norm, help);
    mpz_clear(help);
}

struct GMPmat *GMPmat_appendRow(struct GMPmat *A, mpq_t *row)
{
    assert( A != NULL && row != NULL );
    struct GMPmat *retVal;
    size_t m,n, i, j;
    m = GMPmat_Rows(A) + 1;
    n = GMPmat_Cols(A);
    retVal = GMPmat_create( m, n, 0);
    assert( m != 0 );

    for (i = 0; i < (m-1); ++i)
    {
        for (j = 0; j < n; ++j)
        {
            mpq_init(retVal->data[i*n + j]);
            mpq_set(retVal->data[i*n + j],A->data[i*n + j]);
        }
    }
    for (i = 0; i < n; ++i)
    {
        mpq_init(retVal->data[(m-1)*n + i]);
        mpq_set(retVal->data[(m-1)*n + i], row[i]);
    }
    GMPmat_destroy(A);
    return retVal;
}

void mpz_to_mpq(mpq_t *rop, mpz_t *op, size_t m)
{
    size_t i;
    assert( rop != NULL && op != NULL);
    if (!mpz_sgn(op[0]))
    {
        mpz_t norm;
        mpz_init(norm);
        mpz_norm(norm, op, m);
        for (i = 0; i < m; ++i)
        {
            mpq_set_num(rop[i], op[i]);
            mpq_set_den(rop[i], norm);
            mpq_canonicalize(rop[i]);
        }
        mpz_clear(norm);
    }else{
        for (i = 0; i < m; ++i)
        {
            mpq_set_num(rop[i], op[i]);
            mpq_set_den(rop[i], op[0]);
            mpq_canonicalize(rop[i]);
        }
    }
}

struct GMPmat *dMat2GMPmat(struct dMat *A)
{
    struct GMPmat *retVal;
    retVal = GMPmat_create(dMat_Rows(A),dMat_Cols(A),0);
    for (size_t i = 0; i < dMat_Rows(A); ++i)
    {
        for (size_t j = 0; j < dMat_Cols(A); ++j)
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
    for (size_t i = 0; i < GMPmat_Rows(A); ++i)
    {
        for (size_t j = 0; j < GMPmat_Cols(A); ++j)
        {
            dMat_setValue(retVal, i, j, mpq_get_d(A->data[i*GMPmat_Cols(A)+j]));
        }
    }
    return retVal;
}

void mpz_row_print(mpz_t *row, size_t n)
{
    assert ( row != NULL );
    size_t i;
    for (i = 0; i < n; ++i)
    {
        mpz_out_str(stdout, 10, row[i]);
        fprintf(stdout, " " );
    }
    fprintf(stdout, "\n" );
}

void mpq_row_print(mpq_t *row, size_t n)
{
    assert ( row != NULL );
    size_t i;
    for (i = 0; i < n; ++i)
    {
        mpq_out_str(stdout, 10, row[i]);
        fprintf(stdout, " " );
    }
    fprintf(stdout, "\n" );
}



void mpz_print(mpz_t op)
{
    assert( op != NULL );
    fprintf(stdout, "\n");
    mpz_out_str(stdout, 10, op);
    fprintf(stdout, "\n");
}

void mpq_print(mpq_t op)
{
    assert( op != NULL );
    fprintf(stdout, "\n");
    mpq_out_str(stdout, 10, op);
    fprintf(stdout, "\n");
}

void mpz_print_product(mpz_t numA, mpz_t denA, mpz_t numB, mpz_t denB)
{
    assert( numA != NULL && denA != NULL && numB != NULL && denB != NULL );
    mpq_t a,b, res;
    mpq_init(a);
    mpq_init(b);
    mpq_init(res);
    
    mpq_set_num(a, numA);
    mpq_set_den(a, denA);
    mpq_canonicalize(a);
    
    mpq_set_num(b, numB);
    mpq_set_den(b, denB);
    mpq_canonicalize(b);
    
    mpq_mul(res, a, b);

    mpq_print(res);
    
    mpq_clear(res);
    mpq_clear(b);
    mpq_clear(a);
}
mpq_t *mpq_row_extract(const struct GMPmat *A, size_t r)
{
    assert ( A != NULL && r < GMPmat_Rows(A) );
    mpq_t *retVal;
    size_t i, n; 
    n = GMPmat_Cols(A);
    retVal = malloc( n*sizeof(*retVal) );
    for ( i = 0; i < n; ++i )
    {
        mpq_init( retVal[i] );
        GMPmat_getValue( retVal[i], A, r, i);
    }
    return retVal;

}

// void timeIt(char *proName)
// {
//     long maxRsc;
//     struct rusage rusage;
//     getrusage (RUSAGE_SELF, &rusage);
//     maxRsc = rusage.ru_maxrss;
//     if (maxRsc / 1024 < 10 ){
//     fprintf(stderr, "%s took %0.3f seconds for computation and %0.3f and used a maximum of %ldKb.\n", 
//         proName, 
//         d_time (rusage.ru_utime),
//         d_time (rusage.ru_stime),
//         maxRsc);
//     } else {
//         fprintf(stderr, "%s took %0.3f seconds for computation and %0.3f and used a maximum of %ld.%ldMb.\n", 
//         proName, 
//         d_time (rusage.ru_utime),
//         d_time (rusage.ru_stime),
//         maxRsc/1024,
//         maxRsc%1024);

//     }
// }


#ifdef NOTMATLAB
    void clockout(clock_t in, char *type)
    {
        clock_t out = clock();
        double cpu_time_used = ((double) (out - in)) / CLOCKS_PER_SEC;
        if (cpu_time_used < 60)
        {
            fprintf(stdout, "%s took %f seconds.\n", type, cpu_time_used);
        } else if (cpu_time_used > 60 && cpu_time_used < 3600)
        {
            fprintf(stdout, "%s took %d minutes and %d seconds.\n", type, (int) floor(cpu_time_used/60),
                (int)(cpu_time_used - floor(cpu_time_used/60)*60 ) );
        } else
        {
            fprintf(stdout, "%s took %d hours, %d minutes and %d seconds.\n", type, 
                (int) floor(cpu_time_used/3600),
                (int) (floor((cpu_time_used - floor(cpu_time_used/3600)*3600)/60) ),
                (int) floor(cpu_time_used - (floor((cpu_time_used - floor(cpu_time_used/3600)*3600)/60)*60 ) )
                );
        }
    }
#endif /* NOTMATLAB */

void lprint(long *array, size_t l)
{
    fprintf(stdout, "\n");
    for (size_t i = 0; i < l; i++)
        fprintf(stdout, "%ld ", array[i]);
    fprintf(stdout, "\n");
}