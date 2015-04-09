#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <gmp.h>
#include "lrslib.h"
#include "GMPmat.h"
#include "translation.h"
#include "projection.h"

size_t pN = 20;

int my_lrs_init()
{
    assert ( lrs_mp_init (ZERO, stdin, stdout) );
    lrs_global_count = 0;
    lrs_checkpoint_seconds = 0;
#ifdef SIGNALS
    setup_signals();
#endif
    return 0;
}

// void lrs_print_data(lrs_dic *P)
// {
//     assert ( P != NULL );
//     long m, n;
//     mpq_t frac;
//     mpq_init(frac);
//     fprintf(stdout, "\n" );
//     for (m = 0; m < P->m; ++m)
//     {
//         for (n = 0; n < P->d_orig; ++n)
//         {
//             mpq_set_num(frac, P->A[m][n]);
//             mpq_set_den(frac, P->det);
//             mpq_canonicalize(frac);
//             mpq_out_str(stdout, 10, frac);
//             fprintf(stdout, " ");            
//         }
//         fprintf(stdout, "\n");
//     }
//     fprintf(stdout, "\n");
//     mpq_clear(frac);
// }

struct GMPmat *projection(struct GMPmat *inp, int d)
{

    lrs_dic *Pv, *Ph; /* structure for holding current dictionary and indices  */
    lrs_dat *Qv, *Qh; /* structure for holding static problem data               */
    lrs_mp_vector output; /* one line of output:ray,vertex,facet,linearity */
    lrs_mp_matrix Lin;    /* holds input linearities if any are found      */

      size_t i, j, cols, rows;
      long col;     /* output column index for dictionary            */

      /* Global initialization - done once */

      assert( my_lrs_init () == 0 );

      Qv = lrs_alloc_dat ("LRS globals");
      assert( Qv!= NULL );

      Qv->m = GMPmat_Rows(inp);
      Qv->n = GMPmat_Cols(inp);

      output = lrs_alloc_mp_vector (Qv->n);

      lrs_mp_vector num, den;
      num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
      den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

      Pv = lrs_alloc_dic (Qv);   /* allocate and initialize lrs_dic      */
      assert( Pv != NULL );

      
      struct GMPmat *Helper;
      Helper = GMPmat_create(0, GMPmat_Cols(inp), 1);
      
      mpq_t *curRow;
      curRow = calloc(GMPmat_Cols(inp), sizeof(mpq_t));
      assert( curRow != NULL );
      mpq_row_init(curRow, GMPmat_Cols(inp));



      for (i = 1; i <= GMPmat_Rows(inp); ++i)
      {
        GMPmat_getRow(num, den, inp, i-1);
        lrs_set_row_mp(Pv,Qv,i,num,den,GE);
      }

      assert( lrs_getfirstbasis (&Pv, Qv, &Lin, TRUE) );


      for (col = 0L; col < Qv->nredundcol; col++)  /* print linearity space */
        lrs_printoutput (Qv, Lin[col]); 

      do
        {
          for (col = 0L; col <= Pv->d; col++)
            if (lrs_getsolution (Pv, Qv, output, col)) {
              mpz_to_mpq(curRow, output, GMPmat_Cols(Helper));
              Helper = GMPmat_appendRow(Helper, curRow);
              GMPmal_everyNrows(Helper, pN, "vertices/rays");
            }
        }
        while (lrs_getnextbasis (&Pv, Qv, FALSE));

        mpq_row_clean(curRow, GMPmat_Cols(Helper));
        lrs_clear_mp_vector (output, Qv->n);
        lrs_clear_mp_vector (num, Qv->n);
        lrs_clear_mp_vector (den, Qv->n);
        lrs_free_dic (Pv,Qv);       /* deallocate lrs_dic */
        lrs_free_dat (Qv);          /* deallocate lrs_dat */

        Helper = reducevertices(Helper);

        Qh = lrs_alloc_dat ("LRS globals");
        assert( Qh != NULL );

        Qh->m = GMPmat_Rows(Helper);
        Qh->n = GMPmat_Cols(Helper) - d;

        Qh->hull = TRUE;     /* convex hull problem: facet enumeration      */
        Qh->polytope = TRUE;  /* input is a polytope                         */

        output = lrs_alloc_mp_vector (Qh->n);
        num = lrs_alloc_mp_vector (Qh->n);
        den = lrs_alloc_mp_vector (Qh->n);

        Ph = lrs_alloc_dic (Qh);
        assert( Ph != NULL );
        
        struct GMPmat *retVal;
        retVal = GMPmat_create(0, Qh->n, 0);
        
        rows = GMPmat_Rows (Helper);
        cols = GMPmat_Cols (retVal);

        curRow = calloc(cols, sizeof(mpq_t));
        assert( curRow != NULL );
        mpq_row_init(curRow, cols);

        mpq_t curVal;
        mpq_init(curVal);

       for (i = 0; i < rows; ++i)
       {
         for (j = 0; j < cols; ++j)
          {
            GMPmat_getValue (curVal, Helper, i, j);
            mpz_set (num[j], mpq_numref(curVal));
            mpz_set (den[j], mpq_denref(curVal));
          }
          lrs_set_row_mp (Ph, Qh, i+1 ,num, den, GE);
        }

        mpq_clear(curVal);

        assert( lrs_getfirstbasis (&Ph, Qh, &Lin, TRUE) );

        for (col = 0L; col < Qh->nredundcol; col++)  /* print linearity space */
          lrs_printoutput (Qh, Lin[col]);

        do
        {
          for (col = 0L; col <= Ph->d; col++)
            if (lrs_getsolution (Ph, Qh, output, col)){
              mpz_to_mpq(curRow, output, GMPmat_Cols(retVal));
              retVal = GMPmat_appendRow(retVal, curRow);
              GMPmal_everyNrows(retVal, pN, "inequalities");
            }
        }
        while (lrs_getnextbasis (&Ph, Qh, FALSE));

        GMPmat_destroy(Helper);
        mpq_row_clean(curRow, GMPmat_Cols(retVal));
        lrs_clear_mp_vector (output, Qh->n);
        lrs_clear_mp_vector (num, Qh->n);
        lrs_clear_mp_vector (den, Qh->n);
        lrs_free_dic (Ph,Qh);
        lrs_free_dat (Qh);

        // lrs_close ("lrsTrial:");
        printf ("\n");

  GMPmat_destroy(inp);
  return retVal;
}

struct GMPmat *H2V(struct GMPmat *inp)
{

    lrs_dic *Pv, *Ph; /* structure for holding current dictionary and indices  */
    lrs_dat *Qv, *Qh; /* structure for holding static problem data               */
    lrs_mp_vector output; /* one line of output:ray,vertex,facet,linearity */
    lrs_mp_matrix Lin;    /* holds input linearities if any are found      */

      size_t i, j, cols, rows;
      long col;     /* output column index for dictionary            */

      /* Global initialization - done once */

      assert( my_lrs_init () == 0 );

      Qv = lrs_alloc_dat ("LRS globals");
      assert( Qv!= NULL );

      Qv->m = GMPmat_Rows(inp);
      Qv->n = GMPmat_Cols(inp);
      // Qv->printslack = 1L;
      // Qv->printcobasis = 1L;

      output = lrs_alloc_mp_vector (Qv->n);

      lrs_mp_vector num, den;
      num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
      den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

      Pv = lrs_alloc_dic (Qv);   /* allocate and initialize lrs_dic      */
      assert( Pv != NULL );

      
      struct GMPmat *Helper;
      Helper = GMPmat_create(0, GMPmat_Cols(inp), 1);
      
      mpq_t *curRow;
      curRow = calloc(GMPmat_Cols(inp), sizeof(mpq_t));
      assert( curRow != NULL );
      mpq_row_init(curRow, GMPmat_Cols(inp));



      for (i = 1; i <= GMPmat_Rows(inp); ++i)
      {
        GMPmat_getRow(num, den, inp, i-1);
        lrs_set_row_mp(Pv,Qv,i,num,den,GE);
      }

      assert( lrs_getfirstbasis (&Pv, Qv, &Lin, TRUE) );


      for (col = 0L; col < Qv->nredundcol; col++)  /* print linearity space */
        lrs_printoutput (Qv, Lin[col]); 

      do
        {
          for (col = 0L; col <= Pv->d; col++)
            if (lrs_getsolution (Pv, Qv, output, col)) {
              mpz_to_mpq(curRow, output, GMPmat_Cols(Helper));
              Helper = GMPmat_appendRow(Helper, curRow);
              GMPmal_everyNrows(Helper, pN, "vertices/rays");
              // GMPmat_printRow(Helper, Helper->m -1);
              // print_basis (stdout, Qv);
            }
        }
        while (lrs_getnextbasis (&Pv, Qv, FALSE));

        mpq_row_clean(curRow, GMPmat_Cols(Helper));
        lrs_clear_mp_vector (output, Qv->n);
        lrs_clear_mp_vector (num, Qv->n);
        lrs_clear_mp_vector (den, Qv->n);
        lrs_free_dic (Pv,Qv);       /* deallocate lrs_dic */
        lrs_free_dat (Qv);          /* deallocate lrs_dat */

        // lrs_close ("lrsTrial:");
        GMPmat_destroy(inp);
        return Helper;
}

struct GMPmat *V2H(struct GMPmat *inp) /* This function is untested */
{
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    long i;
    long col;

    assert( my_lrs_init () == 0 );

    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = GMPmat_Rows(inp);
    Q->n = GMPmat_Cols(inp);
    Q->hull = TRUE;
    Q->polytope = TRUE;

    output = lrs_alloc_mp_vector (Q->n);

    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
    den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

    P = lrs_alloc_dic (Q);
    assert ( P != NULL );

    struct GMPmat *retMat;
    retMat = GMPmat_create(0, GMPmat_Cols(inp), 1);
      
    mpq_t *curRow;
    curRow = calloc(GMPmat_Cols(inp), sizeof(mpq_t));
    assert( curRow != NULL );
    mpq_row_init(curRow, GMPmat_Cols(inp));


    for (i = 1; i <= GMPmat_Rows(inp); ++i)
    {
      GMPmat_getRow(num, den, inp, i-1);
      lrs_set_row_mp(P ,Q ,i ,num ,den , GE);
    }

    assert ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) );
    
    for (col = 0L; col < Q->nredundcol; col++)
      lrs_printoutput (Q, Lin[col]);

    do
    {
      for (col = 0; col <= P->d; col++)
        if (lrs_getsolution (P, Q, output, col)){
          mpz_to_mpq(curRow, output, GMPmat_Cols(retMat));
          retMat = GMPmat_appendRow(retMat, curRow);
          GMPmal_everyNrows(retMat, pN, "inequalities");
        }
    }
    while (lrs_getnextbasis (&P, Q, FALSE));

    mpq_row_clean( curRow, GMPmat_Cols(retMat) );
    lrs_clear_mp_vector ( output, Q->n);
    lrs_clear_mp_vector ( num, Q->n);
    lrs_clear_mp_vector ( den, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);
    return retMat;

}

struct GMPmat *reducemat(struct GMPmat *inp)
{
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    long i;
    long col;

    size_t m = GMPmat_Rows(inp);

    size_t *redRows;
    redRows = malloc( m*sizeof(*redRows) );
    assert( redRows != NULL );

    assert( my_lrs_init () == 0 );

    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = m;
    Q->n = GMPmat_Cols(inp);

    output = lrs_alloc_mp_vector (Q->n);

    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
    den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

    P = lrs_alloc_dic (Q);
    assert ( P != NULL );

    struct GMPmat *retMat;
    retMat = GMPmat_create(0, GMPmat_Cols(inp), 1);

    for (i = 1; i <= m; ++i)
    {
      GMPmat_getRow(num, den, inp, i-1);
      lrs_set_row_mp(P ,Q ,i ,num ,den , GE);
    }

    assert ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) );

    size_t lastdv = Q->lastdv;
    size_t d = P->d;
    size_t ineq;
    m = P->m_A;

    for (i = lastdv + 1; i <= m + d; ++i)
    {
      ineq = Q->inequality[i - lastdv] - 1;
      if (!checkindex(P, Q, i))
      {
        retMat = GMPmat_appendRow(retMat, mpq_row_extract(inp, ineq));
        GMPmal_everyNrows(retMat, pN, "irredundant inequalities");
      }
    }

    lrs_clear_mp_vector ( output, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);
    return retMat;
  }

struct GMPmat *reducevertices(struct GMPmat *inp)
{
    lrs_dic *P;
    lrs_dat *Q;
    lrs_mp_vector output;
    lrs_mp_matrix Lin;

    long i;
    long col;

    size_t m = GMPmat_Rows(inp);

    size_t *redRows;
    redRows = malloc( m*sizeof(*redRows) );
    assert( redRows != NULL );

    assert( my_lrs_init () == 0 );

    Q = lrs_alloc_dat ("LRS globals");
    assert ( Q != NULL );
    Q->m = m;
    Q->n = GMPmat_Cols(inp);
    Q->hull = TRUE;
    Q->polytope = TRUE;

    output = lrs_alloc_mp_vector (Q->n);

    lrs_mp_vector num, den;
    num = lrs_alloc_mp_vector(GMPmat_Cols(inp));
    den = lrs_alloc_mp_vector(GMPmat_Cols(inp));

    P = lrs_alloc_dic (Q);
    assert ( P != NULL );

    struct GMPmat *retMat;
    retMat = GMPmat_create(0, GMPmat_Cols(inp), 1);

    for (i = 1; i <= m; ++i)
    {
      GMPmat_getRow(num, den, inp, i-1);
      lrs_set_row_mp(P ,Q ,i ,num ,den , GE);
    }

    assert ( lrs_getfirstbasis (&P, Q, &Lin, TRUE) );

    size_t lastdv = Q->lastdv;
    size_t d = P->d;
    size_t ineq;
    m = P->m_A;

    for (i = lastdv + 1; i <= m + d; ++i)
    {
      ineq = Q->inequality[i - lastdv] - 1;
      if (!checkindex(P, Q, i))
      {
        retMat = GMPmat_appendRow(retMat, mpq_row_extract(inp, ineq));
        GMPmal_everyNrows(retMat, pN, "irredundant vertices/rays");
      }
    }

    lrs_clear_mp_vector ( output, Q->n);
    lrs_free_dic ( P , Q);
    lrs_free_dat ( Q );

    GMPmat_destroy(inp);
    return retMat;
  }
