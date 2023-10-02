////////////////////////////////////////////////////////////////////////////
//
// * FILE: linalg.c
// * DESCRIPTION:
//    Linear algebra. Wrapper around GSL BLAS.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 14/nov/2022
//
//////////////////////////////////////////////////////////////////////////

#include "linalg.h"

// Updates matrix m = m + alpha * x * y^T
int linalg_dger(double alpha, const vector *x, const vector *y, matrix *m){
    return gsl_blas_dger(alpha, x, y, m);
}

int linalg_cholesky_decomp(matrix *m){
    return gsl_linalg_cholesky_decomp1(m);
}

int linalg_cholesky_invert(matrix *cholesky){
    return gsl_linalg_cholesky_invert(cholesky);
}