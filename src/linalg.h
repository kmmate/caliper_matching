/* Header for linalg.c */

#ifndef LINALG_H_
#define LINALG_H_

// includes
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "vector.h"
#include "matrix.h"



// functions
int linalg_dger(double alpha, const vector *x, const vector *y, matrix *m);
int linalg_cholesky_decomp(matrix *m);
int linalg_cholesky_invert(matrix *cholesky);




#endif // LINGALG_H_