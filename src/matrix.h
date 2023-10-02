/* Header for matrix.c */

#ifndef MATRIX_H_
#define MATRIX_H_

// includes
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_blas.h>

#define MATRIX_TEST_DOUBLE_TOLERANCE_NPOWER 9

 // structs
typedef gsl_matrix matrix;  // matrix is type alias for gsl_matrix
typedef gsl_vector_view vector_view;  // vector_view os type alies for gsl_vector_view
// typedef gsl_spmatrix_short spmatrix_short;

// functions
// matrix
matrix *matrix_alloc(size_t n1, size_t n2);
matrix *matrix_calloc(size_t n1, size_t n2);
void matrix_free(matrix *m);
double matrix_get(matrix *m, size_t i, size_t j);
vector_view matrix_row(matrix *m, size_t i);
vector_view matrix_column(matrix *m, size_t j);
void matrix_set(matrix *m, size_t i, size_t j, double x);
int matrix_add(matrix *a, matrix *b);
int matrix_isequal(matrix *a, matrix *b);
// spmatrix_short
// spmatrix_short *spmatrix_short_alloc(size_t n1, size_t n2);
// spmatrix_short *spmatrix_short_alloc_nzmax(size_t n1, size_t n2, size_t nzmax, size_t sptype);
// void spmatrix_short_free(spmatrix_short *m);
// int spmatrix_short_equal(spmatrix_short *a, spmatrix_short *b);
// short spmatrix_short_get(spmatrix_short *m, size_t i, size_t j);
// void spmatrix_short_set(spmatrix_short *m, size_t i, size_t j, short x);
// void spmatrix_short_set_zero(spmatrix_short *m);
// void spmatrix_short_add(spmatrix_short *c, spmatrix_short *a, spmatrix_short *b);
// void spmatrix_short_add_robust(spmatrix_short *c, spmatrix_short *a, spmatrix_short *b);
// void spmatrix_short_csr(spmatrix_short *dest, spmatrix_short *src);
// void spmatrix_short_csc(spmatrix_short *dest, spmatrix_short *src);
// spmatrix_short *spmatrix_short_compress(spmatrix_short *src, size_t sptype);
// int spmatrix_short_equal_robust(spmatrix_short *a, spmatrix_short *b);


// tests
void test_matrix(void);


#endif // MATRIX_H_