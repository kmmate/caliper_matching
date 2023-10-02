////////////////////////////////////////////////////////////////////////////
//
// * FILE: matrix.c
// * DESCRIPTION:
//    Matrices. Curently mostly a wrapper around GSL types and functions.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 04/oct/2022
//
//////////////////////////////////////////////////////////////////////////

#include "matrix.h"


//////////////////      matrix

matrix *matrix_alloc(size_t n1, size_t n2){
    matrix *m = gsl_matrix_alloc(n1, n2);
    return m;
}


matrix *matrix_calloc(size_t n1, size_t n2){
    matrix *m = gsl_matrix_calloc(n1, n2);
    return m;
}


void matrix_free(matrix *m){
    gsl_matrix_free(m);
}

double matrix_get(matrix *m, size_t i, size_t j){
    return gsl_matrix_get(m, i, j);
}


// vector_view of `i`th row of matrix `m`. `&view` is a vector pointer.
vector_view matrix_row(matrix *m, size_t i){
    vector_view view = gsl_matrix_row(m, i);
    return view;
}

// vector_view of `j`th column of matrix `m`. `&view.vector` is a vector pointer.
vector_view matrix_column(matrix *m, size_t j){
    vector_view view = gsl_matrix_column(m, j);
    return view;
}

void matrix_set(matrix *m, size_t i, size_t j, double x){
    gsl_matrix_set(m, i, j, x);
}

// Sets a_{i,j} to a_{i,j}+b_{i,j} for i=0,...,a->size1-1, j=0,...,a->size2-1.
int matrix_add(matrix *a, matrix *b){
    return gsl_matrix_add(a, b);
}

int matrix_isequal(matrix *a, matrix *b){
    if ((a->size1 != b->size1) || (a->size2 != b->size2)){  // not equal if they're of different sizes
        return 0;
    } else {
        int equal = 1;
        for (int i=0; i<a->size1; i++){
            for (int j=0; j<a->size2; j++){
                equal *= fabs(matrix_get(a, i, j) - matrix_get(b, i, j)) <= pow(10, -MATRIX_TEST_DOUBLE_TOLERANCE_NPOWER);
            }
        }
        return equal;
    }
}

void test_matrix_isequal(void){
    // test 1: different row number
    matrix *a1 = matrix_alloc(3, 2);
    matrix *b1 = matrix_alloc(4, 2);
    assert(matrix_isequal(a1, b1) == 0);
    matrix_free(a1); matrix_free(b1);
    
    // test 2: different column number
    matrix *a2 = matrix_alloc(3, 2);
    matrix *b2 = matrix_alloc(3, 5);
    assert(matrix_isequal(a2, b2) == 0);
    matrix_free(a2); matrix_free(b2);
    
    // test 3: different row and column number
    matrix *a3 = matrix_alloc(4, 2);
    matrix *b3 = matrix_alloc(3, 5);
    assert(matrix_isequal(a3, b3) == 0);
    matrix_free(a3); matrix_free(b3);

    // test 4: matching dimensions, different values
    matrix *a = matrix_alloc(3, 2);
    matrix *b = matrix_alloc(3, 2);
    matrix_set(a, 0, 0, 10); matrix_set(b, 0, 0, -10);
    matrix_set(a, 0, 1, 10); matrix_set(b, 0, 1, 10);
    matrix_set(a, 1, 0, 10); matrix_set(b, 1, 0, 10);
    matrix_set(a, 1, 1, 10); matrix_set(b, 1, 1, 10);
    matrix_set(a, 2, 0, 10); matrix_set(b, 2, 0, 10);
    matrix_set(a, 2, 1, 10); matrix_set(b, 2, 1, 10);
    assert(matrix_isequal(a, b) == 0);

    // test 4: matching dimensions, same values
    matrix_set(a, 0, 0, 10); matrix_set(b, 0, 0, 10);
    matrix_set(a, 0, 1, 10); matrix_set(b, 0, 1, 10);
    matrix_set(a, 1, 0, 10); matrix_set(b, 1, 0, 10);
    matrix_set(a, 1, 1, 10); matrix_set(b, 1, 1, 10);
    matrix_set(a, 2, 0, 10); matrix_set(b, 2, 0, 10);
    matrix_set(a, 2, 1, 10); matrix_set(b, 2, 1, 10);
    assert(matrix_isequal(a, b) == 1);

    printf("Tests for `matrix_isequal` completed: no issues found.\n");
}



//////////////////      spmatrix

// spmatrix_short *spmatrix_short_alloc(size_t n1, size_t n2){
//     spmatrix_short *m = gsl_spmatrix_short_alloc(n1, n2);
//     return m;
// };

// spmatrix_short *spmatrix_short_alloc_nzmax(size_t n1, size_t n2, size_t nzmax, size_t sptype){
//     spmatrix_short *m = gsl_spmatrix_short_alloc_nzmax(n1, n2, nzmax, sptype);
//     return m;
// }

// void spmatrix_short_free(spmatrix_short *m){
//     gsl_spmatrix_short_free(m);
// };

// int spmatrix_short_equal(spmatrix_short *a, spmatrix_short *b){
//     return gsl_spmatrix_short_equal(a, b);
// }

// void spmatrix_short_set_zero(spmatrix_short *m){
//     gsl_spmatrix_short_set_zero(m);
// }

// short spmatrix_short_get(spmatrix_short *m, size_t i, size_t j){
//     return gsl_spmatrix_short_get(m, i, j);
// };

// void spmatrix_short_set(spmatrix_short *m, size_t i, size_t j, short x){
//     gsl_spmatrix_short_set(m, i, j, x);
// };

// // ERROR: gsl_spmatrix_short_add does not work correctly, do not use this!!!
// void spmatrix_short_add(spmatrix_short *c, spmatrix_short *a, spmatrix_short *b){
//     gsl_spmatrix_short_add(c, a, b);
// }

// void test_spmatrix_short_add(void){
//     // test 1
//     spmatrix_short *c = spmatrix_short_alloc_nzmax(4, 4, 4, GSL_SPMATRIX_CRS);
//     spmatrix_short *a = spmatrix_short_alloc_nzmax(4, 4, 4, GSL_SPMATRIX_COO);
//     spmatrix_short *b = spmatrix_short_alloc_nzmax(4, 4, 4, GSL_SPMATRIX_COO);
//     spmatrix_short_set(a, 0, 0, 2); spmatrix_short_set(a, 0, 1, 10);
//     spmatrix_short_set(b, 0, 0, 5); 
//     spmatrix_short_set(b, 1, 1, 17);
//     spmatrix_short *res = spmatrix_short_alloc_nzmax(4, 4, 4, GSL_SPMATRIX_COO);
//     spmatrix_short_set(res, 0, 0, 5); spmatrix_short_set(res, 0, 1, 0);
//     spmatrix_short_set(res, 1, 1, 17);
//     spmatrix_short_add(c, c, spmatrix_short_compress(b, GSL_SPMATRIX_CRS));
//     assert(spmatrix_short_equal_robust(c, res));
//     // test 2
//     spmatrix_short_set(res, 0, 0, 7); spmatrix_short_set(res, 0, 1, 10);
//     spmatrix_short_set(res, 1, 1, 17);
//     spmatrix_short_add(c, spmatrix_short_compress(a, GSL_SPMATRIX_CRS), spmatrix_short_compress(b, GSL_SPMATRIX_CRS));
//     assert(spmatrix_short_equal_robust(c, res));
//     // test 3: repeated addition to self. This test fails!!!
//     spmatrix_short *c3 = spmatrix_short_alloc_nzmax(4, 4, 4, GSL_SPMATRIX_CSR);
//     spmatrix_short *b3 = spmatrix_short_alloc_nzmax(4, 4, 4, GSL_SPMATRIX_COO);
//     spmatrix_short *res3 = spmatrix_short_alloc_nzmax(4, 4, 4, GSL_SPMATRIX_COO);
//     spmatrix_short_set(b3, 0, 0, 5); spmatrix_short_set(b3, 0, 2, 20); 
//     spmatrix_short_set(b3, 1, 1, 17);
//     spmatrix_short_set(b3, 3, 1, 40);
//     spmatrix_short_add(c3, c3, spmatrix_short_compress(b3, GSL_SPMATRIX_CSR));
//     spmatrix_short_set(res3, 0, 0, 1*5); spmatrix_short_set(res3, 0, 2, 1*20);
//     spmatrix_short_set(res3, 1, 1, 1*17);
//     spmatrix_short_set(res3, 3, 1, 1*40);
//     for (int i=0; i<4; i++){
//         for(int j=0; j<4; j++){
//             printf("\nb3(%d,%d)=%hd", i, j, spmatrix_short_get(b3, i, j));
//             printf("\nc3(%d,%d)=%hd", i, j, spmatrix_short_get(c3, i, j));
//         }
//     }
//     assert(spmatrix_short_equal_robust(c3, res3));
//     printf("First assert passed.\n");
//     spmatrix_short_add(c3, c3, spmatrix_short_compress(b3, GSL_SPMATRIX_CRS));
//     spmatrix_short_set(res3, 0, 0, 2*5); spmatrix_short_set(res3, 0, 2, 2*20);
//     spmatrix_short_set(res3, 1, 1, 2*17);
//     spmatrix_short_set(res3, 3, 1, 2*40);
//     assert(spmatrix_short_equal_robust(c3, res3));
//     printf("Tests for `spmatrix_short_add` completed: no issues found.\n");
// }

// /*
// Elementwise addition of sparse matrices `a` and `b` of  with result loaded into COO-type sparse matrix `c`.
// */
// void spmatrix_short_add_robust(spmatrix_short *c, spmatrix_short *a, spmatrix_short *b){
//     assert((a->size1 == b->size1) && (b->size1 == c->size1) && (a->size2 == b->size2) && (b->size2 == c->size2));
//     for (int i=0; i<a->size1; i++){
//         for (int j=0; j<a->size2; j++){
//             spmatrix_short_set(c, i, j, spmatrix_short_get(a, i, j) + spmatrix_short_get(b, i, j));
//         }
//     }
// }

// void spmatrix_short_csr(spmatrix_short *dest, spmatrix_short *src){
//     gsl_spmatrix_short_csr(dest, src);
// }

// void spmatrix_short_csc(spmatrix_short *dest, spmatrix_short *src){
//     gsl_spmatrix_short_csc(dest, src);
// }

// spmatrix_short *spmatrix_short_compress(spmatrix_short *src, size_t sptype){
//     return gsl_spmatrix_short_compress(src, sptype);
// }

// // Comparison of two spmatrix_short, robust to sptype.
// int spmatrix_short_equal_robust(spmatrix_short *a, spmatrix_short *b){
//     int equal = 1;
//     assert((a->size1 == b->size1) && (a->size2 == b->size2));
//     for (int i=0; i<a->size1; i++){
//         for (int j=0; j<a->size2; j++){
//             equal *= spmatrix_short_get(a, i, j) == spmatrix_short_get(b, i, j);
//         }
//     }
//     return equal;
// }

void test_matrix(void){
    test_matrix_isequal();
    // test_spmatrix_short_add();
    printf("=== Tests: all tests for `matrix` completed: no issues found.\n");
}