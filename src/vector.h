/* Header for vector.c */

#ifndef VECTOR_H_
#define VECTOR_H_

// includes
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_short.h>


#define VECTOR_TEST_DOUBLE_TOLERANCE_NPOWER 9


 // structs
typedef gsl_vector vector;  // vector is type alias for gsl_vector
typedef gsl_vector_short vector_short;
typedef gsl_vector_int vector_int;

// functions
// doubles
vector *vector_alloc(size_t n);
vector *vector_calloc(size_t n);
void vector_set(vector *v, size_t i, double x);
double vector_get(vector *v, size_t i);
double *vector_ptr(vector *v, size_t i);
void vector_free(vector *v);
int vector_is_zero_one(vector *v);
double vector_sum(vector *v);
int vector_add(vector *a, vector *b);
int vector_sub(vector *a, vector *b);
int vector_isequal(vector *a, vector *b);

// shorts
vector_short *vector_short_alloc(size_t n);
vector_short *vector_short_calloc(size_t n);
void vector_short_set(vector_short *v, size_t i, short x);
short vector_short_get(vector_short *v, size_t i);
short *vector_short_ptr(vector_short *v, size_t i);
void vector_short_free(vector_short *v);
int vector_short_is_zero_one(vector_short *v);
int vector_short_sum(vector_short *v);
double vector_short_mean(vector_short *v);


// ints
vector_int *vector_int_alloc(size_t n);
vector_int *vector_int_calloc(size_t n);
void vector_int_set(vector_int *v, size_t i, int x);
int vector_int_get(vector_int *v, size_t i);
int *vector_int_ptr(vector_int *v, size_t i);
void vector_int_free(vector_int *v);
int vector_int_equal(vector_int *v1, vector_int *v2);
int vector_int_is_zero_one(vector_int *v);
int vector_int_sum(vector_int *v);
int vector_int_add(vector_int *a, vector_int *b);

// tests
void test_vector(void);


#endif // VECTOR_H_