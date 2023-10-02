//////////////////////////////////////////////////////////////////////////
//
// * FILE: vector.c
// * DESCRIPTION:
//    Vectors.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 04/oct/2022
//
//////////////////////////////////////////////////////////////////////////

#include "vector.h"


////////////////////////  double

vector *vector_alloc(size_t n){
    vector *v = gsl_vector_alloc(n);
    return v;
}

vector *vector_calloc(size_t n){
    vector *v = gsl_vector_calloc(n);
    return v;
}

void vector_free(vector *v){
    gsl_vector_free(v);
}

void vector_set(vector *v, size_t i, double x){
    gsl_vector_set(v, i, x);
}

double vector_get(vector *v, size_t i){
    return gsl_vector_get(v, i);
}

// Pointer to `i`th element of a vector `v`.
double *vector_ptr(vector *v, size_t i){
    return gsl_vector_ptr(v, i);
}

// Checks if a vector only contains 0-1 entries
int vector_is_zero_one(vector *v){
    int zero_one_counter = 0;
    for (size_t i=0; i<v->size; i++){
        if ((vector_get(v, i) == 0) || (vector_get(v, i) == 1)){
            zero_one_counter++;
        }
    }
    return (zero_one_counter == v->size);
}

void test_vector_is_zero_one(void){
    vector *v = vector_alloc(3);
    // test 1: 0-1 vector
    vector_set(v, 0, 0); vector_set(v, 1, 1); vector_set(v, 2, 0);
    assert(vector_is_zero_one(v) == 1);
    // test 2: not 0-1 vector
    vector_set(v, 0, 0); vector_set(v, 1, 1); vector_set(v, 2, 3);
    assert(vector_is_zero_one(v) == 0);
    vector_free(v);
    printf("Tests for `vector_is_zero_one` completed: no issues found.\n");
}

// Sums elements of a vector
double vector_sum(vector *v){
    double sum = 0.0;
    for (size_t i=0; i<v->size; i++){
        sum += vector_get(v, i);
    }
    return sum;
}

void test_vector_sum(void){
    vector *v = vector_alloc(2);
    vector_set(v, 0, 10); vector_set(v, 1, 3);
    assert(vector_sum(v) == 13.0);
    vector_set(v, 0, 10); vector_set(v, 1, -3);
    assert(vector_sum(v) == 7.0);
    vector_free(v);
    printf("Tests for `vector_sum` completed: no issues found.\n");
}

// Tests where all elements of vector `a` and `b` are the same within given error.
int vector_isequal(vector *a, vector *b){
    if (a->size != b->size){
        return 0;
    } else {
        int equal = 1;
        for (int i=0; i<a->size; i++){
            equal *= fabs(vector_get(a, i) - vector_get(b, i)) <= pow(10, -VECTOR_TEST_DOUBLE_TOLERANCE_NPOWER);
        }
        return equal;
    }
}

void test_vector_isequal(void){
    // test 1: different length
    vector *a1 = vector_alloc(3);
    vector *b1 = vector_alloc(4);
    assert(vector_isequal(a1, b1) == 0);
    vector_free(a1); vector_free(b1);

    // test 2: same length different values
    vector *a2 = vector_alloc(3);
    vector *b2 = vector_alloc(3);
    vector_set(a2, 0, 10);   vector_set(b2, 0, -10);
    vector_set(a2, 1, 10);   vector_set(b2, 1, 10);
    vector_set(a2, 2, 10);   vector_set(b2, 2, 10);
    assert(vector_isequal(a2, b2) == 0);

    // test 3: same length, same values
    vector_set(a2, 0, 10);   vector_set(b2, 0, 10);
    vector_set(a2, 1, 10);   vector_set(b2, 1, 10);
    vector_set(a2, 2, 10);   vector_set(b2, 2, 10);
    assert(vector_isequal(a2, b2) == 1);

    printf("Tests for `vector_isequal` completed: no issues found.\n");
}


// Sets a_i to a_i+b_i for i=0,...,a->size-1
int vector_add(vector *a, vector *b){
    return gsl_vector_add(a, b);
}

// Sets a_i to a_i-b_i for i=0,...,a->size-1
int vector_sub(vector *a, vector *b){
    return gsl_vector_sub(a, b);
}


////////////////////////  short


vector_short *vector_short_alloc(size_t n){
    vector_short *v = gsl_vector_short_alloc(n);
    return v;
}

vector_short *vector_short_calloc(size_t n){
    vector_short *v = gsl_vector_short_calloc(n);
    return v;
}

void vector_short_free(vector_short *v){
    gsl_vector_short_free(v);
}

void vector_short_set(vector_short *v, size_t i, short x){
    gsl_vector_short_set(v, i, x);
}

short vector_short_get(vector_short *v, size_t i){
    return gsl_vector_short_get(v, i);
}

// Pointer to `i`th element of a vector `v`.
short *vector_short_ptr(vector_short *v, size_t i){
    return gsl_vector_short_ptr(v, i);
}

// Checks if a vector only contains 0-1 entries
int vector_short_is_zero_one(vector_short *v){
    int zero_one_counter = 0;
    for (size_t i=0; i<v->size; i++){
        if ((vector_short_get(v, i) == 0) || (vector_short_get(v, i) == 1)){
            zero_one_counter++;
        }
    }
    return (zero_one_counter == v->size);
}

void test_vector_short_is_zero_one(void){
    vector_short *v = vector_short_alloc(3);
    // test 1: 0-1 vector
    vector_short_set(v, 0, 0); vector_short_set(v, 1, 1); vector_short_set(v, 2, 0);
    assert(vector_short_is_zero_one(v) == 1);
    // test 2: not 0-1 vector
    vector_short_set(v, 0, 0); vector_short_set(v, 1, 1); vector_short_set(v, 2, 3);
    assert(vector_short_is_zero_one(v) == 0);
    vector_short_free(v);
    printf("Tests for `vector_short_is_zero_one` completed: no issues found.\n");
}

// Sums elements of the vector `v`.
int vector_short_sum(vector_short *v){
    int sum = 0;
    for (size_t i=0; i<v->size; i++){
        sum += vector_short_get(v, i);
    }
    return sum;
}


void test_vector_short_sum(void){
    vector_short *v = vector_short_alloc(2);
    vector_short_set(v, 0, 10); vector_short_set(v, 1, 3);
    assert(vector_short_sum(v) == 13.0);
    vector_short_set(v, 0, 10); vector_short_set(v, 1, -3);
    assert(vector_short_sum(v) == 7.0);
    vector_short_free(v);
    printf("Tests for `vector_short_sum` completed: no issues found.\n");
}

// Computes the mean of elements in `v`.
double vector_short_mean(vector_short *v){
    double mean = 0;
    double denom = (v->size + 0.0);
    for (size_t i=0; i<v->size; i++){
        mean += vector_short_get(v, i) / denom;
    }
    return mean;
}

void test_vector_short_mean(void){
    size_t n = 5;
    vector_short *v = vector_short_alloc(n);
    short v0, v1, v2, v3, v4;
    // test 1
    v0 = 0; 
    v1 = 0;
    v2 = 0;
    v3 = 0;
    v4 = 0;
    vector_short_set(v, 0, v0);
    vector_short_set(v, 1, v1);
    vector_short_set(v, 2, v2);
    vector_short_set(v, 3, v3);
    vector_short_set(v, 4, v4);
    double mean = (v0 + v1 + v2 + v3 + v4) / (n + 0.0);
    double result = vector_short_mean(v);
    assert(fabs(mean - result) <= pow(10, VECTOR_TEST_DOUBLE_TOLERANCE_NPOWER));
    // cleanup
    vector_short_free(v);
    printf("Tests for `vector_short_mean` completed: no issues found.\n");
}

////////////////////////  int



vector_int *vector_int_alloc(size_t n){
    vector_int *v = gsl_vector_int_alloc(n);
    return v;
}

vector_int *vector_int_calloc(size_t n){
    vector_int *v = gsl_vector_int_calloc(n);
    return v;
}

void vector_int_free(vector_int *v){
    gsl_vector_int_free(v);
}

void vector_int_set(vector_int *v, size_t i, int x){
    gsl_vector_int_set(v, i, x);
}

int vector_int_get(vector_int *v, size_t i){
    return gsl_vector_int_get(v, i);
}

// Pointer to `i`th element of a vector `v`.
int *vector_int_ptr(vector_int *v, size_t i){
    return gsl_vector_int_ptr(v, i);
}

int vector_int_equal(vector_int *v1, vector_int *v2){
    return gsl_vector_int_equal(v1, v2);
}

// Checks if a vector only contains 0-1 entries
int vector_int_is_zero_one(vector_int *v){
    int zero_one_counter = 0;
    for (size_t i=0; i<v->size; i++){
        if ((vector_int_get(v, i) == 0) || (vector_int_get(v, i) == 1)){
            zero_one_counter++;
        }
    }
    return (zero_one_counter == v->size);
}

void test_vector_int_is_zero_one(void){
    vector_int *v = vector_int_alloc(3);
    // test 1: 0-1 vector
    vector_int_set(v, 0, 0); vector_int_set(v, 1, 1); vector_int_set(v, 2, 0);
    assert(vector_int_is_zero_one(v) == 1);
    // test 2: not 0-1 vector
    vector_int_set(v, 0, 0); vector_int_set(v, 1, 1); vector_int_set(v, 2, 3);
    assert(vector_int_is_zero_one(v) == 0);
    vector_int_free(v);
}

// Sums elements of a vector
int vector_int_sum(vector_int *v){
    int sum = 0;
    for (size_t i=0; i<v->size; i++){
        sum += vector_int_get(v, i);
    }
    return sum;
}

void test_vector_int_sum(void){
    vector_int *v = vector_int_alloc(2);
    vector_int_set(v, 0, 10); vector_int_set(v, 1, 3);
    assert(vector_int_sum(v) == 13.0);
    vector_int_set(v, 0, 10); vector_int_set(v, 1, -3);
    assert(vector_int_sum(v) == 7.0);
    vector_int_free(v);
}

// Sets a_i to a_i+b_i for i=0,...,a->size-1
int vector_int_add(vector_int *a, vector_int *b){
    return gsl_vector_int_add(a, b);
}


////////////////////////  test




void test_vector(void){
    test_vector_is_zero_one();
    test_vector_sum();
    test_vector_isequal();
    test_vector_short_is_zero_one();
    test_vector_short_sum();
    test_vector_int_is_zero_one();
    test_vector_int_sum();
    test_vector_short_mean();
    printf("=== Tests: all tests for `vector` completed: no issues found.\n");
}