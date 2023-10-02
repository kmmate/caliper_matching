/* Header for dynamicarray.c */

#ifndef DYNAMICARRAY_H_
#define DYNAMICARRAY_H_

// includes
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


 // structs
 typedef struct dynamicarray_int {
    size_t capacity;    // allocated size, i.e. capacity
    size_t usage;  // actually used size, i.e. number of objects currently contained
    int *data;
 } dynamicarray_int;

// functions
dynamicarray_int *dynamicarray_int_alloc(size_t capacity);
void dynamicarray_int_free(dynamicarray_int *array);
int dynamicarray_int_get(dynamicarray_int *array, size_t i);
int dynamicarray_int_isempty(dynamicarray_int *a);
void dynamicarray_int_append(dynamicarray_int *array, int x);
int dynamicarray_int_equal(dynamicarray_int *a, dynamicarray_int *b);
int dynamicarray_int_permute_equal(dynamicarray_int *a, dynamicarray_int *b);
dynamicarray_int *dynamicarray_int_unique_elements(dynamicarray_int *a);
int dynamicarray_int_contains(dynamicarray_int *a, int x);

// tests
void test_dynamicarray(void);


#endif // DYNAMICARRAY_H_