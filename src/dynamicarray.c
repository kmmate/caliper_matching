////////////////////////////////////////////////////////////////////////////
//
// * FILE: dynamicarray.c
// * DESCRIPTION:
//    Dynamic arrays.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 31/oct/2022
//
//////////////////////////////////////////////////////////////////////////

#include "dynamicarray.h"

dynamicarray_int *dynamicarray_int_alloc(size_t capacity){
    if (capacity == 0){
        fprintf(stderr, "%s: line %d: `dynamicarray_int_alloc` capacity must be at least one, now is zero.\n", __FILE__, __LINE__);
        return NULL;
    }
    dynamicarray_int *array = malloc(sizeof(dynamicarray_int));
    int *data = malloc(capacity * sizeof(int));
    array->capacity = capacity;
    array->usage = 0;
    array->data = data;
    return array;
}

void dynamicarray_int_free(dynamicarray_int *array){
    free(array->data);
    array->capacity = 0;
    array->usage = 0;
    free(array);
}

int dynamicarray_int_get(dynamicarray_int *array, size_t i){
    if (i >= array->capacity){
        fprintf(stderr, "%s: line %d: `dynamicarray_int_get` index %zu is out of range for array of capacity %zu.\n", __FILE__, __LINE__, i, array->capacity);
        exit(EXIT_FAILURE);
    }
    return array->data[i];
}

int dynamicarray_int_isempty(dynamicarray_int *a){
    return a->usage == 0;
}

// Appends `x` to the end of array `array`. Expands array if necessary.
void dynamicarray_int_append(dynamicarray_int *array, int x){
    if (array->usage == array->capacity){ // if capacity is full...
        array->capacity *= 2;  // double capacity
        array->data = realloc(array->data, (array->capacity) * sizeof(int));
    }
    (array->usage)++; // increase usage by one
    array->data[(array->usage) - 1] = x;  // assign new value
}

void test_dynamicarray_int_append(void){
    // test 1
    dynamicarray_int *a = dynamicarray_int_alloc(1);
    assert(a->capacity == 1);
    assert(a->usage == 0);
    // test 2
    dynamicarray_int_append(a, 20);
    assert(a->usage == 1);
    // test 3
    assert(dynamicarray_int_get(a, 0)==20);
    // test 4
    dynamicarray_int *array = dynamicarray_int_alloc(2);
    dynamicarray_int_append(array, 20);
    dynamicarray_int_append(array, 30);
    dynamicarray_int_append(array, 50);
    // for (int i=0; i<3; i++){
    //     printf("dynamicarray_get(%d) = %d.\n", i , dynamicarray_int_get(array, i));
    // }
    assert((dynamicarray_int_get(array, 0) == 20) && (dynamicarray_int_get(array, 1) == 30) && (dynamicarray_int_get(array, 2) == 50));
    printf("Tests for `dynamicarray_int_append` completed: no issues found.\n");
}

// Checks if `a` contains `x` in `a->data` up to index`a->usage-1`
int dynamicarray_int_contains(dynamicarray_int *a, int x){
    if (a->usage==0){
        return 0;
    } else {
        int contains = 0;
        for (int i=0; i<a->usage; i++){
            // printf("a->data[%d] = %d.\n", i, a->data[i]);
            contains += a->data[i] == x;
        }
        return contains;
    }
}

void test_dynamicarray_int_contains(void){
    // test 1
    dynamicarray_int *a = dynamicarray_int_alloc(1);
    assert(dynamicarray_int_contains(a, 10)==0);
    // test 2
    dynamicarray_int_append(a, 20);
    assert(dynamicarray_int_contains(a, 10)==0);
    // test 3
    assert(dynamicarray_int_contains(a, 20)==1);
    // test 4
    dynamicarray_int_append(a, 30);
    assert(dynamicarray_int_contains(a, 20) && dynamicarray_int_contains(a, 30));
    printf("Tests for `dynamicarray_int_contains` completed: no issues found.\n");
}


/* 
Returns one if and only if arrays `a` and `b` are equal; zero otherwise.
They are considered equal if and only if both
- they have the usage >= 1; and
- all their data entries are equal, up until index usage-1,
  that is, if a->data[i] = b->data[i] for all i=0,...,usage-1.
*/
int dynamicarray_int_equal(dynamicarray_int *a, dynamicarray_int *b){
    if (a->usage != b->usage){
        return 0;
    } else if (a->usage==0){
        return 0;  // not equal if they do not hold any data
    } else {
        int equal = 1;
        for (int i=0; i<a->usage; i++){
            equal *= dynamicarray_int_get(a, i) == dynamicarray_int_get(b, i); //(a->data[i] == b->data[i]);
        }
        return equal; 
    }
}

void test_dynamicarray_equal(void){
    dynamicarray_int *a = dynamicarray_int_alloc(1);
    dynamicarray_int *b = dynamicarray_int_alloc(1);
    // test 1:  same capacity, no data
    assert(dynamicarray_int_equal(a, b)==0);
    // test 2: same capacity, same data
    dynamicarray_int_append(b, 10);
    dynamicarray_int_append(a, 10);
    assert(dynamicarray_int_equal(a, b)==1);
    // test 3: same "extended capacity", same data
    dynamicarray_int_append(a, 20);
    dynamicarray_int_append(b, 20);
    assert(dynamicarray_int_equal(a, b)==1);
    // test 4: same "extended capacity", different data
    dynamicarray_int_append(a, 12);
    dynamicarray_int_append(b, 10);
    assert(dynamicarray_int_equal(a, b)==0);
    // test 5: different capacity
    assert(dynamicarray_int_equal(dynamicarray_int_alloc(1), dynamicarray_int_alloc(2))==0);
    // test 6: same capacity, different usage
    a = dynamicarray_int_alloc(2);
    b = dynamicarray_int_alloc(2);
    dynamicarray_int_append(a, 12);
    assert(dynamicarray_int_equal(a, b)==0);
    printf("Tests for `dynamicarray_int_equal` completed: no issues found.\n");
}



dynamicarray_int *dynamicarray_int_unique_elements(dynamicarray_int *a){
    dynamicarray_int *unique_elements = dynamicarray_int_alloc(1);
    if (a->usage>0){
        dynamicarray_int_append(unique_elements, dynamicarray_int_get(a, 0));
        for (int i=1; i<a->usage; i++){
            // if `i`th element has not been added to unique elements, ad it now
            if (!dynamicarray_int_contains(unique_elements, dynamicarray_int_get(a, i))){
                dynamicarray_int_append(unique_elements, dynamicarray_int_get(a, i));
            }
        }
    }
    return unique_elements;
}

void test_dynamicarray_int_unique_elements(void){
    // test 1
    dynamicarray_int *a = dynamicarray_int_alloc(1);
    assert(dynamicarray_int_unique_elements(a)->usage == 0); // empty array
    // test 2
    dynamicarray_int_append(a, 10);
    dynamicarray_int_append(a, 10);  
    dynamicarray_int_append(a, 42);
    assert(dynamicarray_int_unique_elements(a)->usage == 2); // two unique elements
    // test 3
    assert((dynamicarray_int_get(dynamicarray_int_unique_elements(a), 0) == 10) && dynamicarray_int_get(dynamicarray_int_unique_elements(a), 1) == 42);
    printf("Tests for `dynamicarray_int_unique_elements` completed: no issues found.\n");
}


/* 
Returns one if and only if arrays `a` and `b` are equal up to permutations; zero otherwise.
They are considered equal up to permutations if and only if both
- they have the same usage >= 1; and
- a->data and b->data are permutations of one another.
*/
int dynamicarray_int_permute_equal(dynamicarray_int *a, dynamicarray_int *b){
    if (a->usage != b->usage){
        return 0;
    } else if (a->usage == 0){
        return 0;  // not equal if they do not hold any data
    } else if (a->usage == 1){
        return dynamicarray_int_get(a, 0) == dynamicarray_int_get(b, 0);
    } else {
        dynamicarray_int *unique_elements_a = dynamicarray_int_unique_elements(a);
        dynamicarray_int *unique_elements_b = dynamicarray_int_unique_elements(a);
        if (unique_elements_a->usage != unique_elements_b->usage){  // different number of unique elements
            printf("im here!!!");
            return 0;
        } else {
            int equal = 1;
            // count number of a times a unique element of `a` is present in `a` and `b`
            for (int i=0; i<unique_elements_a->usage; i++){
                int count_a, count_b;
                count_a = count_b = 0;
                for (int j=0; j<a->usage; j++){
                    count_a += dynamicarray_int_get(a, j) == dynamicarray_int_get(unique_elements_a, i);
                    count_b += dynamicarray_int_get(b, j) == dynamicarray_int_get(unique_elements_a, i);                    
                }
                equal *= count_a == count_b;
            }
            return equal; 
        }
    }
}


void test_dynamicarray_int_permute_equal(void){
    dynamicarray_int *a = dynamicarray_int_alloc(1);
    dynamicarray_int *b = dynamicarray_int_alloc(1);
    // test 1:  same capacity, no data
    assert(dynamicarray_int_permute_equal(a, b)==0);
    // test 2: same capacity, same data
    dynamicarray_int_append(b, 10);
    dynamicarray_int_append(a, 10);
    assert(dynamicarray_int_permute_equal(a, b)==1);
    // test 3: same capacity, different usage
    dynamicarray_int_append(a, 30);
    assert(dynamicarray_int_permute_equal(a, b)==0);
    // test 4: same capacity, same data in same order
    a = dynamicarray_int_alloc(1);
    b = dynamicarray_int_alloc(1);
    dynamicarray_int_append(b, 10);
    dynamicarray_int_append(a, 10);
    dynamicarray_int_append(b, 40);
    dynamicarray_int_append(a, 40);
    dynamicarray_int_append(b, 40);
    dynamicarray_int_append(a, 40);
    assert(dynamicarray_int_permute_equal(a, b)==1);
    // test 5: same capacity, same data in different order
    dynamicarray_int_append(b, 60);
    dynamicarray_int_append(a, 70);
    dynamicarray_int_append(b, 70);
    dynamicarray_int_append(a, 60);
    assert(dynamicarray_int_permute_equal(a, b)==1);
    // test 6
    a = dynamicarray_int_alloc(4);
    b = dynamicarray_int_alloc(4);
    dynamicarray_int_append(a, 2);
    dynamicarray_int_append(a, 3);
    dynamicarray_int_append(b, 3);
    dynamicarray_int_append(b, 2);
    assert(dynamicarray_int_permute_equal(a, b)==1);
    printf("Tests for `dynamicarray_int_permute_equal` completed: no issues found.\n");
}


void test_dynamicarray(void){
    test_dynamicarray_int_append();
    test_dynamicarray_int_contains();
    test_dynamicarray_int_unique_elements();
    test_dynamicarray_equal();
    test_dynamicarray_int_permute_equal();
    printf("=== Tests: all tests for `dynamicarray` completed: no issues found.\n");
}