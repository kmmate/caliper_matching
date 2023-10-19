//////////////////////////////////////////////////////////////////////////
//
// * FILE: cm.c
// * DESCRIPTION:
//    Caliper matching library.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 06/oct/2023
// * COMPILE:
//   Mac: gcc -pthread -lgsl cm.c vector.c matrix.c propscore.c -o cm
//  * NOTE: To add new propensity score models, follow the instructions
//          indicated by the comment `EDIT_TO_ADD`. In particular,
//          (0) modify `propscore.c` and `propscore.h` as described
//              therein; and
//          (1) adjust `cm_number_of_modeltypes`; and
//          (2) add the name of the model to the
//          `CM_FOREACH_MODEL` definition; and
//          (3) add the new model to the g/ginv/gderiv/ginvderiv arrays.
//
//////////////////////////////////////////////////////////////////////////

#include "cm.h"

//////////////////      Propensity score utilities

// EDIT_TO_ADD: increase number by J if J new models are added
const int cm_number_of_modeltypes = 2;

// EDIT_TO_ADD: add model name
#define CM_FOREACH_MODEL(MODEL) \
    MODEL(logit)             \
    MODEL(probit)            \
    // MODEL (mynewmodel) \  // comment out this line to add mynewmodel

#define CM_GENERATE_ENUM(ENUM) ENUM,
#define CM_GENERATE_STRING(STRING) #STRING,

#define CM_STR(x) #x
#define CM_XSTR(x) CM_STR(x) // use to print an enum item as string, e.g. CM_XSTR(logit) prints "logit".

enum CM_MODEL_TYPES{
    CM_FOREACH_MODEL(CM_GENERATE_ENUM)
};

const char *CM_MODEL_TYPES_CHAR[] = {
    CM_FOREACH_MODEL(CM_GENERATE_STRING)};

typedef double (*fnc_ptr)(double);
// EDIT_TO_ADD: add [mynewmodel] = propscore_<g/ginv/gderiv/ginvderiv>_mynewmodel to
//              each of the arrays below.
// order can be arbitrary, [modeltype]=. syntax makes sure it matches those in CM_MODEL_TYPES
fnc_ptr propscore_g[cm_number_of_modeltypes] = {[logit] = propscore_g_logit, [probit] = propscore_g_probit};
fnc_ptr propscore_ginv[cm_number_of_modeltypes] = {[logit] = propscore_ginv_logit, [probit] = propscore_ginv_probit};
fnc_ptr propscore_gderiv[cm_number_of_modeltypes] = {[logit] = propscore_gderiv_logit, [probit] = propscore_gderiv_probit};
fnc_ptr propscore_ginvderiv[cm_number_of_modeltypes] = {[logit] = propscore_ginvderiv_logit, [probit] = propscore_ginvderiv_probit};

// const char *modeltypes[] = {"logit", "probit"};
//  note: order must correspond to that in *modeltypes[]
// enum mdltypes {logit, probit};
// enum mdltypes mdltype;

// possible errors in CMModel
enum CMModel_ERRORS{
    CMModel_ERROR_DIMENSIONMISMATCH = 1, // reserve zero for no error
    CMModel_ERROR_MODELTYPE,
    CMModel_ERROR_THETALENGTH,
    CMModel_ERROR_CALIPER,
    CMModel_ERROR_TREATMENTZEROONE,
    CMModel_ERROR_TREATMENTIMBALANCE,
    CMModel_ERROR_BETA,
    CMModel_ERROR_ALPHA,
    CMModel_ERROR_BETAALPHA,
    CMModel_ERROR_KAPPAA,
    CMModel_ERROR_KAPPAGAMMA,
    CMModel_ERROR_KAPPAGAMMADERIVATIVE
};

// possible errors in CMModelKnownPropscore
enum CMModelKnownPropscore_ERRORS{
    CMModelKnownPropscore_ERROR_DIMENSIONMISMATCH = 1, // reserve zero for no error
    CMModelKnownPropscore_ERROR_CALIPER,
    CMModelKnownPropscore_ERROR_TREATMENTZEROONE,
    CMModelKnownPropscore_ERROR_TREATMENTIMBALANCE,
    CMModelKnownPropscore_ERROR_PROPSCOREOUTOFRANGE,
    CMModelKnownPropscore_ERROR_BETA,
    CMModelKnownPropscore_ERROR_ALPHA,
    CMModelKnownPropscore_ERROR_BETAALPHA,
    CMModelKnownPropscore_ERROR_KAPPAA,
    CMModelKnownPropscore_ERROR_KAPPAGAMMA
};

pthread_mutex_t mutex_lock = PTHREAD_MUTEX_INITIALIZER;

// Finds index of propensity score model type `modeltype` in CM_MODEL_TYPES_CHAR.
int cm_find_modeltype_index(char *modeltype){
    for (size_t i = 0; i < cm_number_of_modeltypes; i++){
        if (!strcmp(modeltype, CM_MODEL_TYPES_CHAR[i])){
            return i;
        }
    }
    return cm_number_of_modeltypes + 1; // out of range
}



//////////////////      Auxiliary functions: generic functions


// Minimum of two integers.
int cm_int_min(int a, int b){
    return a <= b ? a : b;
}

void test_cm_int_min(void){
    assert(cm_int_min(0, 0) == 0);
    assert(cm_int_min(3, 4) == 3);
    assert(cm_int_min(3, -4) == -4);
    assert(cm_int_min(-10, 5) == -10);
    assert(cm_int_min(-3, -4) == -4);
    printf("Tests for `cm_int_min` completed: no issues found.\n");
}

// Maximum of two integers.
int cm_int_max(int a, int b){
    return a <= b ? b : a;
}

void test_cm_int_max(void){
    assert(cm_int_max(3, 4) == 4);
    assert(cm_int_max(3, -4) == 3);
    assert(cm_int_max(-10, 5) == 5);
    assert(cm_int_max(5, -10) == 5);
    assert(cm_int_max(-3, -4) == -3);
    printf("Tests for `cm_int_max` completed: no issues found.\n");
}

// Minimum of two doubles.
double cm_d_min(double a, double b){
    return a <= b ? a : b;
}

void test_cm_d_min(void){
    assert(fabs(cm_d_min(0.0, 0.0) - 0.0) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_min(3.4, 4.5) - 3.4) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_min(3, -4.3) - (-4.3)) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_min(-10.0, 5.5) - (-10.0)) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_min(-3, -4.42) - (-4.42)) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    printf("Tests for `cm_d_min` completed: no issues found.\n");
}

// Maximum of two doubles.
double cm_d_max(double a, double b){
    return a <= b ? b : a;
}


void test_cm_d_max(void){
    assert(fabs(cm_d_max(0.0, 0) - 0.0) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_max(3.4, 4.5) - 4.5) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_max(3.0, -4.3) - 3.0) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_max(-10.0, 5.5) - 5.5) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(cm_d_max(-3.0, -4.42) - (-3.0)) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    printf("Tests for `cm_d_max` completed: no issues found.\n");
}

// Computes absolute difference between two doubles.
double cm_abs_distance(double x, double y){
    double out = fabs((x) - (y));
    return out;
}

void test_cm_abs_distance(void){
    double x, y;
    // test 1
    x = -0.5;
    y = 0.8;
    assert((cm_abs_distance(x, y) == cm_abs_distance(y, x)) && (cm_abs_distance(x, y) >= 0));
    // test 2
    x = 0.5;
    y = 0.8;
    assert((cm_abs_distance(x, y) == cm_abs_distance(y, x)) && (cm_abs_distance(x, y) >= 0));
    // test 3
    x = -0.5;
    y = -0.8;
    assert((cm_abs_distance(x, y) == cm_abs_distance(y, x)) && (cm_abs_distance(x, y) >= 0));
    // test 4
    x = 1.5;
    y = -0.8;
    assert((cm_abs_distance(x, y) == cm_abs_distance(y, x)) && (cm_abs_distance(x, y) >= 0));
    printf("Tests for `cm_abs_distance` completed: no issues found.\n");
}


// Returns whether `x` is in the closed interval [`low`,`high`].
int cm_inrange(double x, double low, double high){
    assert(low <= high);
    return ((low <= x) && (x <= high));
}

void test_cm_inrange(void){
    assert(cm_inrange(2, -1, 4));
    assert(!cm_inrange(2, -4, -1));
    assert(cm_inrange(2, 1, 4));
    assert(!cm_inrange(-2, -1, 4));
    assert(!cm_inrange(-2, -34, -3));
}


/*
Comparison function of two double pointers. Used by `vector_sort_index_ptr`.
Whichever pointer points to a larger value is considered larger.
*/
int vector_sort_index_ptr_compare(const void *p_left, const void *p_right){
    const double *left = *(const double **)p_left;
    const double *right = *(const double **)p_right;
    return (*left > *right) - (*left < *right);
}

/*
Returns an array `sort_index` of double pointers such that larger-index-elements in `sort_index`
point to larger elements in `v`. That is, sort_index[0] is a pointer to the smallest element in `v`,
so that the smallest element in `v` has value *sort_index[0], and sort_index[v->size-1] is a pointer
to the largest element in `v`, so that the largest element in `v` has value *sort_index[v->size-1].
The `v` itself is not changed.

Useful when `sort_index` is only used to access elements in `v`.
Memory efficient, does not allocate extra memory.

See also: `vector_sort_index_ptr`.
*/
double **vector_sort_index_ptr(vector *v){
    double **sort_index = malloc((v->size) * (sizeof(double *)));
    for (int i = 0; i < v->size; i++){
        sort_index[i] = vector_ptr(v, i);
    }
    qsort(sort_index, v->size, sizeof(sort_index[0]), vector_sort_index_ptr_compare);
    return sort_index;
}

void test_vector_sort_index_ptr(void){
    double v0, v1, v2, v3;
    int cond1, cond2, cond3, cond4, condv1, condv2, condv3, condv4;
    vector *v = vector_alloc(4);
    v0 = 1.2;
    v1 = 1.0;
    v2 = 1.3;
    v3 = 1.1;
    vector_set(v, 0, v0);
    vector_set(v, 1, v1);
    vector_set(v, 2, v2);
    vector_set(v, 3, v3);
    double **sort_index = vector_sort_index_ptr(v);
    // for (int i=0; i<4; i++){
    //     printf("*sort_index[%d] = %f\n", i, *sort_index[i]);
    //     printf("v[%d] = %f\n", i, vector_get(v, i));
    // }
    cond1 = (*sort_index[0] == v1); // first element should point to least element in v
    cond4 = (*sort_index[1] == v3);
    cond4 = (*sort_index[2] == v0);
    cond4 = (*sort_index[3] == v2);    // last element should point to largest element in v
    condv1 = (vector_get(v, 0) == v0); // check v is unchanged
    condv2 = (vector_get(v, 1) == v1);
    condv3 = (vector_get(v, 2) == v2);
    condv4 = (vector_get(v, 3) == v3);
    assert(cond1 && cond2 && cond3 && cond4 && condv1 && condv2 && condv3 && condv4);
    printf("Tests for `vector_sort_index_ptr` completed: no issues found.\n");
}

/*
Comparison function of two size_t objects.  Used by `vector_sort_index`.
Whichever pointer indexes a larger element in `v` is considered the larger.
*/
int vector_sort_index_compare(void *vv, const void *p_left, const void *p_right){
    vector *v = (vector *)vv;
    const size_t left = *(const size_t *)p_left;
    const size_t right = *(const size_t *)p_right;
    double v_left, v_right;
    v_left = vector_get(v, left);
    v_right = vector_get(v, right);
    return (v_left > v_right) - (v_left < v_right);
}

/*
Return a permutation of the indicies of `v` as an array `sort_index`. The permutation sorts `v`
in ascending order. That is, sort_index[0] is the index of the smallest element in `v`,
so that the smallest element in `v` has value vector_get(v, sort_index[0]);
and sort_index[v->size-1] is the index of the largest element in `v`, so that the largest element
in `v` has value vector_get(v, sort_index[v->size-1]). The `v` itself is not changed.

Useful when `sort_index` is not only used to access elements in `v` but also to index other vectors,
e.g. to permute a data set based on the sorted `v`.
Allocates extra memory of the order (v->size)*sizeof(size_t).

See also: `vector_sort_index_ptr`.
*/
size_t *vector_sort_index(vector *v){
    size_t *sort_index = malloc((v->size) * sizeof(size_t));
    for (int i = 0; i < v->size; i++){
        sort_index[i] = i;
    }
    qsort_r(sort_index, v->size, sizeof(sort_index[0]), (void *)v, vector_sort_index_compare);
    return sort_index;
}

void test_vector_sort_index(void){
    double v0, v1, v2, v3;
    int cond1, cond2, cond3, cond4, condv1, condv2, condv3, condv4;
    vector *v = vector_alloc(4);
    v0 = 1.2;
    v1 = 1.0;
    v2 = 1.3;
    v3 = 1.1;
    vector_set(v, 0, v0);
    vector_set(v, 1, v1);
    vector_set(v, 2, v2);
    vector_set(v, 3, v3);
    size_t *sort_index = vector_sort_index(v);
    // for (int i=0; i<4; i++){
    //     printf("sort_index[%d] = %zu\n", i, sort_index[i]);
    //     printf("v[%d] = %f\n", i, vector_get(v, i));
    // }
    cond1 = (vector_get(v, sort_index[0]) == v1); // first element should index the least element in v
    cond4 = (vector_get(v, sort_index[1]) == v3);
    cond4 = (vector_get(v, sort_index[2]) == v0);
    cond4 = (vector_get(v, sort_index[3]) == v2); // last element should index the largest element in v
    condv1 = (vector_get(v, 0) == v0);            // check v is unchanged
    condv2 = (vector_get(v, 1) == v1);
    condv3 = (vector_get(v, 2) == v2);
    condv4 = (vector_get(v, 3) == v3);
    assert(cond1 && cond2 && cond3 && cond4 && condv1 && condv2 && condv3 && condv4);
    printf("Tests for `vector_sort_index` completed: no issues found.\n");
}


//////////////////      Auxiliary functions: specific functions


// Check if before-instantiation-inputs to CMModel are correct; return zero if and only if all inputs are correct.
int cm_cmmodel_init_check_values(CMModel *cmm, int test_mode){
    // if test_mode==1, cm_cm_input_check does not print error messages, only returns error code
    // data dimension check
    if (!((cmm->d->size == cmm->y->size) && (cmm->y->size == cmm->x->size1))){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: dimension mismatch between `y`, `d`, `x`.\n", __FILE__, __LINE__);
        }
        return CMModel_ERROR_DIMENSIONMISMATCH;
    }
    // modeltype check
    int is_modeltype_in_modeltypes = 1;
    for (int i = 0; i < cm_number_of_modeltypes; i++){
        is_modeltype_in_modeltypes *= strcmp(cmm->modeltype, CM_MODEL_TYPES_CHAR[i]);
    }
    if (is_modeltype_in_modeltypes){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: `modeltype` '%s' is not supported for propensity score.\n", __FILE__, __LINE__, cmm->modeltype);
            printf("Available options are:\n");
            for (int i = 0; i < cm_number_of_modeltypes; i++){
                printf("%s\n", CM_MODEL_TYPES_CHAR[i]);
            }
        }
        return CMModel_ERROR_MODELTYPE;
    }
    // theta length check
    if (cmm->theta->size != (cmm->x->size2 + 1)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: length of `theta` must be larger than that of `x` exactly by one; current length of `theta` is %zu.\n", __FILE__, __LINE__, cmm->theta->size);
        }
        return CMModel_ERROR_THETALENGTH;
    }
    // caliper check
    if (cmm->delta < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: caliper `delta` must be zero to use default value, or be srictly positive; now it is %f.\n", __FILE__, __LINE__, cmm->delta);
        }
        return CMModel_ERROR_CALIPER;
    }
    // treatment vector: zero/one check
    for (int i = 0; i < cmm->d->size; i++){
        if (((int)vector_short_get(cmm->d, i) != 0) && ((int)vector_short_get(cmm->d, i) != 1)){
            if (!test_mode)
            {
                fprintf(stderr, "%s: line %d: CMModel input check failed: treatment vector `d` must contain elements zero or one; d[%d] has value %d.\n", __FILE__, __LINE__, i, (int)vector_short_get(cmm->d, i));
            }
            return CMModel_ERROR_TREATMENTZEROONE;
        }
    }
    // treatment vector: not enough treated or control units
    int nr_treated = (int)vector_short_sum(cmm->d);
    if (nr_treated == 0 || nr_treated == cmm->d->size){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: imbalanced number of treated and control units; there are %d treated units, %lu control units.\n", __FILE__, __LINE__, nr_treated, (cmm->d->size) - nr_treated);
        }
        return CMModel_ERROR_TREATMENTIMBALANCE;
    }
    // variance estimation: bandwidth exponent
    if ((cmm->beta < 0) || (cmm->beta >= 1.0 / 4)){
        if (!test_mode)
        {
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `beta` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->beta);
        }
        return CMModel_ERROR_BETA;
    }
    // variance estimation: truncation sequence exponent
    if ((cmm->alpha < 0) || (cmm->alpha >= 1.0 / 4)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `alpha` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->alpha);
        }
        return CMModel_ERROR_ALPHA;
    }
    // variance estimation: bandwidth and truncation sequence exponent
    if ((cmm->alpha > 0) && (cmm->beta > 0) && (cmm->alpha >= cmm->beta)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 < `alpha` < `beta` must be or at least one of `alpha` or `beta` must be zero; now alpha is %f, beta is %f.\n", __FILE__, __LINE__, cmm->alpha, cmm->beta);
        }
        return CMModel_ERROR_BETAALPHA;
    }
    // variance estimation: truncation sequence scale
    if (cmm->kappa_a < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `kappa_a` must be; now kappa_a is %f.\n", __FILE__, __LINE__, cmm->kappa_a);
        }
        return CMModel_ERROR_KAPPAA;
    }
    // variance estimation: bandwidth sequence scale
    if (cmm->kappa_gamma < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `kappa_gamma` must be; now kappa_gamma is %f.\n", __FILE__, __LINE__, cmm->kappa_gamma);
        }
        return CMModel_ERROR_KAPPAGAMMA;
    }
    // variance estimation: bandwidth sequence scale for derivatives
    if (cmm->kappa_gamma_derivative < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `kappa_gamma_derivative` must be; now kappa_gamma_derivative is %f.\n", __FILE__, __LINE__, cmm->kappa_gamma_derivative);
        }
        return CMModel_ERROR_KAPPAGAMMADERIVATIVE;
    }
    // all checks passed
    return 0;
}


// Check if after-instantiation-inputs to CMModel are correct; return zero if and only if all inputs are correct.
int cm_cmmodel_check_values(CMModel *cmm, int test_mode){
    // if test_mode==1, cm_cm_input_check does not print error messages, only returns error code
    // data dimension check
    if (!((cmm->d->size == cmm->y->size) && (cmm->y->size == cmm->x->size1))){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: dimension mismatch between `y`, `d`, `x`.\n", __FILE__, __LINE__);
        }
        return CMModel_ERROR_DIMENSIONMISMATCH;
    }
    // modeltype check
    int is_modeltype_in_modeltypes = 1;
    for (int i = 0; i < cm_number_of_modeltypes; i++){
        is_modeltype_in_modeltypes *= strcmp(cmm->modeltype, CM_MODEL_TYPES_CHAR[i]);
    }
    if (is_modeltype_in_modeltypes){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: `modeltype` '%s' is not supported for propensity score.\n", __FILE__, __LINE__, cmm->modeltype);
            printf("Available options are:\n");
            for (int i = 0; i < cm_number_of_modeltypes; i++){
                printf("%s\n", CM_MODEL_TYPES_CHAR[i]);
            }
        }
        return CMModel_ERROR_MODELTYPE;
    }
    // theta length check
    if (cmm->theta->size != (cmm->x->size2 + 1)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: length of `theta` must be larger than that of `x` exactly by one; current length of `theta` is %zu.\n", __FILE__, __LINE__, cmm->theta->size);
        }
        return CMModel_ERROR_THETALENGTH;
    }
    // caliper check
    if (cmm->delta <= 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: caliper `delta` must be strictly positive; now it is %f.\n", __FILE__, __LINE__, cmm->delta);
        }
        return CMModel_ERROR_CALIPER;
    }
    // treatment vector: zero/one check
    for (int i = 0; i < cmm->d->size; i++){
        if (((int)vector_short_get(cmm->d, i) != 0) && ((int)vector_short_get(cmm->d, i) != 1)){
            if (!test_mode)
            {
                fprintf(stderr, "%s: line %d: CMModel input check failed: treatment vector `d` must contain elements zero or one; d[%d] has value %d.\n", __FILE__, __LINE__, i, (int)vector_short_get(cmm->d, i));
            }
            return CMModel_ERROR_TREATMENTZEROONE;
        }
    }
    // treatment vector: not enough treated or control units
    int nr_treated = (int)vector_short_sum(cmm->d);
    if (nr_treated == 0 || nr_treated == cmm->d->size){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: imbalanced number of treated and control units; there are %d treated units, %lu control units.\n", __FILE__, __LINE__, nr_treated, (cmm->d->size) - nr_treated);
        }
        return CMModel_ERROR_TREATMENTIMBALANCE;
    }
    // variance estimation: bandwidth exponent
    if ((cmm->beta < 0) || (cmm->beta >= 1.0 / 4)){
        if (!test_mode)
        {
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `beta` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->beta);
        }
        return CMModel_ERROR_BETA;
    }
    // variance estimation: truncation sequence exponent
    if ((cmm->alpha < 0) || (cmm->alpha >= 1.0 / 4)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `alpha` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->alpha);
        }
        return CMModel_ERROR_ALPHA;
    }
    // variance estimation: bandwidth and truncation sequence exponent
    if ((cmm->alpha > 0) && (cmm->beta > 0) && (cmm->alpha >= cmm->beta)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 < `alpha` < `beta` must be or at least one of `alpha` or `beta` must be zero; now alpha is %f, beta is %f.\n", __FILE__, __LINE__, cmm->alpha, cmm->beta);
        }
        return CMModel_ERROR_BETAALPHA;
    }
    // variance estimation: truncation sequence scale
    if (cmm->kappa_a < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `kappa_a` must be; now kappa_a is %f.\n", __FILE__, __LINE__, cmm->kappa_a);
        }
        return CMModel_ERROR_KAPPAA;
    }
    // variance estimation: bandwidth sequence scale
    if (cmm->kappa_gamma < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `kappa_gamma` must be; now kappa_gamma is %f.\n", __FILE__, __LINE__, cmm->kappa_gamma);
        }
        return CMModel_ERROR_KAPPAGAMMA;
    }
    // variance estimation: bandwidth sequence scale for derivatives
    if (cmm->kappa_gamma_derivative < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModel input check failed: 0 <= `kappa_gamma_derivative` must be; now kappa_gamma_derivative is %f.\n", __FILE__, __LINE__, cmm->kappa_gamma_derivative);
        }
        return CMModel_ERROR_KAPPAGAMMADERIVATIVE;
    }
    // all checks passed
    return 0;
}

void test_cm_cmmodel_check_values(void){
    int test_mode = 1; // if 1, cm_cm_input_check does not print error messages, only returns error code
    CMModel *cm_model = malloc(sizeof(CMModel));
    // Data dimension mismatch
    int ret_code;
    int err_code = CMModel_ERROR_DIMENSIONMISMATCH;
    int n = 10;
    int k = 2;
    double delta = 0.1;
    char *modeltype = "probit";
    vector *theta = vector_alloc(k + 1);
    double beta = 0.0;
    double alpha = 0.0;
    double kappa_a = 0.0;
    double kappa_gamma = 0.0;
    double kappa_gamma_derivative = 0.0;
    // test 1
    vector *y = vector_alloc(n + 1);
    vector_short *d = vector_short_alloc(n);
    matrix *x = matrix_alloc(n, k);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    // test 2
    y = vector_alloc(n);
    d = vector_short_alloc(n + 1);
    x = matrix_alloc(n, k);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    // test 3
    y = vector_alloc(n);
    d = vector_short_alloc(n);
    x = matrix_alloc(n + 1, k);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    vector_free(theta);

    // Modeltype check
    err_code = CMModel_ERROR_MODELTYPE;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "lpm";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_alloc(n);
    x = matrix_alloc(n, k);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    vector_free(theta);

    // Theta dimension check
    err_code = CMModel_ERROR_THETALENGTH;
    n = 10;
    delta = 0.1;
    modeltype = "logit";
    k = 2;
    // test 1
    theta = vector_alloc(k);
    y = vector_alloc(n);
    d = vector_short_alloc(n);
    x = matrix_alloc(n, k);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(theta);
    // test 2
    theta = vector_alloc(k - 1);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    vector_free(theta);

    // Caliper check
    err_code = CMModel_ERROR_CALIPER;
    n = 10;
    k = 2;
    delta = -0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_alloc(n);
    x = matrix_alloc(n, k);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    vector_free(theta);

    // Treatment vector check
    // test 1: is vector zero/one
    err_code = CMModel_ERROR_TREATMENTZEROONE;
    n = 3;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 2);
    x = matrix_alloc(n, k);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_short_free(d);
    // test 2: balance of treated and control units
    err_code = CMModel_ERROR_TREATMENTIMBALANCE;
    d = vector_short_calloc(n);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_short_free(d);
    // test 3: balance of treated and control units
    err_code = CMModel_ERROR_TREATMENTIMBALANCE;
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    vector_short_set(d, 1, 1);
    vector_short_set(d, 2, 1);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    vector_free(theta);

    // Variance estimation: bandwidth exponent check
    // test 1: beta too small
    err_code = CMModel_ERROR_BETA;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = -1.0;
    alpha = 0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: beta too large
    err_code = CMModel_ERROR_BETA;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 1.0;
    alpha = 0.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: truncation exponent check
    // test 1: alpha too small
    err_code = CMModel_ERROR_ALPHA;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0;
    alpha = -1.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: alpha too large
    err_code = CMModel_ERROR_ALPHA;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.0;
    alpha = 1.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: bandwidth versus truncation exponent check
    // test 1: alpha too large compared to beta!=0
    err_code = CMModel_ERROR_BETAALPHA;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.1;
    alpha = 0.2;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // test 2: no error for beta=0
    err_code = 0;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0;
    alpha = 0.2;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // test 3: no error for alpha=0
    err_code = 0;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.1;
    alpha = 0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: truncation sequence scale
    // test 1: scale is negative
    err_code = CMModel_ERROR_KAPPAA;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.2;
    alpha = 0.1;
    kappa_a = -1.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: no error for kappa_a = 0
    err_code = 0;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.2;
    alpha = 0.1;
    kappa_a = 0.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: bandwidth sequence scale
    // test 1: kappa_gamma is negative
    err_code = CMModel_ERROR_KAPPAGAMMA;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.2;
    alpha = 0.1;
    kappa_a = 1.0;
    kappa_gamma = -1.3;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: no error for kappa_gamma = 0
    err_code = 0;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.2;
    alpha = 0.1;
    kappa_a = 0.0;
    kappa_gamma = 0.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: bandwidth sequence scale for derivative
    // test 1: kappa_gamma_derivative is negative
    err_code = CMModel_ERROR_KAPPAGAMMADERIVATIVE;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.2;
    alpha = 0.1;
    kappa_a = 1.0;
    kappa_gamma = 1.3;
    kappa_gamma_derivative = -4.5;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: no error for kappa_gamma_derivative = 0
    err_code = 0;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.2;
    alpha = 0.1;
    kappa_a = 0.0;
    kappa_gamma = 0.0;
    kappa_gamma_derivative = 0.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);


    // All checks pass
    err_code = 0;
    n = 10;
    k = 2;
    delta = 0.1;
    modeltype = "probit";
    theta = vector_alloc(k + 1);
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    x = matrix_alloc(n, k);
    beta = 0.2;
    alpha = 0.1;
    kappa_a = 0.01;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_model->modeltype = modeltype;
    cm_model->theta = theta;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    ret_code = cm_cmmodel_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    matrix_free(x);
    vector_free(theta);
    printf("Tests for `cm_cmmodel_check_values` completed: no issues found.\n");
}


// Check if before-instantiation-inputs to CMModelKnownPropscore are correct; return zero if and only if all inputs are correct.
int cm_cmmodelknownpropscore_init_check_values(CMModelKnownPropscore *cmm, int test_mode){
    // if test_mode==1, cm_cm_input_check does not print error messages, only returns error code
    // data dimension check
    if (!((cmm->d->size == cmm->y->size) && (cmm->y->size == cmm->propscore->size))){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: input check failed: dimension mismatch between `y`, `d`, `propscore` or length unequal to `n`.\n", __FILE__, __LINE__);
        }
        return CMModelKnownPropscore_ERROR_DIMENSIONMISMATCH;
    }
    // caliper check
    if (cmm->delta < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: caliper `delta` must be zero to use default value, or be srictly positive; now it is %f.\n", __FILE__, __LINE__, cmm->delta);
        }
        return CMModelKnownPropscore_ERROR_CALIPER;
    }
    // treatment vector: zero/one check
    for (int i = 0; i < cmm->d->size; i++){
        if (((int)vector_short_get(cmm->d, i) != 0) && ((int)vector_short_get(cmm->d, i) != 1)){
            if (!test_mode)
            {
                fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: treatment vector `d` must contain elements zero or one; d[%d] has value %d.\n", __FILE__, __LINE__, i, (int)vector_short_get(cmm->d, i));
            }
            return CMModelKnownPropscore_ERROR_TREATMENTZEROONE;
        }
    }
    // treatment vector: not enough treated or control units
    int nr_treated = (int)vector_short_sum(cmm->d);
    if (nr_treated == 0 || nr_treated == cmm->d->size){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: imbalanced number of treated and control units; there are %d treated units, %lu control units.\n", __FILE__, __LINE__, nr_treated, (cmm->d->size) - nr_treated);
        }
        return CMModelKnownPropscore_ERROR_TREATMENTIMBALANCE;
    }
    // propensity score: values not in  interval (0,1)
    for (int i = 0; i < cmm->propscore->size; i++){
        if (!((0 < vector_get(cmm->propscore, i)) && (vector_get(cmm->propscore, i) < 1)))
        {
            if (!test_mode)
            {
                fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: propensity score `propscore` must contain elements strictly between zero or one; propscore[%d] has value %f.\n", __FILE__, __LINE__, i, vector_get(cmm->propscore, i));
            }
            return CMModelKnownPropscore_ERROR_PROPSCOREOUTOFRANGE;
        }
    }
    // variance estimation: bandwidth exponent
    if ((cmm->beta < 0) || (cmm->beta >= 1.0 / 4)){
        if (!test_mode)
        {
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: 0 <= `beta` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->beta);
        }
        return CMModelKnownPropscore_ERROR_BETA;
    }
    // variance estimation: truncation sequence exponent
    if ((cmm->alpha < 0) || (cmm->alpha >= 1.0 / 4)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: 0 <= `alpha` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->alpha);
        }
        return CMModelKnownPropscore_ERROR_ALPHA;
    }
    // variance estimation: bandwidth and truncation sequence exponent
    if ((cmm->alpha > 0) && (cmm->beta > 0) && (cmm->alpha >= cmm->beta)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d:  CMModelKnownPropscore input check failed: 0 < `alpha` < `beta` must be or at least one of `alpha` or `beta` must be zero; now alpha is %f, beta is %f.\n", __FILE__, __LINE__, cmm->alpha, cmm->beta);
        }
        return CMModelKnownPropscore_ERROR_BETAALPHA;
    }
    // variance estimation: truncation sequence scale
    if (cmm->kappa_a < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnwonPropscore input check failed: 0 <= `kappa_a` must be; now kappa_a is %f.\n", __FILE__, __LINE__, cmm->kappa_a);
        }
        return CMModelKnownPropscore_ERROR_KAPPAA;
    }
    // variance estimation: bandwidth sequence scale
    if (cmm->kappa_gamma < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnwonPropscore input check failed: 0 <= `kappa_gamma` must be; now kappa_gamma is %f.\n", __FILE__, __LINE__, cmm->kappa_gamma);
        }
        return CMModelKnownPropscore_ERROR_KAPPAGAMMA;
    }
    // all checks passed
    return 0;
}


// Check if after-instantiation-inputs to CMModelKnownPropscore are correct; return zero if and only if all inputs are correct.
int cm_cmmodelknownpropscore_check_values(CMModelKnownPropscore *cmm, int test_mode){
    // if test_mode==1, cm_cm_input_check does not print error messages, only returns error code
    // data dimension check
    if (!((cmm->d->size == cmm->y->size) && (cmm->d->size == cmm->n) && (cmm->y->size == cmm->propscore->size))){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: input check failed: dimension mismatch between `y`, `d`, `propscore` or length unequal to `n`.\n", __FILE__, __LINE__);
        }
        return CMModelKnownPropscore_ERROR_DIMENSIONMISMATCH;
    }
    // caliper check
    if (cmm->delta <= 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: caliper `delta` must be strictly positive; now it is %f.\n", __FILE__, __LINE__, cmm->delta);
        }
        return CMModelKnownPropscore_ERROR_CALIPER;
    }
    // treatment vector: zero/one check
    for (int i = 0; i < cmm->d->size; i++){
        if (((int)vector_short_get(cmm->d, i) != 0) && ((int)vector_short_get(cmm->d, i) != 1)){
            if (!test_mode)
            {
                fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: treatment vector `d` must contain elements zero or one; d[%d] has value %d.\n", __FILE__, __LINE__, i, (int)vector_short_get(cmm->d, i));
            }
            return CMModelKnownPropscore_ERROR_TREATMENTZEROONE;
        }
    }
    // treatment vector: not enough treated or control units
    int nr_treated = (int)vector_short_sum(cmm->d);
    if (nr_treated == 0 || nr_treated == cmm->d->size){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: imbalanced number of treated and control units; there are %d treated units, %lu control units.\n", __FILE__, __LINE__, nr_treated, (cmm->d->size) - nr_treated);
        }
        return CMModelKnownPropscore_ERROR_TREATMENTIMBALANCE;
    }
    // propensity score: values not in  interval (0,1)
    for (int i = 0; i < cmm->propscore->size; i++){
        if (!((0 < vector_get(cmm->propscore, i)) && (vector_get(cmm->propscore, i) < 1)))
        {
            if (!test_mode)
            {
                fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: propensity score `propscore` must contain elements strictly between zero or one; propscore[%d] has value %f.\n", __FILE__, __LINE__, i, vector_get(cmm->propscore, i));
            }
            return CMModelKnownPropscore_ERROR_PROPSCOREOUTOFRANGE;
        }
    }
    // variance estimation: bandwidth exponent
    if ((cmm->beta < 0) || (cmm->beta >= 1.0 / 4)){
        if (!test_mode)
        {
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: 0 <= `beta` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->beta);
        }
        return CMModelKnownPropscore_ERROR_BETA;
    }
    // variance estimation: truncation sequence exponent
    if ((cmm->alpha < 0) || (cmm->alpha >= 1.0 / 4)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: 0 <= `alpha` < 1/4 must be; now it is %f.\n", __FILE__, __LINE__, cmm->alpha);
        }
        return CMModelKnownPropscore_ERROR_ALPHA;
    }
    // variance estimation: bandwidth and truncation sequence exponent
    if ((cmm->alpha > 0) && (cmm->beta > 0) && (cmm->alpha >= cmm->beta)){
        if (!test_mode){
            fprintf(stderr, "%s: line %d:  CMModelKnownPropscore input check failed: 0 < `alpha` < `beta` must be or at least one of `alpha` or `beta` must be zero; now alpha is %f, beta is %f.\n", __FILE__, __LINE__, cmm->alpha, cmm->beta);
        }
        return CMModelKnownPropscore_ERROR_BETAALPHA;
    }
    // variance estimation: truncation sequence scale
    if (cmm->kappa_a < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnwonPropscore input check failed: 0 <= `kappa_a` must be; now kappa_a is %f.\n", __FILE__, __LINE__, cmm->kappa_a);
        }
        return CMModelKnownPropscore_ERROR_KAPPAA;
    }
    // variance estimation: bandwidth sequence scale
    if (cmm->kappa_gamma < 0){
        if (!test_mode){
            fprintf(stderr, "%s: line %d: CMModelKnownPropscore input check failed: 0 <= `kappa_gamma` must be; now kappa_gamma is %f.\n", __FILE__, __LINE__, cmm->kappa_gamma);
        }
        return CMModelKnownPropscore_ERROR_KAPPAGAMMA;
    }
    // all checks passed
    return 0;
}

void test_cm_cmmodelknownpropscore_check_values(void){
    int test_mode = 1; // if 1, cm_cm_input_check does not print error messages, only returns error code
    CMModelKnownPropscore *cm_model = malloc(sizeof(CMModelKnownPropscore));
    // Data dimension mismatch
    int ret_code;
    int err_code = CMModelKnownPropscore_ERROR_DIMENSIONMISMATCH;
    int n = 10;
    double delta = 0.1; // default value that should not cause errors
    double beta = 0.0;  // default value that should not cause errors
    double alpha = 0.0;  // default value that should not cause errors
    double kappa_a = 0.0;  // default value that should not cause errors
    double kappa_gamma = 0.0;  // default value that should not cause errors
    // test 1
    vector *y = vector_alloc(n + 1);
    vector_short *d = vector_short_alloc(n);
    vector *propscore = vector_alloc(n);
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    vector_free(propscore);
    // test 2
    y = vector_alloc(n);
    d = vector_short_alloc(n + 1);
    propscore = vector_alloc(n);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    vector_free(propscore);
    // test 3
    y = vector_alloc(n);
    d = vector_short_alloc(n);
    propscore = vector_alloc(n + 1);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    vector_free(propscore);

    // Caliper check
    err_code = CMModelKnownPropscore_ERROR_CALIPER;
    n = 10;
    delta = -0.1;
    y = vector_alloc(n);
    d = vector_short_alloc(n);
    propscore = vector_alloc(n);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    // printf("ret_code = %d.\n", ret_code);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    vector_free(propscore);

    // Treatment vector check
    // test 1: is vector zero/one
    err_code = CMModelKnownPropscore_ERROR_TREATMENTZEROONE;
    n = 3;
    delta = 0.1;
    cm_model->n = n;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 2);
    propscore = vector_alloc(n);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_short_free(d);
    // test 2: balance of treated and control units
    err_code = CMModelKnownPropscore_ERROR_TREATMENTIMBALANCE;
    d = vector_short_calloc(n);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_short_free(d);
    // test 3: balance of treated and control units
    err_code = CMModelKnownPropscore_ERROR_TREATMENTIMBALANCE;
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    vector_short_set(d, 1, 1);
    vector_short_set(d, 2, 1);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    vector_free(propscore);

    // Propensity score
    // test 1
    err_code = CMModelKnownPropscore_ERROR_PROPSCOREOUTOFRANGE;
    n = 3;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    vector_set(propscore, 0, 2);
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(propscore);
    // test 2
    propscore = vector_alloc(n);
    vector_set(propscore, 0, 0);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(propscore);
    // test 3
    propscore = vector_alloc(n);
    vector_set(propscore, 0, 1);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(propscore);
    // test 4
    propscore = vector_alloc(n);
    vector_set(propscore, 0, -2);
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(propscore);

    // Variance estimation: bandwidth exponent check
    // test 1: beta too small
    err_code = CMModelKnownPropscore_ERROR_BETA;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = -1.0;
    alpha = 0;
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: beta too large
    err_code = CMModelKnownPropscore_ERROR_BETA;
    beta = 1.0;
    alpha = 0.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: truncation exponent check
    // test 1: alpha too small
    err_code = CMModelKnownPropscore_ERROR_ALPHA;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = 0;
    alpha = -1.0;
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: alpha too large
    err_code = CMModelKnownPropscore_ERROR_ALPHA;
    beta = 0.0;
    alpha = 1.0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: bandwidth versus truncation exponent check
    // test 1: alpha too large compared to beta!=0
    err_code = CMModelKnownPropscore_ERROR_BETAALPHA;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = 0.1;
    alpha = 0.2;   
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: no error for beta=0
    err_code = 0;
    beta = 0;
    alpha = 0.2;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 3: no error for alpha=0
    err_code = 0;
    beta = 0.1;
    alpha = 0;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: truncation sequence scale
    // test 1: kappa_a is negative
    err_code = CMModelKnownPropscore_ERROR_KAPPAA;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = 0.1;
    alpha = 0;
    kappa_a = -2.0;   
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: no error for kappa_a = 0
    err_code = 0;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = 0.1;
    alpha = 0; 
    kappa_a = 0.0;  
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);

    // Variance estimation: bandwidth sequence scale
    // test 1: kappa_a is negative
    err_code = CMModelKnownPropscore_ERROR_KAPPAGAMMA;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = 0.1;
    alpha = 0;
    kappa_a = 0.0;
    kappa_gamma = -3.14;   
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    // test 2: no error for kappa_a = 0
    err_code = 0;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = 0.1;
    alpha = 0; 
    kappa_a = 0.0;  
    kappa_gamma = 0.0;
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma = kappa_gamma;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);


    // All checks pass
    err_code = 0;
    n = 10;
    delta = 0.1;
    y = vector_alloc(n);
    d = vector_short_calloc(n);
    vector_short_set(d, 0, 1);
    propscore = vector_alloc(n);
    for (int i = 0; i < n; i++){
        vector_set(propscore, i, 0.1);
    }
    beta = 0.2;
    alpha = 0.1;
    kappa_a = 0.01;
    cm_model->n = n;
    cm_model->y = y;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_a = kappa_a;
    ret_code = cm_cmmodelknownpropscore_check_values(cm_model, test_mode);
    assert(ret_code == err_code);
    vector_free(y);
    vector_short_free(d);
    vector_free(propscore);
    printf("Tests for `cm_cmmodelknownpropscore_check_values` completed: no issues found.\n");
}


// Compute the single-index theta^T*x_i for i=0,...,n-1.
vector *cm_compute_singleindex_score(matrix *x, vector *theta){
    vector *score = vector_alloc(x->size1);
    double score_i; // single-index
    for (int i = 0; i < x->size1; i++){
        score_i = vector_get(theta, theta->size - 1); // intercept
        for (int j = 0; j < x->size2; j++){
            score_i += matrix_get(x, i, j) * vector_get(theta, j);
        }
        vector_set(score, i, score_i);
    }
    return score;
}

void test_cm_compute_singleindex_score(void){
    size_t n, k, i;
    n = 4;
    k = 3;
    matrix *x = matrix_alloc(n, k);
    vector *theta = vector_alloc(k+1);
    vector *score = vector_alloc(n);
    double score0, score1, score2, score3;
    // test 1
    vector_set(theta, 0, -7.1);
    vector_set(theta, 1, -2.3);
    vector_set(theta, 2, 5.1);
    vector_set(theta, 3, 0.001);
    matrix_set(x, 0, 0, -2.3);  matrix_set(x, 0, 1, 3.6);       matrix_set(x, 0, 2, 6.2);
    matrix_set(x, 1, 0, 2.3);   matrix_set(x, 1, 1, 50.6);      matrix_set(x, 1, 2, 1.5);
    matrix_set(x, 2, 0, 96.1);  matrix_set(x, 2, 1, 10.3);      matrix_set(x, 2, 2, 3.6);
    matrix_set(x, 3, 0, 97.1);  matrix_set(x, 3, 1, 20.78);     matrix_set(x, 3, 2, 4.7);
    i = 0;
    score0 = vector_get(theta, 3) + matrix_get(x, i, 0) * vector_get(theta, 0) +
                                    matrix_get(x, i, 1) * vector_get(theta, 1) +
                                    matrix_get(x, i, 2) * vector_get(theta, 2);
    i = 1;
    score1 = vector_get(theta, 3) + matrix_get(x, i, 0) * vector_get(theta, 0) +
                                    matrix_get(x, i, 1) * vector_get(theta, 1) +
                                    matrix_get(x, i, 2) * vector_get(theta, 2);
    i = 2;
    score2 = vector_get(theta, 3) + matrix_get(x, i, 0) * vector_get(theta, 0) +
                                    matrix_get(x, i, 1) * vector_get(theta, 1) +
                                    matrix_get(x, i, 2) * vector_get(theta, 2);
    i = 3;
    score3 = vector_get(theta, 3) + matrix_get(x, i, 0) * vector_get(theta, 0) +
                                    matrix_get(x, i, 1) * vector_get(theta, 1) +
                                    matrix_get(x, i, 2) * vector_get(theta, 2);                                                                                       
    vector_set(score, 0, score0);
    vector_set(score, 1, score1);
    vector_set(score, 2, score2);
    vector_set(score, 3, score3);
    vector *result = cm_compute_singleindex_score(x, theta);
    assert(vector_isequal(score, result));
    // cleanup
    matrix_free(x);
    vector_free(theta);
    printf("Tests for `cm_compute_singleindex_score` completed: no issues found.\n");
}

// Compute the propensity score for the whole sample
vector *cm_compute_propscore(vector *singleindex_score, int modeltype_index){
    vector *propscore = vector_alloc(singleindex_score->size);
    double (*g)(double);
    g = propscore_g[modeltype_index];
    for (int i = 0; i < singleindex_score->size; i++)
    {
        vector_set(propscore, i, g(vector_get(singleindex_score, i)));
    }
    return propscore;
}

void test_cm_compute_propscore(void){
    size_t n, i;
    n = 4;
    int modeltype_index = 0;
    double (*g)(double);
    g = propscore_g[modeltype_index];
    vector *singleindex_score = vector_alloc(n);
    vector *propscore = vector_alloc(n);
    // test 1
    vector_set(singleindex_score, 0, -2.3);
    vector_set(singleindex_score, 1, 2.3);
    vector_set(singleindex_score, 2, 96.1);
    vector_set(singleindex_score, 3, 97.1);                                                                                      
    vector_set(propscore, 0, g(vector_get(singleindex_score, 0)));
    vector_set(propscore, 1, g(vector_get(singleindex_score, 1)));
    vector_set(propscore, 2, g(vector_get(singleindex_score, 2)));
    vector_set(propscore, 3, g(vector_get(singleindex_score, 3)));
    vector *result = cm_compute_propscore(singleindex_score, modeltype_index);
    assert(vector_isequal(propscore, result));

    // test 2
    vector_set(singleindex_score, 0, 0.003);
    vector_set(singleindex_score, 1, -2.3);
    vector_set(singleindex_score, 2, 6.1);
    vector_set(singleindex_score, 3, 1.1);                                                                                      
    vector_set(propscore, 0, g(vector_get(singleindex_score, 0)));
    vector_set(propscore, 1, g(vector_get(singleindex_score, 1)));
    vector_set(propscore, 2, g(vector_get(singleindex_score, 2)));
    vector_set(propscore, 3, g(vector_get(singleindex_score, 3)));
    result = cm_compute_propscore(singleindex_score, modeltype_index);
    assert(vector_isequal(propscore, result));
    // cleanup
    vector_free(singleindex_score);
    vector_free(propscore);
    vector_free(result);
    printf("Tests for `cm_compute_propscore` completed: no issues found.\n");
}


// // Compute the propensity score for the whole sample
// vector *cm_compute_propscore(matrix *x, int modeltype_index, vector *theta){
//     vector *propscore = vector_alloc(x->size1);
//     double score_i;  // single-index
//     double (*g)(double);
//     g = propscore_g[modeltype_index];
//     for (int i = 0; i<x->size1; i++){
//         score_i = vector_get(theta, theta->size - 1);  // intercept
//         for (int j = 0; j<x->size2; j++){
//             score_i += matrix_get(x, i, j) * vector_get(theta, j);
//         }
//         vector_set(propscore, i, g(score_i));
//     }
//     return propscore;
// }


typedef struct CMThreadArgumentsMatches{
    long threadid;
    CMModelKnownPropscore *cm_model;
    size_t *pscore_sortindex;    // index array which sorts propensity score is ascending order
    size_t pscore_sortindex_low; // a thread processes observations in pscore_sortindex in range [pscore_sortindex_low,pscore_sortindex_high)
    size_t pscore_sortindex_high;
    dynamicarray_int **match_indices; // array of dynamic arrays of indicies of matched observations
    vector_int *matches;              // number of matches
} CMThreadArgumentsMatches;

void *cm_matches_thread(void *thread_arguments){
    CMThreadArgumentsMatches *thread_args = (CMThreadArgumentsMatches *)thread_arguments;
    double dist;
    size_t n = thread_args->cm_model->d->size;
    size_t nr_matches = 5 * cm_int_max((int) log(n), 1); // roughly n * delta matches for each unit
    // allocate new local variables if and only if NUM_THREAD>1; otherwise use those in thread_args directly
    vector_int *matches_local = (NUM_THREADS > 1 ? vector_int_calloc(n) : thread_args->matches);
    dynamicarray_int **match_indices_local = (NUM_THREADS > 1 ? malloc(n * sizeof(dynamicarray_int)) : thread_args->match_indices);
    // dynamicarray_int **match_indices_local = thread_args->match_indices;
    // loop through propscore starting with smallest one
    for (size_t i = thread_args->pscore_sortindex_low; i < thread_args->pscore_sortindex_high; i++){
        // allocate dynamic list which will hold indices of matches
        match_indices_local[thread_args->pscore_sortindex[i]] = dynamicarray_int_alloc((int) nr_matches);
        // look to the left until left-neighbourhood is exhausted
        for (int j = i; 0 <= j; j--){
            // only look at different treatment group
            if (vector_short_get(thread_args->cm_model->d, thread_args->pscore_sortindex[j]) != vector_short_get(thread_args->cm_model->d, thread_args->pscore_sortindex[i])){
                // distance between propensity scores
                dist = cm_abs_distance(vector_get(thread_args->cm_model->propscore, thread_args->pscore_sortindex[j]),
                                       vector_get(thread_args->cm_model->propscore, thread_args->pscore_sortindex[i]));
                if (dist <= thread_args->cm_model->delta){ // match them if within delta
                    // update match indices
                    // increase number of matches by one
                    dynamicarray_int_append(match_indices_local[thread_args->pscore_sortindex[i]], thread_args->pscore_sortindex[j]);
                    vector_int_set(matches_local, thread_args->pscore_sortindex[i], vector_int_get(matches_local, thread_args->pscore_sortindex[i]) + 1);
                } else { // there will be no matches to the left from here, neighbourhood is exhausted, break inner loop
                    break;
                }
            }
        }
        // look to the right until right-neighbourhood is exhausted
        for (int j = i; j < n; j++){
            // only look at different treatment group
            if (vector_short_get(thread_args->cm_model->d, thread_args->pscore_sortindex[j]) != vector_short_get(thread_args->cm_model->d, thread_args->pscore_sortindex[i])){
                // distance between propensity scores
                dist = cm_abs_distance(vector_get(thread_args->cm_model->propscore, thread_args->pscore_sortindex[j]),
                                       vector_get(thread_args->cm_model->propscore, thread_args->pscore_sortindex[i]));
                if (dist <= thread_args->cm_model->delta){ // match them if within delta
                    // update match indices
                    // increase number of matches by one
                    dynamicarray_int_append(match_indices_local[thread_args->pscore_sortindex[i]], thread_args->pscore_sortindex[j]);
                    vector_int_set(matches_local, thread_args->pscore_sortindex[i], vector_int_get(matches_local, thread_args->pscore_sortindex[i]) + 1);
                } else { // there will be no matches to the right from here, neighbourhood is exhausted, break inner loop
                    break;
                }
            }
        }
    }
    // merge in results
    if (NUM_THREADS > 1){
        pthread_mutex_lock(&mutex_lock);
        //  assign pointer for observations dealt by a given thread
        for (size_t i = thread_args->pscore_sortindex_low; i < thread_args->pscore_sortindex_high; i++){
            thread_args->match_indices[thread_args->pscore_sortindex[i]] = match_indices_local[thread_args->pscore_sortindex[i]];
        }
        // addition works because threads work on disjoint set of observation indices
        vector_int_add(thread_args->matches, matches_local);
        pthread_mutex_unlock(&mutex_lock);
    }
    void *ptr;
    return ptr;
}

/*
Computes the indices of matches for each observation.
The indices are written to `results->match_indices`, so that
`results->match_indices[i]` is a dynamic array containing the
indices of observations in the data set matched to observation i.

Also computes number of matches for each unit.
This vector is loaded into `results->number_of_matches`.
*/
void cm_matches(CMModelKnownPropscore *cm_model, CMResults *results){
    size_t n = cm_model->propscore->size;
    // to hold number of matches for each unit
    vector_int *matches = vector_int_calloc(n);
    // to hold indices of matches for every unit
    dynamicarray_int **match_indices = malloc(n * sizeof(dynamicarray_int));
    // sort propensity scores if not called internally
    size_t *pscore_sortindex = cm_model->_called_internally ? cm_model->_pscore_sortindex : vector_sort_index(cm_model->propscore);
    // traverse neighbourhood of propscores
    if (NUM_THREADS > 1){ // multithreaded version
        CMThreadArgumentsMatches thread_args[NUM_THREADS];
        int obs_per_thread = n / NUM_THREADS; // integer division gives floor(n/NUM_THREADS);
        // printf("Multithreading started, obs_per_thread = %d.\n", obs_per_thread);
        //  create threads
        int result_code;
        pthread_t threads[NUM_THREADS];
        for (long t = 0; t < NUM_THREADS; t++){
            if (CM_VERBOSE_MODE){
                printf("`cm_matches`: creating thread %ld.\n", t);
            }
            // arguments to threads
            thread_args[t].threadid = t;
            thread_args[t].cm_model = cm_model;
            thread_args[t].pscore_sortindex = pscore_sortindex;
            thread_args[t].pscore_sortindex_low = t * obs_per_thread;
            thread_args[t].pscore_sortindex_high = (t == NUM_THREADS - 1 ? n : (t + 1) * obs_per_thread); // take care of remaining observations
            thread_args[t].match_indices = match_indices;
            thread_args[t].matches = matches;
            result_code = pthread_create(&threads[t], NULL, cm_matches_thread, (void *)&thread_args[t]);
            if (CM_VERBOSE_MODE){
                printf("`cm_matches`: thread %ld pscore_sortindex_high is %zu.\n", t, thread_args[t].pscore_sortindex_high);
            }
            assert(!result_code);
        }
        // wait for threads to finish
        for (long t = 0; t < NUM_THREADS; t++){
            result_code = pthread_join(threads[t], NULL);
            assert(!result_code);
            if (CM_VERBOSE_MODE){
                printf("`cm_matches`: thread %ld has finished.\n", t);
            }
        }
    } else { // single thread version
        CMThreadArgumentsMatches *thread_args = malloc(sizeof(CMThreadArgumentsMatches));
        thread_args->cm_model = cm_model;
        thread_args->pscore_sortindex = pscore_sortindex;
        thread_args->pscore_sortindex_low = 0;
        thread_args->pscore_sortindex_high = n;
        thread_args->match_indices = match_indices;
        thread_args->matches = matches;
        cm_matches_thread((void *)thread_args);
    }
    results->match_indices = match_indices;
    results->number_of_matches = matches;
}

void test_cm_matches(void){
    double delta;
    CMModelKnownPropscore *cm_model = malloc(sizeof(CMModelKnownPropscore));
    CMResults *results = malloc(sizeof(CMResults));
    size_t n = 4;
    vector_short *d = vector_short_alloc(4);
    vector *propscore = vector_alloc(4);
    dynamicarray_int **match_indices = malloc(4 * sizeof(dynamicarray_int)); // true indices
    size_t capacity = 4;
    vector_int *number_of_matches = vector_int_alloc(4);
    // test 1
    // matchmatrix =
    //      c1 c2 t1 t2
    //      _  _  _  _
    // c1 | 0, 0, 1, 1,
    // c2 | 0, 0, 1, 1,
    // t1 | 1, 1, 0, 0,
    // t2 | 1, 1, 0, 0
    delta = 1.0;
    vector_short_set(d, 0, 0);
    vector_set(propscore, 0, 0.1); // c1: first control unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.5); // c2: second control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.9); // t1: first treated unit
    vector_short_set(d, 3, 1);
    vector_set(propscore, 3, 0.2); // t2: second treated unit
    vector_int_set(number_of_matches, 0, 2);
    vector_int_set(number_of_matches, 1, 2);
    vector_int_set(number_of_matches, 2, 2);
    vector_int_set(number_of_matches, 3, 2);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 2);
    dynamicarray_int_append(match_indices[0], 3); // match indices for c1
    dynamicarray_int_append(match_indices[1], 2);
    dynamicarray_int_append(match_indices[1], 3); // match indices for c2
    dynamicarray_int_append(match_indices[2], 0);
    dynamicarray_int_append(match_indices[2], 1); // match indices for t1
    dynamicarray_int_append(match_indices[3], 0);
    dynamicarray_int_append(match_indices[3], 1); // match indices for t2
    cm_model->n = n;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_matches(cm_model, results);
    for (int i = 0; i < 4; i++){
        // for (int j=0; j<results->match_indices[i]->usage; j++){
        //     printf("match_indices[%d][%d] = %d.\n", i, j, dynamicarray_int_get(match_indices[i], j));
        //     printf("results->match_indices[%d][%d] = %d.\n", i, j, dynamicarray_int_get(results->match_indices[i], j));
        // }
        assert(dynamicarray_int_permute_equal(results->match_indices[i], match_indices[i]));
    }
    assert(vector_int_equal(results->number_of_matches, number_of_matches));

    // test 2
    // matchmatrix =
    //      c1 c2 t1 t2
    //      _  _  _  _
    // c1 | 0, 0, 0, 1,
    // c2 | 0, 0, 0, 1,
    // t1 | 0, 0, 0, 0,
    // t2 | 1, 1, 0, 0
    delta = 0.3;
    vector_short_set(d, 0, 0);
    vector_set(propscore, 0, 0.1); // c1: first control unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.5); // c2: second control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.9); // t1: first treated unit
    vector_short_set(d, 3, 1);
    vector_set(propscore, 3, 0.2); // t2: second treated unit
    vector_int_set(number_of_matches, 0, 1);
    vector_int_set(number_of_matches, 1, 1);
    vector_int_set(number_of_matches, 2, 0);
    vector_int_set(number_of_matches, 3, 2);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 3); // match indices for c1
    dynamicarray_int_append(match_indices[1], 3); // match indices for c2
    // .    no match indices for t1
    dynamicarray_int_append(match_indices[3], 0);
    dynamicarray_int_append(match_indices[3], 1); // match indices for t2
    cm_model->n = n;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_matches(cm_model, results);
    for (int i = 0; i < 4; i++){
        if (i != 2){
            assert(dynamicarray_int_permute_equal(results->match_indices[i], match_indices[i]));
        } else {
            assert(dynamicarray_int_isempty(results->match_indices[i]));
        }
    }
    assert(vector_int_equal(results->number_of_matches, number_of_matches));

    // test 3
    // matchmatrix =
    //      t1 c1 t2 c2
    //      _  _  _  _
    // t1 | 0, 0, 0, 1,
    // c1 | 0, 0, 1, 0,
    // t2 | 0, 1, 0, 1,
    // c2 | 1, 0, 1, 0
    delta = 0.4;
    vector_short_set(d, 0, 1);
    vector_set(propscore, 0, 0.9); // t1: first treated unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.1); // c1: first control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.2); // t2: second treated unit
    vector_short_set(d, 3, 0);
    vector_set(propscore, 3, 0.5); // c2: second control unit
    vector_int_set(number_of_matches, 0, 1);
    vector_int_set(number_of_matches, 1, 1);
    vector_int_set(number_of_matches, 2, 2);
    vector_int_set(number_of_matches, 3, 2);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 3); // match indices for t1
    dynamicarray_int_append(match_indices[1], 2); // match indices for c1
    dynamicarray_int_append(match_indices[2], 1);
    dynamicarray_int_append(match_indices[2], 3); // match indices for t2
    dynamicarray_int_append(match_indices[3], 0);
    dynamicarray_int_append(match_indices[3], 2); // match indices for c2
    cm_model->n = n;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_matches(cm_model, results);
    for (int i = 0; i < 4; i++){
        assert(dynamicarray_int_permute_equal(results->match_indices[i], match_indices[i]));
    }
    assert(vector_int_equal(results->number_of_matches, number_of_matches));

    // test 4
    // matchmatrix =
    //      c1 c2 t1 c3
    //      _  _  _  _
    // c1 | 0, 0, 1, 0,
    // c2 | 0, 0, 0, 0,
    // t1 | 1, 0, 0, 0,
    // c3 | 0, 0, 0, 0
    delta = 0.1;
    vector_short_set(d, 0, 0);
    vector_set(propscore, 0, 0.9); // c1: first control unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.1); // c2: second control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.8); // t1: first treated unit
    vector_short_set(d, 3, 0);
    vector_set(propscore, 3, 0.5); // c3: thrid control unit
    vector_int_set(number_of_matches, 0, 1);
    vector_int_set(number_of_matches, 1, 0);
    vector_int_set(number_of_matches, 2, 1);
    vector_int_set(number_of_matches, 3, 0);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 2); // match indices for c1
    // . no match indices for c2
    dynamicarray_int_append(match_indices[2], 0); // match indices for t1
                                                  // .  no match indices for t2
    cm_model->n = n;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_matches(cm_model, results);
    assert(dynamicarray_int_permute_equal(results->match_indices[0], match_indices[0]) && dynamicarray_int_permute_equal(results->match_indices[2], match_indices[2]));
    assert(dynamicarray_int_isempty(results->match_indices[1]) && dynamicarray_int_isempty(results->match_indices[3]));
    assert(vector_int_equal(results->number_of_matches, number_of_matches));

    // test 5
    // matchmatrix =
    //      t1 c1 t2 t3
    //      _  _  _  _
    // t1 | 0, 0, 0, 0,
    // c1 | 0, 0, 0, 1,
    // t2 | 0, 0, 0, 0,
    // t3 | 0, 1, 0, 0
    delta = 0.1;
    vector_short_set(d, 0, 1);
    vector_set(propscore, 0, 0.9); // t1: first treated unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.2); // c1: first control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.8); // t2: second treated unit
    vector_short_set(d, 3, 1);
    vector_set(propscore, 3, 0.1); // t3: third treated unit
    vector_int_set(number_of_matches, 0, 0);
    vector_int_set(number_of_matches, 1, 1);
    vector_int_set(number_of_matches, 2, 0);
    vector_int_set(number_of_matches, 3, 1);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    // . no match indices for c1
    dynamicarray_int_append(match_indices[1], 3); // match indices for c2
    // . no match indices for t1
    dynamicarray_int_append(match_indices[3], 1); // match indices for t2
    cm_model->n = n;
    cm_model->d = d;
    cm_model->propscore = propscore;
    cm_model->delta = delta;
    cm_matches(cm_model, results);
    for (int i = 0; i < 4; i++){
        if ((i != 0) && (i != 2)){
            assert(dynamicarray_int_permute_equal(results->match_indices[i], match_indices[i]));
        } else{
            assert(dynamicarray_int_isempty(results->match_indices[i]));
        }
    }
    assert(vector_int_equal(results->number_of_matches, number_of_matches));

    // cleanup
    vector_short_free(d);
    vector_free(propscore);
    for (int i = 0; i < 4; i++){
        dynamicarray_int_free(match_indices[i]);
    }
    free(match_indices);
    vector_int_free(number_of_matches);
    printf("Tests for `cm_matches` completed: no issues found.\n");
}

typedef struct CMThreadArgumentsTE{
    long threadid;
    CMModelKnownPropscore *cm_model;
    CMResults *results;
    size_t n_treated;
    size_t index_low; // thread processes observations in range [index_low, index_high)
    size_t index_high;
    double *ate_hat;
    double *att_hat;
} CMThreadArgumentsTE;

void *cm_te_estimator_thread(void *thread_arguments){
    CMThreadArgumentsTE *thread_args = (CMThreadArgumentsTE *)thread_arguments;
    size_t n = thread_args->cm_model->d->size;
    size_t n_treated = thread_args->n_treated;
    double ate_hat_contribution, att_hat_contribution;
    ate_hat_contribution = att_hat_contribution = 0;
    double w_i;
    // compute contribution of observations in range [index_low, index_high) to estimators
    for (int i = thread_args->index_low; i < thread_args->index_high; i++){
        if (thread_args->results->match_indices[i]->usage){ // only consider observations with matches
            w_i = 0;
            // compute w_i = (sum of 1/M_j for j in the match set of i)
            for (int k = 0; k < thread_args->results->match_indices[i]->usage; k++){                                                                                                                                          // k-th element in i's match set
                w_i += 1.0 / vector_int_get(thread_args->results->number_of_matches, dynamicarray_int_get(thread_args->results->match_indices[i], k)); //  1/M_j
            }
            // compute estimators
            ate_hat_contribution += (2 * vector_short_get(thread_args->cm_model->d, i) - 1) * (1 + w_i) * vector_get(thread_args->cm_model->y, i) / n;
            att_hat_contribution += (vector_short_get(thread_args->cm_model->d, i) - (1 - vector_short_get(thread_args->cm_model->d, i)) * w_i) * vector_get(thread_args->cm_model->y, i) / n_treated;
        }
    }
    // merge in results
    pthread_mutex_lock(&mutex_lock);
    *(thread_args->ate_hat) += ate_hat_contribution;
    *(thread_args->att_hat) += att_hat_contribution;
    pthread_mutex_unlock(&mutex_lock);
    void *ptr;
    return ptr;
}

void test_cm_te_estimator_thread(void){
    // setup
    const size_t n = 4;
    size_t n_treated;
    double ate_hat, att_hat;
    double result_ate_hat, result_att_hat;
    CMModelKnownPropscore *cm_model = malloc(sizeof(CMModelKnownPropscore));
    CMResults *results = malloc(sizeof(CMResults));
    CMThreadArgumentsTE *thread_args = malloc(sizeof(CMThreadArgumentsTE));
    vector_short *d = vector_short_alloc(n);
    vector *y = vector_alloc(n);
    double *yy = malloc(n * sizeof(double));  // for less typing
    dynamicarray_int **match_indices = malloc(4 * sizeof(dynamicarray_int)); // true indices
    size_t capacity = 4;
    vector_int *number_of_matches = vector_int_alloc(4);
    
    // test 1
    // matchmatrix =
    //      c1 c2 t1 t2
    //      _  _  _  _
    // c1 | 0, 0, 1, 1,
    // c2 | 0, 0, 1, 1,
    // t1 | 1, 1, 0, 0,
    // t2 | 1, 1, 0, 0
    vector_short_set(d, 0, 0);  yy[0] = 1000.1;   vector_set(y, 0, yy[0]);
    vector_short_set(d, 1, 0);  yy[1] = -30.3;   vector_set(y, 1, yy[1]); 
    vector_short_set(d, 2, 1);  yy[2] = 40.1;    vector_set(y, 2, yy[2]);
    vector_short_set(d, 3, 1);  yy[3] = 0.5;     vector_set(y, 3, yy[3]);
    vector_int_set(number_of_matches, 0, 2);
    vector_int_set(number_of_matches, 1, 2);
    vector_int_set(number_of_matches, 2, 2);
    vector_int_set(number_of_matches, 3, 2);
    for (int i = 0; i < n; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 2);
    dynamicarray_int_append(match_indices[0], 3); // match indices for c1
    dynamicarray_int_append(match_indices[1], 2);
    dynamicarray_int_append(match_indices[1], 3); // match indices for c2
    dynamicarray_int_append(match_indices[2], 0);
    dynamicarray_int_append(match_indices[2], 1); // match indices for t1
    dynamicarray_int_append(match_indices[3], 0);
    dynamicarray_int_append(match_indices[3], 1); // match indices for t2
    // true estimates
    n_treated = 2;
    ate_hat = (1.0 / n) * ((1.0 / 2) * (yy[2] + yy[3]) - yy[0] +
                           (1.0 / 2) * (yy[2] + yy[3]) - yy[1] +
                           yy[2] - (1.0 / 2) * (yy[0] + yy[1]) +
                           yy[3] - (1.0 / 2) * (yy[0] + yy[1]));
    att_hat = (1.0 / n_treated) * (yy[2] - (1.0 / 2) * (yy[0] + yy[1]) +
                                   yy[3] - (1.0 / 2) * (yy[0] + yy[1]));;  
    // estimates from function
    result_ate_hat = result_att_hat = 0;
    cm_model->d = d;
    cm_model->y = y;
    results->match_indices = match_indices;
    results->number_of_matches = number_of_matches;
    thread_args->cm_model = cm_model;
    thread_args->results = results;
    thread_args->index_low = 0;
    thread_args->index_high = n;
    thread_args->n_treated = n_treated;
    thread_args->ate_hat = &result_ate_hat;
    thread_args->att_hat = &result_att_hat;
    cm_te_estimator_thread((void *) thread_args);
    results->ate_hat = result_ate_hat;
    results->att_hat = result_att_hat;
    // printf("ate_hat = %f, results->ate_hat = %f\n", ate_hat, results->ate_hat);
    // printf("att_hat = %f, results->att_hat = %f\n", att_hat, results->att_hat);
    assert(fabs(results->ate_hat - ate_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(results->att_hat - att_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));

    // test 2
    // matchmatrix =
    //      c1 c2 t1 t2
    //      _  _  _  _
    // c1 | 0, 0, 0, 1,
    // c2 | 0, 0, 0, 1,
    // t1 | 0, 0, 0, 0,
    // t2 | 1, 1, 0, 0
    vector_short_set(d, 0, 0);
    vector_short_set(d, 1, 0);
    vector_short_set(d, 2, 1);
    vector_short_set(d, 3, 1);
    vector_int_set(number_of_matches, 0, 1);
    vector_int_set(number_of_matches, 1, 1);
    vector_int_set(number_of_matches, 2, 0);
    vector_int_set(number_of_matches, 3, 2);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 3); // match indices for c1
    dynamicarray_int_append(match_indices[1], 3); // match indices for c2
    // .    no match indices for t1
    dynamicarray_int_append(match_indices[3], 0);
    dynamicarray_int_append(match_indices[3], 1); // match indices for t2
    // true estimates
    n_treated = 2;
    ate_hat = (1.0 / n) * ((1.0 / 1) * (yy[3]) - yy[0] +
                           (1.0 / 1) * (yy[3]) - yy[1] +
                           yy[3] - (1.0 / 2) * (yy[0] + yy[1]));
    att_hat = (1.0 / n_treated) * (yy[3] - (1.0 / 2) * (yy[0] + yy[1]));;  
    // estimates from function
    result_ate_hat = result_att_hat = 0;
    cm_model->d = d;
    cm_model->y = y;
    results->match_indices = match_indices;
    results->number_of_matches = number_of_matches;
    thread_args->cm_model = cm_model;
    thread_args->results = results;
    thread_args->index_low = 0;
    thread_args->index_high = n;
    thread_args->n_treated = n_treated;
    thread_args->ate_hat = &result_ate_hat;
    thread_args->att_hat = &result_att_hat;
    cm_te_estimator_thread((void *) thread_args);
    results->ate_hat = result_ate_hat;
    results->att_hat = result_att_hat;
    // printf("ate_hat = %f, results->ate_hat = %f\n", ate_hat, results->ate_hat);
    // printf("att_hat = %f, results->att_hat = %f\n", att_hat, results->att_hat);
    assert(fabs(results->ate_hat - ate_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(results->att_hat - att_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));

    // test 3
    // matchmatrix =
    //      t1 c1 t2 c2
    //      _  _  _  _
    // t1 | 0, 0, 0, 1,
    // c1 | 0, 0, 1, 0,
    // t2 | 0, 1, 0, 1,
    // c2 | 1, 0, 1, 0
    vector_short_set(d, 0, 1);
    vector_short_set(d, 1, 0);
    vector_short_set(d, 2, 1);
    vector_short_set(d, 3, 0);
    vector_int_set(number_of_matches, 0, 1);
    vector_int_set(number_of_matches, 1, 1);
    vector_int_set(number_of_matches, 2, 2);
    vector_int_set(number_of_matches, 3, 2);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 3); // match indices for t1
    dynamicarray_int_append(match_indices[1], 2); // match indices for c1
    dynamicarray_int_append(match_indices[2], 1);
    dynamicarray_int_append(match_indices[2], 3); // match indices for t2
    dynamicarray_int_append(match_indices[3], 0);
    dynamicarray_int_append(match_indices[3], 2); // match indices for c2
    // true estimates
    n_treated = 2;
    ate_hat = (1.0 / n) * (yy[0] - (1.0 / 1) * (yy[3]) +
                           (1.0 / 1) * (yy[2]) - yy[1] +
                           yy[2] - (1.0 / 2) * (yy[1] + yy[3]) +
                           (1.0 / 2) * (yy[0] + yy[2]) - yy[3]);
    att_hat = (1.0 / n_treated) * (yy[0] - (1.0 / 1) * (yy[3]) +
                                   yy[2] - (1.0 / 2) * (yy[1] + yy[3]));
    // estimates from function
    result_ate_hat = result_att_hat = 0;
    cm_model->d = d;
    cm_model->y = y;
    results->match_indices = match_indices;
    results->number_of_matches = number_of_matches;
    thread_args->cm_model = cm_model;
    thread_args->results = results;
    thread_args->index_low = 0;
    thread_args->index_high = n;
    thread_args->n_treated = n_treated;
    thread_args->ate_hat = &result_ate_hat;
    thread_args->att_hat = &result_att_hat;
    cm_te_estimator_thread((void *) thread_args);
    results->ate_hat = result_ate_hat;
    results->att_hat = result_att_hat;
    // printf("ate_hat = %f, results->ate_hat = %f\n", ate_hat, results->ate_hat);
    // printf("att_hat = %f, results->att_hat = %f\n", att_hat, results->att_hat);
    assert(fabs(results->ate_hat - ate_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(results->att_hat - att_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));    
    
    // test 4
    // matchmatrix =
    //      c1 c2 t1 c3
    //      _  _  _  _
    // c1 | 0, 0, 1, 0,
    // c2 | 0, 0, 0, 0,
    // t1 | 1, 0, 0, 0,
    // c3 | 0, 0, 0, 0
    vector_short_set(d, 0, 0);
    vector_short_set(d, 1, 0);
    vector_short_set(d, 2, 1);
    vector_short_set(d, 3, 0);
    vector_int_set(number_of_matches, 0, 1);
    vector_int_set(number_of_matches, 1, 0);
    vector_int_set(number_of_matches, 2, 1);
    vector_int_set(number_of_matches, 3, 0);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    dynamicarray_int_append(match_indices[0], 2); // match indices for c1
    // . no match indices for c2
    dynamicarray_int_append(match_indices[2], 0); // match indices for t1
    // .  no match indices for t2    
    // true estimates
    n_treated = 1;
    ate_hat = (1.0 / n) * ((1.0 / 1) * (yy[2]) - yy[0] +
                           yy[2] - (1.0 / 1) * (yy[0]));
    att_hat = (1.0 / n_treated) * (yy[2] - (1.0 / 1) * (yy[0]));
    // estimates from function
    result_ate_hat = result_att_hat = 0;
    cm_model->d = d;
    cm_model->y = y;
    results->match_indices = match_indices;
    results->number_of_matches = number_of_matches;
    thread_args->cm_model = cm_model;
    thread_args->results = results;
    thread_args->index_low = 0;
    thread_args->index_high = n;
    thread_args->n_treated = n_treated;
    thread_args->ate_hat = &result_ate_hat;
    thread_args->att_hat = &result_att_hat;
    cm_te_estimator_thread((void *) thread_args);
    results->ate_hat = result_ate_hat;
    results->att_hat = result_att_hat;
    // printf("ate_hat = %f, results->ate_hat = %f\n", ate_hat, results->ate_hat);
    // printf("att_hat = %f, results->att_hat = %f\n", att_hat, results->att_hat);
    assert(fabs(results->ate_hat - ate_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(results->att_hat - att_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));    


    // test 5
    // matchmatrix =
    //      t1 c1 t2 t3
    //      _  _  _  _
    // t1 | 0, 0, 0, 0,
    // c1 | 0, 0, 0, 1,
    // t2 | 0, 0, 0, 0,
    // t3 | 0, 1, 0, 0
    vector_short_set(d, 0, 1);
    vector_short_set(d, 1, 0);
    vector_short_set(d, 2, 1);
    vector_short_set(d, 3, 1);
    vector_int_set(number_of_matches, 0, 0);
    vector_int_set(number_of_matches, 1, 1);
    vector_int_set(number_of_matches, 2, 0);
    vector_int_set(number_of_matches, 3, 1);
    for (int i = 0; i < 4; i++){
        match_indices[i] = dynamicarray_int_alloc(capacity);
    }
    // . no match indices for c1
    dynamicarray_int_append(match_indices[1], 3); // match indices for c2
    // . no match indices for t1
    dynamicarray_int_append(match_indices[3], 1); // match indices for t2
    // true estimates
    n_treated = 3;
    ate_hat = (1.0 / n) * ((1.0 / 1) * (yy[3]) - yy[1] +
                           yy[3] - (1.0 / 1) * (yy[1]));
    att_hat = (1.0 / n_treated) * (yy[3] - (1.0 / 1) * (yy[1]));
    // estimates from function
    result_ate_hat = result_att_hat = 0;
    cm_model->d = d;
    cm_model->y = y;
    results->match_indices = match_indices;
    results->number_of_matches = number_of_matches;
    thread_args->cm_model = cm_model;
    thread_args->results = results;
    thread_args->index_low = 0;
    thread_args->index_high = n;
    thread_args->n_treated = n_treated;
    thread_args->ate_hat = &result_ate_hat;
    thread_args->att_hat = &result_att_hat;
    cm_te_estimator_thread((void *) thread_args);
    results->ate_hat = result_ate_hat;
    results->att_hat = result_att_hat;
    // printf("ate_hat = %f, results->ate_hat = %f\n", ate_hat, results->ate_hat);
    // printf("att_hat = %f, results->att_hat = %f\n", att_hat, results->att_hat);
    assert(fabs(results->ate_hat - ate_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(results->att_hat - att_hat) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER)); 
    // cleanup
    for (int i=0; i<4; i++){
        dynamicarray_int_free(match_indices[i]);
    }
    free(match_indices);
    vector_free(y);
    vector_short_free(d); 
    vector_int_free(number_of_matches);
    free(yy);
    free(results);
    free(cm_model);
    free(thread_args);
    printf("Tests for `cm_te_estimator_thread` completed: no issues found.\n");
}

// Estimates ATE and ATT, using matches in `results` and data in `cm_model`, and writes them to `results`.
void cm_te_estimator(CMModelKnownPropscore *cm_model, CMResults *results){
    size_t n = cm_model->d->size;
    size_t n_treated = vector_short_sum(cm_model->d);
    double ate_hat, att_hat;
    ate_hat = att_hat = 0;
    if (NUM_THREADS > 1){ // multithreading
        pthread_t threads[NUM_THREADS];
        CMThreadArgumentsTE thread_args[NUM_THREADS];
        int obs_per_thread = n / NUM_THREADS;
        int result_code;
        // creating threads
        for (long t = 0; t < NUM_THREADS; t++){
            if (CM_VERBOSE_MODE){
                printf("`cm_te_estimator`: creating thread %ld.\n", t);
            }
            thread_args[t].cm_model = cm_model;
            thread_args[t].results = results;
            thread_args[t].n_treated = n_treated;
            thread_args[t].index_low = t * obs_per_thread;
            thread_args[t].index_high = (t == NUM_THREADS - 1 ? n : (t + 1) * obs_per_thread);
            thread_args[t].ate_hat = &ate_hat;
            thread_args[t].att_hat = &att_hat;
            if (CM_VERBOSE_MODE){
                printf("`cm_te_estimator`: thread %ld index_high is %zu.\n", t, thread_args[t].index_high);
            }
            result_code = pthread_create(&threads[t], NULL, cm_te_estimator_thread, (void *)&thread_args[t]);
            assert(!result_code);
        }
        // wait for threads to finish
        for (long t = 0; t < NUM_THREADS; t++){
            result_code = pthread_join(threads[t], NULL);
            assert(!result_code);
            if (CM_VERBOSE_MODE){
                printf("`cm_te_estimator`: thread %ld has finished.\n", t);
            }
        }
    } else { // single thread
        CMThreadArgumentsTE *thread_args = malloc(sizeof(CMThreadArgumentsTE));
        thread_args->cm_model = cm_model;
        thread_args->results = results;
        thread_args->n_treated = n_treated;
        thread_args->index_low = 0;
        thread_args->index_high = cm_model->d->size;
        thread_args->ate_hat = &ate_hat;
        thread_args->att_hat = &att_hat;
        cm_te_estimator_thread((void *)thread_args);
    }
    results->ate_hat = ate_hat;
    results->att_hat = att_hat;
}


// Returns one if unit `i` is treated, zero otherwise.
int cm_kernel_filter_func(size_t i, void *filter_args){
    vector_short *d = (vector_short *)filter_args;
    return vector_short_get(d, i);
}


// Compute third variance components (q_d1-q_d0)^T * inv(fisher_info) * (q_d1-q_d0) and
// (q_d1_treated-q_d0_treated)^T * inv(fisher_info) * (q_d1_treated-q_d0_treated)
// in-place: it changes q_d1, q_d1_treated and fisher_info.
// The factor (1/p1_hat) is assumed to be included in q_d0_treated and q_d1_treated.
double *cm_compute_var_estpi(vector *q_d0, vector *q_d0_treated, vector *q_d1, vector *q_d1_treated, matrix *fisher_info){
    double *v_estpi = calloc(2, sizeof(double));
    // in-place difference q_d1-q_d0
    vector_sub(q_d1, q_d0);
    vector_sub(q_d1_treated, q_d0_treated);
    // invert fisher_info in-place with Cholesky decomposition
    linalg_cholesky_decomp(fisher_info);
    linalg_cholesky_invert(fisher_info);
    // compute quadratic form
    double *work = calloc(q_d0->size, sizeof(double));  // temporary array
    double *work_treated = calloc(q_d0->size, sizeof(double));  // temporary array
    // inv(fisher_info)*q
    for (int j = 0; j < q_d0->size; j++){
        for (int l = 0; l < q_d0->size; l++){
            work[j] += matrix_get(fisher_info, j, l) * vector_get(q_d1, l);  // fisher_info is already inverted in-place
            work_treated[j] += matrix_get(fisher_info, j, l) * vector_get(q_d1_treated, l);
        }
        // q^T*inv(fisher_info)*q
        v_estpi[0] += vector_get(q_d1, j) * work[j];
        v_estpi[1] += vector_get(q_d1_treated, j) * work_treated[j];
    }
    // cleanup
    free(work);
    free(work_treated);
    return v_estpi;
}

void test_cm_compute_var_estpi(void){
    // test 1: identity fisher info
    int k = 2;
    double *v_estpi = malloc(2 * sizeof(double));
    vector *q_d0 = vector_calloc(k);
    vector *q_d1 = vector_calloc(k);
    vector *q_d0_treated = vector_calloc(k);
    vector *q_d1_treated = vector_calloc(k);
    matrix *fisher_info = matrix_alloc(k, k);
    vector_set(q_d1, 0, 1);
    vector_set(q_d1, 1, 1);
    vector_set(q_d1_treated, 0, 1);
    vector_set(q_d1_treated, 1, 1);
    matrix_set(fisher_info, 0, 0, 1);
    matrix_set(fisher_info, 0, 1, 0);
    matrix_set(fisher_info, 1, 0, 0);
    matrix_set(fisher_info, 1, 1, 1);
    v_estpi[0] = 2; // (1, 1) * inv(fisher) * (1, 1)^T
    v_estpi[1] = 2; // (1, 1) * inv(fisher) * (1, 1)^T
    double *result = cm_compute_var_estpi(q_d0, q_d0_treated, q_d1, q_d1_treated, fisher_info);
    // printf("result[0] = %f, result[1] = %f\n", result[0], result[1]);
    assert(fabs(result[0] - v_estpi[0]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(result[1] - v_estpi[1]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    // test 2: diagonal fisher info
    matrix_set(fisher_info, 0, 0, 0.5);
    matrix_set(fisher_info, 0, 1, 0);
    matrix_set(fisher_info, 1, 0, 0);
    matrix_set(fisher_info, 1, 1, 0.5);
    v_estpi[0] = 4; // (1, 1) * inv(fisher) * (1, 1)^T
    v_estpi[1] = 4; // (1, 1) * inv(fisher) * (1, 1)^T
    result = cm_compute_var_estpi(q_d0, q_d0_treated, q_d1, q_d1_treated, fisher_info);
    assert(fabs(result[0] - v_estpi[0]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(result[1] - v_estpi[1]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    // test 3: fisher_info =
    //  ---------
    // | 2  1
    // | 1  2
    //  ---------
    //
    // inv(fisher_info) =
    //  ---------
    // | 2/3    -1/3
    // | -1/3   2/3
    //  ---------
    matrix_set(fisher_info, 0, 0, 2);
    matrix_set(fisher_info, 0, 1, 1);
    matrix_set(fisher_info, 1, 0, 1);
    matrix_set(fisher_info, 1, 1, 2);
    v_estpi[0] = (1 * 2.0 / 3 + 1 * (-1.0 / 3)) + (1 * (-1.0 / 3) + 1 * 2.0 / 3); // (1, 1) * inv(fisher) * (1, 1)^T
    v_estpi[1] = (1 * 2.0 / 3 + 1 * (-1.0 / 3)) + (1 * (-1.0 / 3) + 1 * 2.0 / 3); // (1, 1) * inv(fisher) * (1, 1)^T
    result = cm_compute_var_estpi(q_d0, q_d0_treated, q_d1, q_d1_treated, fisher_info);
    assert(fabs(result[0] - v_estpi[0]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(result[1] - v_estpi[1]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    // test 4: different q's, fisher_info =
    //  ---------
    // | 2  1
    // | 1  2
    //  ---------
    //
    // inv(fisher_info) =
    //  ---------
    // | 2/3    -1/3
    // | -1/3   2/3
    //  ---------
    vector_set(q_d1, 0, 1);
    vector_set(q_d1, 1, 1);
    vector_set(q_d1_treated, 0, 1);
    vector_set(q_d1_treated, 1, 0);
    matrix_set(fisher_info, 0, 0, 2);
    matrix_set(fisher_info, 0, 1, 1);
    matrix_set(fisher_info, 1, 0, 1);
    matrix_set(fisher_info, 1, 1, 2);
    v_estpi[0] = (1 * 2.0 / 3 + 1 * (-1.0 / 3)) + (1 * (-1.0 / 3) + 1 * 2.0 / 3); // (1, 1) * inv(fisher) * (1, 1)^T
    v_estpi[1] = 2.0 / 3;                                                         // (1, 0) * inv(fisher) * (1, 0)^T
    result = cm_compute_var_estpi(q_d0, q_d0_treated, q_d1, q_d1_treated, fisher_info);
    assert(fabs(result[0] - v_estpi[0]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    assert(fabs(result[1] - v_estpi[1]) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    printf("Tests for `cm_compute_var_estpi` completed: no issues found.\n");
}

// Update Fisher-information matrix: adds i-th term to `fisher_info`.
int cm_update_fisher_info(size_t i, matrix *fisher_info, CMModelKnownPropscore *cm_model){
    // link function and derivative for appropriate propscore model
    double (*g)(double);
    double (*g_deriv)(double);
    g = propscore_g[cm_model->_modeltype_index];
    g_deriv = propscore_gderiv[cm_model->_modeltype_index];
    // score (single-index, not log-likelihood derivative) for observation i
    // double score = vector_get(cm_model->_theta, cm_model->_theta->size - 1); // intercept
    // for (int j = 0; j < cm_model->_x->size2; j++){
    //     score += matrix_get(cm_model->_x, i, j) * vector_get(cm_model->_theta, j);
    // }
    double score = vector_get(cm_model->_singleindex_score , i);
    double weight = pow(g_deriv(score), 2) / ((cm_model->propscore->size) * vector_get(cm_model->propscore, i) * (1 - vector_get(cm_model->propscore, i)));
    // update fisher information
    double update, x_j, x_l;
    for (int j = 0; j < fisher_info->size2; j++){
        for (int l = 0; l < fisher_info->size2; l++){
            // adjust for intercept
            x_j = j < fisher_info->size2 - 1 ? matrix_get(cm_model->_x, i, j) : 1.0;
            x_l = l < fisher_info->size2 - 1 ? matrix_get(cm_model->_x, i, l) : 1.0;
            update = weight * x_j * x_l;
            matrix_set(fisher_info, j, l, matrix_get(fisher_info, j, l) + update);
        }
    }
    return 0;
}


void test_cm_update_fisher_info(void){
    // setup
    size_t n, i, k, l0, l1;
    n = 4;
    k = 2;
    double score, update, weight;
    matrix *x = matrix_alloc(n, k);
    vector *theta = vector_alloc(k+1);
    vector *singleindex_score = vector_alloc(n);
    vector *propscore = vector_alloc(n);
    matrix *fisher_info = matrix_calloc(k+1, k+1);
    matrix *result = matrix_calloc(k+1, k+1);
    CMModelKnownPropscore *cm_model = malloc(sizeof(CMModelKnownPropscore));
    double (*g)(double);
    double (*g_deriv)(double);
    int modeltype_index = 0;  //
    g = propscore_g[modeltype_index];
    g_deriv = propscore_gderiv[modeltype_index];
    vector_set(theta, 0, 0.5);
    vector_set(theta, 1, 0.22);
    vector_set(theta, 2, 0.06);
    matrix_set(x, 0, 0, 5.0);   matrix_set(x, 0, 1, 2.0);
    matrix_set(x, 1, 0, -3.3);  matrix_set(x, 1, 1, -1.43);
    matrix_set(x, 2, 0, -10.1); matrix_set(x, 2, 1, 1.96);
    matrix_set(x, 3, 0, 1.3);   matrix_set(x, 3, 1, 2.8);
    // score and propscore for observations
    for (int j=0; j<4; j++){
        // score for obervation j
        score = vector_get(theta, k) + vector_get(theta, 0) * matrix_get(x, j, 0) + vector_get(theta, 1) * matrix_get(x, j, 1);
        vector_set(singleindex_score, j, score);
        vector_set(propscore, j, g(score));
    }
    cm_model->_x = x;
    cm_model->_modeltype_index = modeltype_index;
    cm_model->_singleindex_score = singleindex_score;
    cm_model->propscore = propscore;
    
    // test 0: i=0
    i = 0;
    weight = pow(g_deriv(vector_get(singleindex_score, i)), 2) / (n * vector_get(propscore, i) * (1 - vector_get(propscore, i)));
    // --- (0, 0)
    l0 = 0; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 1)
    l0 = 0; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 2)
    l0 = 0; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 0)
    l0 = 1; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 1)
    l0 = 1; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 2)
    l0 = 1; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 0)
    l0 = 2; l1 = 0;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 1)
    l0 = 2; l1 = 1;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 2)
    l0 = 2; l1 = 2;
    update = weight * 1.0 * 1.0;  // intercept-intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);    
    cm_update_fisher_info(i, result, cm_model);
    // for (int q0 =0; q0<k+1; q0++){
    //     for (int q1=0; q1<k+1; q1++){
    //         printf("fisher_info[%d, %d] = %f, result[%d, %d] = %f\n", q0, q1, matrix_get(fisher_info, q0, q1), q0, q1, matrix_get(result, q0, q1));
    //     }
    // }
    assert(matrix_isequal(result, fisher_info));

    // test 1: i=1
    i = 1;
    weight = pow(g_deriv(vector_get(singleindex_score, i)), 2) / (n * vector_get(propscore, i) * (1 - vector_get(propscore, i)));
    // --- (0, 0)
    l0 = 0; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 1)
    l0 = 0; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 2)
    l0 = 0; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 0)
    l0 = 1; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 1)
    l0 = 1; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 2)
    l0 = 1; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 0)
    l0 = 2; l1 = 0;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 1)
    l0 = 2; l1 = 1;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 2)
    l0 = 2; l1 = 2;
    update = weight * 1.0 * 1.0;  // intercept-intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);    
    cm_update_fisher_info(i, result, cm_model);
    // for (int q0 =0; q0<k+1; q0++){
    //     for (int q1=0; q1<k+1; q1++){
    //         printf("fisher_info[%d, %d] = %f, result[%d, %d] = %f\n", q0, q1, matrix_get(fisher_info, q0, q1), q0, q1, matrix_get(result, q0, q1));
    //     }
    // }
    assert(matrix_isequal(result, fisher_info));

    // test 2: i=2
    i = 2;
    weight = pow(g_deriv(vector_get(singleindex_score, i)), 2) / (n * vector_get(propscore, i) * (1 - vector_get(propscore, i)));
    // --- (0, 0)
    l0 = 0; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 1)
    l0 = 0; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 2)
    l0 = 0; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 0)
    l0 = 1; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 1)
    l0 = 1; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 2)
    l0 = 1; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 0)
    l0 = 2; l1 = 0;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 1)
    l0 = 2; l1 = 1;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 2)
    l0 = 2; l1 = 2;
    update = weight * 1.0 * 1.0;  // intercept-intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);    
    cm_update_fisher_info(i, result, cm_model);
    // for (int q0 =0; q0<k+1; q0++){
    //     for (int q1=0; q1<k+1; q1++){
    //         printf("fisher_info[%d, %d] = %f, result[%d, %d] = %f\n", q0, q1, matrix_get(fisher_info, q0, q1), q0, q1, matrix_get(result, q0, q1));
    //     }
    // }
    assert(matrix_isequal(result, fisher_info));

    // test 3: i=3
    i = 3;
    weight = pow(g_deriv(vector_get(singleindex_score, i)), 2) / (n * vector_get(propscore, i) * (1 - vector_get(propscore, i)));
    // --- (0, 0)
    l0 = 0; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 1)
    l0 = 0; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (0, 2)
    l0 = 0; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 0)
    l0 = 1; l1 = 0;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 1)
    l0 = 1; l1 = 1;
    update = weight * matrix_get(x, i, l0) * matrix_get(x, i, l1);
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (1, 2)
    l0 = 1; l1 = 2;
    update = weight * matrix_get(x, i, l0) * 1.0;  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 0)
    l0 = 2; l1 = 0;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 1)
    l0 = 2; l1 = 1;
    update = weight * 1.0 * matrix_get(x, i, l1);  // intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);
    // --- (2, 2)
    l0 = 2; l1 = 2;
    update = weight * 1.0 * 1.0;  // intercept-intercept
    matrix_set(fisher_info, l0, l1, matrix_get(fisher_info, l0, l1) + update);    
    cm_update_fisher_info(i, result, cm_model);
    // for (int q0 =0; q0<k+1; q0++){
    //     for (int q1=0; q1<k+1; q1++){
    //         printf("fisher_info[%d, %d] = %f, result[%d, %d] = %f\n", q0, q1, matrix_get(fisher_info, q0, q1), q0, q1, matrix_get(result, q0, q1));
    //     }
    // }
    assert(matrix_isequal(result, fisher_info));

    printf("Tests for `cm_update_fisher_info` completed: no issues found.\n");
}

// Update q's in the third variance component: adds i-th term to `q_d0, q_d0_treated, q_d1, q_d1_treated`.
// On output, the factor (1/p1_hat) is included in q_d0_treated and q_d1_treated.
int cm_update_qs(size_t i, vector *q_d0, vector *q_d0_treated, vector *q_d1, vector *q_d1_treated, CMModelKnownPropscore *cm_model, double p1_hat, double gamma_n, double gamma_derivative_n, double *h_d, double *qmu_d){
    // link function and its derivative for appropriate propscore model
    size_t n = cm_model->d->size;
    double (*g)(double);
    double (*g_deriv)(double);
    double (*g_inv_deriv)(double);
    g = propscore_g[cm_model->_modeltype_index];
    g_deriv = propscore_gderiv[cm_model->_modeltype_index];
    g_inv_deriv = propscore_ginvderiv[cm_model->_modeltype_index];
    // score (=single-index) for observation i
    // double score = vector_get(cm_model->_theta, cm_model->_theta->size - 1); // intercept
    // for (int j = 0; j < cm_model->_x->size2; j++){
    //     score += matrix_get(cm_model->_x, i, j) * vector_get(cm_model->_theta, j);
    // }
    double score = vector_get(cm_model->_singleindex_score, i);
    // Derivates w.r.t. to p:  (1/n) * (deriv mu_d / deriv p)(theta, g(theta^T*x_i))*g'(theta^T*x_i)
    // ---- (deriv qmu_d / deriv p):  (-1/(n*gamma_n^2))sum_{j:D_j=d} Y_j * K'((propscore_j-propscore_i)/gamma_n)
    double *qmu_d_deriv_p = nonpara_qk_deriv_eval_filter(vector_get(cm_model->propscore, i),
                                                         cm_model->propscore->data,
                                                         cm_model->propscore->size,
                                                         cm_model->propscore->stride,
                                                         cm_model->y->data,
                                                         cm_model->y->size,
                                                         cm_model->y->stride,
                                                         1, // exponent
                                                         nonpara_kernelderiv_gauss,
                                                         gamma_n,
                                                         1, // scale
                                                         0,
                                                         (void *)cm_model->d,
                                                         cm_kernel_filter_func);
    // ---- (deriv h_d / deriv p):  (-1/(n*gamma_n^2))sum_{j:D_j=d} K'((propscore_j-propscore_i)/gamma_n)
    double *h_d_deriv_p = nonpara_h_deriv_eval_filter(vector_get(cm_model->propscore, i),
                                                      cm_model->propscore->data,
                                                      cm_model->propscore->size,
                                                      cm_model->propscore->stride,
                                                      nonpara_kernelderiv_gauss,
                                                      gamma_n,
                                                      1, // scale
                                                      0,
                                                      (void *)cm_model->d,
                                                      cm_kernel_filter_func);
    // ---- (1/n) * (deriv mu_d / deriv p)
    double *mu_d_deriv_p = malloc(2 * sizeof(double));
    mu_d_deriv_p[0] = (qmu_d_deriv_p[0] * h_d[0] - qmu_d[0] * h_d_deriv_p[0]) / (n * pow(h_d[0], 2));
    mu_d_deriv_p[1] = (qmu_d_deriv_p[1] * h_d[1] - qmu_d[1] * h_d_deriv_p[1]) / (n * pow(h_d[1], 2));
    // ---- altogether
    double *deriv_p = malloc(2 * sizeof(double));
    deriv_p[0] = mu_d_deriv_p[0] * g_deriv(score);
    deriv_p[1] = mu_d_deriv_p[1] * g_deriv(score);
    // printf("h_d_deriv_p[0] = %f; h_d_deriv_p[1] = %f.\n", h_d_deriv_p[0], h_d_deriv_p[1]);
    // printf("qmu_d_deriv_p[0] = %f; qmu_d_deriv_p[1] = %f.\n", qmu_d_deriv_p[0], qmu_d_deriv_p[1]);
    // printf("mu_d_deriv_p[0] = %f; mu_d_deriv_p[1] = %f.\n", mu_d_deriv_p[0], mu_d_deriv_p[1]);
    // printf("deriv_p[0] = %f; deriv_p[1] = %f.\n", deriv_p[0], deriv_p[1]);


    // Derivates w.r.t. theta_j: (1/n) * ((deriv mu_d / deriv theta_j))_{j=0,...,k+1} ; update q's.
    double *qmu_d_deriv_thetaj = malloc(2 * sizeof(double));
    double *h_d_deriv_thetaj = malloc(2 * sizeof(double));
    double *mu_d_deriv_thetaj = malloc(2 * sizeof(double));
    double update_q_d0, update_q_d0_treated, update_q_d1, update_q_d1_treated, x_ij;
    for (int j = 0; j < q_d0->size; j++){
        if (j < q_d0->size - 1){ // non-intercepts
            x_ij = matrix_get(cm_model->_x, i, j);
            vector_view x_j = matrix_column(cm_model->_x, j); // j'th covariate for all observations
            // ---- (deriv h_d / deriv theta_j):   ((g^{-1})'(propscore_i)/(n*gamma_derivative_n^2))sum_{j:D_j=d} X_{jk} * K'((singleindex_score_j-singleindex_score_i)/gamma_derivative_n)
            h_d_deriv_thetaj = nonpara_qk_deriv_eval_filter(vector_get(cm_model->_singleindex_score, i),
                                                            cm_model->_singleindex_score->data,
                                                            cm_model->_singleindex_score->size,
                                                            cm_model->_singleindex_score->stride,
                                                            (&x_j.vector)->data,
                                                            (&x_j.vector)->size,
                                                            (&x_j.vector)->stride,
                                                            1,  // exponent
                                                            nonpara_kernelderiv_gauss,
                                                            gamma_derivative_n,
                                                            -g_inv_deriv(vector_get(cm_model->propscore, i)), // scale
                                                            0,
                                                            (void *)cm_model->d,
                                                            cm_kernel_filter_func);
            // ---- (deriv qmu_d / deriv theta_j):   ((g^{-1})'(propscore_i)/(n*gamma_derivative_n^2))sum_{j:D_j=d} X_{jk} * Y_j * K'((singleindex_score_j-singleindex_score_i)/gamma_derivative_n)
            qmu_d_deriv_thetaj = nonpara_qk_deriv_weighted_eval_filter(vector_get(cm_model->_singleindex_score, i),
                                                                       cm_model->_singleindex_score->data,
                                                                       cm_model->_singleindex_score->size,
                                                                       cm_model->_singleindex_score->stride,
                                                                       cm_model->y->data,
                                                                       cm_model->y->size,
                                                                       cm_model->y->stride,
                                                                       (&x_j.vector)->data, // x_{jk} as weight
                                                                       (&x_j.vector)->size,
                                                                       (&x_j.vector)->stride,
                                                                       1,   // exponent
                                                                       nonpara_kernelderiv_gauss,
                                                                       gamma_derivative_n,
                                                                       -g_inv_deriv(vector_get(cm_model->propscore, i)), // scale
                                                                       0,
                                                                       (void *)cm_model->d,
                                                                       cm_kernel_filter_func);
        } else { // intercept
            x_ij = 1;
            // ---- (deriv h_d / deriv theta_{k+1}):   ((g^{-1})'(propscore_i)/(n*gamma_derivative_n^2))sum_{j:D_j=d} K'((singleindex_score_j-singleindex_score_i)/gamma_derivative_n)
            h_d_deriv_thetaj = nonpara_h_deriv_eval_filter(vector_get(cm_model->_singleindex_score, i),
                                                           cm_model->_singleindex_score->data,
                                                           cm_model->_singleindex_score->size,
                                                           cm_model->_singleindex_score->stride,
                                                           nonpara_kernelderiv_gauss,
                                                           gamma_derivative_n,
                                                           -g_inv_deriv(vector_get(cm_model->propscore, i)),
                                                           0,
                                                           (void *)cm_model->d,
                                                           cm_kernel_filter_func);
            // ---- (deriv qmu_d / deriv theta_{k+1}):   ((g^{-1})'(propscore_i)/(n*gamma_derivative_n^2))sum_{j:D_j=d} Y_j * K'((singleindex_score_j-singleindex_score_i)/gamma_derivative_n)
            qmu_d_deriv_thetaj = nonpara_qk_deriv_eval_filter(vector_get(cm_model->_singleindex_score, i),
                                                              cm_model->_singleindex_score->data,
                                                              cm_model->_singleindex_score->size,
                                                              cm_model->_singleindex_score->stride,
                                                              cm_model->y->data,
                                                              cm_model->y->size,
                                                              cm_model->y->stride,
                                                              1,    // exponent
                                                              nonpara_kernelderiv_gauss,
                                                              gamma_derivative_n,
                                                              -g_inv_deriv(vector_get(cm_model->propscore, i)), // scale
                                                              0,
                                                              (void *)cm_model->d,
                                                              cm_kernel_filter_func);
        }
        // derivatives w.r.t. j-th element in theta normalised by n
        mu_d_deriv_thetaj[0] = (qmu_d_deriv_thetaj[0] * h_d[0] - qmu_d[0] * h_d_deriv_thetaj[0]) / (n * pow(h_d[0], 2));
        mu_d_deriv_thetaj[1] = (qmu_d_deriv_thetaj[1] * h_d[1] - qmu_d[1] * h_d_deriv_thetaj[1]) / (n * pow(h_d[1], 2));
        // update q_d0
        update_q_d0 = mu_d_deriv_thetaj[0] + deriv_p[0] * x_ij; //  Lambda0_ij / n
        vector_set(q_d0, j, vector_get(q_d0, j) + update_q_d0); // update of coordinate j by element i estimates
        // update q_d1
        update_q_d1 = mu_d_deriv_thetaj[1] + deriv_p[1] * x_ij; // Lambda1_ij / n
        vector_set(q_d1, j, vector_get(q_d1, j) + update_q_d1); // update of coordinate j by element i estimates
        if (vector_short_get(cm_model->d, i)){ // update treated if and only if d_i=1
            // update q_d0_treated
            update_q_d0_treated = update_q_d0 / p1_hat; // Lambda0_ij / (n * p1_hat)
            vector_set(q_d0_treated, j, vector_get(q_d0_treated, j) + update_q_d0_treated);
            // update q_d1_treated
            update_q_d1_treated = update_q_d1 / p1_hat; // Lambda1_ij / (n * p1_hat)
            vector_set(q_d1_treated, j, vector_get(q_d1_treated, j) + update_q_d1_treated);
        }
    }
    // cleanup
    free(qmu_d_deriv_p);
    free(h_d_deriv_p);
    free(mu_d_deriv_p);
    free(deriv_p);
    free(qmu_d_deriv_thetaj);
    free(h_d_deriv_thetaj);
    return 0;
}


void test_cm_update_qs(void){
    // setup
    size_t n, i, k, j;
    n = 4;
    k = 2;
    double score, gamma_n, gamma_derivative_n;
    gamma_n = pow(n, -1.0/5);
    gamma_derivative_n = 0.1 * pow(n, -1.0/5);
    double y_i, score_i, propscore_i, x_ij;
    short d_i;
    vector *y = vector_alloc(n);
    vector_short *d = vector_short_alloc(n);
    matrix *x = matrix_alloc(n, k);
    vector *theta = vector_alloc(k+1);
    vector *singleindex_score = vector_alloc(n);
    vector *propscore = vector_alloc(n);
    double (*g)(double);
    double (*g_deriv)(double);
    double (*g_invderiv)(double);
    int modeltype_index = 0;  //
    g = propscore_g[modeltype_index];
    g_deriv = propscore_gderiv[modeltype_index];
    g_invderiv = propscore_ginvderiv[modeltype_index];
    double (*kernel_deriv)(double);
    kernel_deriv = nonpara_kernelderiv_gauss;
    vector_short_set(d, 0, 0);      vector_set(y, 0, -40.42);
    vector_short_set(d, 1, 0);      vector_set(y, 1, 1000.3);
    vector_short_set(d, 2, 1);      vector_set(y, 2, 30000.4);
    vector_short_set(d, 3, 1);      vector_set(y, 3, 20.96);
    vector_set(theta, 0, 0.5);
    vector_set(theta, 1, 0.22);
    vector_set(theta, 2, 0.06);
    matrix_set(x, 0, 0, 5.0);   matrix_set(x, 0, 1, 2.0);
    matrix_set(x, 1, 0, -3.3);  matrix_set(x, 1, 1, -1.43);
    matrix_set(x, 2, 0, -10.1); matrix_set(x, 2, 1, 1.96);
    matrix_set(x, 3, 0, 1.3);   matrix_set(x, 3, 1, 2.8);
    double p1_hat = vector_short_mean(d);
    // score and propscore for observations
    for (int j=0; j<n; j++){
        // score for observation j
        score = vector_get(theta, k) + vector_get(theta, 0) * matrix_get(x, j, 0) + vector_get(theta, 1) * matrix_get(x, j, 1);
        vector_set(singleindex_score, j, score);
        vector_set(propscore, j, g(score));
    }
    double *h_d = malloc(2*sizeof(double));
    double *qmu_d = malloc(2*sizeof(double));
    double *h_d_deriv_p = malloc(2*sizeof(double));
    double *qmu_d_deriv_p = malloc(2*sizeof(double));
    double *mu_d_deriv_p = malloc(2*sizeof(double));
    double *deriv_p = malloc(2*sizeof(double));
    double *h_d_deriv_thetaj = malloc(2*sizeof(double));
    double *qmu_d_deriv_thetaj = malloc(2*sizeof(double));
    double *mu_d_deriv_thetaj = malloc(2*sizeof(double));
    double *deriv_thetaj = malloc(2*sizeof(double));
    double update_q0, update_q0_treated, update_q1, update_q1_treated;
    vector *q_d0 = vector_calloc(k+1);
    vector *q_d0_treated = vector_calloc(k+1);
    vector *q_d1 = vector_calloc(k+1);
    vector *q_d1_treated = vector_calloc(k+1);
    vector *result_q_d0 = vector_calloc(k+1);
    vector *result_q_d0_treated = vector_calloc(k+1);
    vector *result_q_d1 = vector_calloc(k+1);
    vector *result_q_d1_treated = vector_calloc(k+1);
    CMModelKnownPropscore *cm_model = malloc(sizeof(CMModelKnownPropscore));
    cm_model->d = d;
    cm_model->y = y;
    cm_model->_x = x;
    cm_model->_modeltype_index = modeltype_index;
    cm_model->_singleindex_score = singleindex_score;
    cm_model->propscore = propscore;
    cm_model->_theta = theta;

    for (int i=0; i<n; i++){ // test i=0,...,3
        // test i
        h_d[0] = 0.001; h_d[1] = 0.1;  // made-up values: correct values are irrelevant for the test
        qmu_d[0] = 10.34; qmu_d[1] = 3.2;  // made-up values: correct values are irrelevant for the test
        d_i = vector_short_get(d, i);
        score_i = vector_get(singleindex_score, i);
        propscore_i = vector_get(propscore, i);
        // true values
        h_d_deriv_p[0] = (-1.0 / (n * pow(gamma_n, 2))) * (kernel_deriv((vector_get(propscore, 0) - propscore_i) / gamma_n) + // controls have index 0,1
                                                        kernel_deriv((vector_get(propscore, 1) - propscore_i) / gamma_n));
        h_d_deriv_p[1] = (-1.0 / (n * pow(gamma_n, 2))) * (kernel_deriv((vector_get(propscore, 2) - propscore_i) / gamma_n) + // treated have index 2,3
                                                        kernel_deriv((vector_get(propscore, 3) - propscore_i) / gamma_n));
        qmu_d_deriv_p[0] = (-1.0 / (n * pow(gamma_n, 2))) * (vector_get(y, 0) * kernel_deriv((vector_get(propscore, 0) - propscore_i) / gamma_n) + // controls have index 0,1
                                                            vector_get(y, 1) * kernel_deriv((vector_get(propscore, 1) - propscore_i) / gamma_n));
        qmu_d_deriv_p[1] = (-1.0 / (n * pow(gamma_n, 2))) * (vector_get(y, 2) * kernel_deriv((vector_get(propscore, 2) - propscore_i) / gamma_n) + // treated have index 2,3
                                                            vector_get(y, 3) * kernel_deriv((vector_get(propscore, 3) - propscore_i) / gamma_n));
        mu_d_deriv_p[0] = (qmu_d_deriv_p[0] * h_d[0] - qmu_d[0] * h_d_deriv_p[0]) / (n * pow(h_d[0], 2));
        mu_d_deriv_p[1] = (qmu_d_deriv_p[1] * h_d[1] - qmu_d[1] * h_d_deriv_p[1]) / (n * pow(h_d[1], 2));
        // printf("TEST: h_d_deriv_p[0] = %f, h_d_deriv_p[1] = %f \n", h_d_deriv_p[0], h_d_deriv_p[1]);
        // printf("TEST: qmu_d_deriv_p[0] = %f; qmu_d_deriv_p[1] = %f.\n", qmu_d_deriv_p[0], qmu_d_deriv_p[1]);
        // printf("TEST: mu_d_deriv_p[0] = %f; mu_d_deriv_p[1] = %f.\n", mu_d_deriv_p[0], mu_d_deriv_p[1]);
        //  entry 0
        j = 0;
        x_ij = matrix_get(x, i, j);
        deriv_p[0] = mu_d_deriv_p[0] * g_deriv(score_i) * x_ij;
        deriv_p[1] = mu_d_deriv_p[1] * g_deriv(score_i) * x_ij;
        h_d_deriv_thetaj[0] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (matrix_get(x, 0, j) * kernel_deriv((vector_get(singleindex_score, 0) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                            matrix_get(x, 1, j) * kernel_deriv((vector_get(singleindex_score, 1) - score_i) / gamma_derivative_n));
        h_d_deriv_thetaj[1] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (matrix_get(x, 2, j) * kernel_deriv((vector_get(singleindex_score, 2) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                            matrix_get(x, 3, j) * kernel_deriv((vector_get(singleindex_score, 3) - score_i) / gamma_derivative_n));
        qmu_d_deriv_thetaj[0] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (vector_get(y, 0) * matrix_get(x, 0, j) * kernel_deriv((vector_get(singleindex_score, 0) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                vector_get(y, 1) * matrix_get(x, 1, j) * kernel_deriv((vector_get(singleindex_score, 1) - score_i) / gamma_derivative_n));
        qmu_d_deriv_thetaj[1] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (vector_get(y, 2) * matrix_get(x, 2, j) * kernel_deriv((vector_get(singleindex_score, 2) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                vector_get(y, 3) * matrix_get(x, 3, j) * kernel_deriv((vector_get(singleindex_score, 3) - score_i) / gamma_derivative_n));    
        mu_d_deriv_thetaj[0] = (qmu_d_deriv_thetaj[0] * h_d[0] - qmu_d[0] * h_d_deriv_thetaj[0]) / (n * pow(h_d[0], 2));
        mu_d_deriv_thetaj[1] = (qmu_d_deriv_thetaj[1] * h_d[1] - qmu_d[1] * h_d_deriv_thetaj[1]) / (n * pow(h_d[1], 2));
        // printf("TEST: h_d_deriv_thetaj[0] = %f; h_d_deriv_thetaj[1] = %f.\n", h_d_deriv_thetaj[0], h_d_deriv_thetaj[1]);
        // printf("TEST: qmu_d_deriv_thetaj[0] = %f; qmu_d_deriv_thetaj[1] = %f.\n", qmu_d_deriv_thetaj[0], qmu_d_deriv_thetaj[1]);
        // printf("TEST: mu_d_deriv_thetaj[0] = %f, mu_d_deriv_thetaj[1] = %f \n", mu_d_deriv_thetaj[0], mu_d_deriv_thetaj[1]);
        deriv_thetaj[0] = mu_d_deriv_thetaj[0];
        deriv_thetaj[1] = mu_d_deriv_thetaj[1];
        update_q0 = deriv_p[0] + deriv_thetaj[0];
        update_q1 = deriv_p[1] + deriv_thetaj[1];
        update_q0_treated = d_i * (deriv_p[0] + deriv_thetaj[0]) / p1_hat;
        update_q1_treated = d_i * (deriv_p[1] + deriv_thetaj[1]) / p1_hat;
        vector_set(q_d0, j, vector_get(q_d0, j) + update_q0);
        vector_set(q_d1, j, vector_get(q_d1, j) + update_q1);
        vector_set(q_d0_treated, j, vector_get(q_d0_treated, j) + update_q0_treated);
        vector_set(q_d1_treated, j, vector_get(q_d1_treated, j) + update_q1_treated);
        //  entry 1
        j = 1;
        x_ij = matrix_get(x, i, j);
        deriv_p[0] = mu_d_deriv_p[0] * g_deriv(score_i) * x_ij;
        deriv_p[1] = mu_d_deriv_p[1] * g_deriv(score_i) * x_ij;
        h_d_deriv_thetaj[0] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (matrix_get(x, 0, j) * kernel_deriv((vector_get(singleindex_score, 0) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                            matrix_get(x, 1, j) * kernel_deriv((vector_get(singleindex_score, 1) - score_i) / gamma_derivative_n));
        h_d_deriv_thetaj[1] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (matrix_get(x, 2, j) * kernel_deriv((vector_get(singleindex_score, 2) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                            matrix_get(x, 3, j) * kernel_deriv((vector_get(singleindex_score, 3) - score_i) / gamma_derivative_n));
        qmu_d_deriv_thetaj[0] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (vector_get(y, 0) * matrix_get(x, 0, j) * kernel_deriv((vector_get(singleindex_score, 0) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                vector_get(y, 1) * matrix_get(x, 1, j) * kernel_deriv((vector_get(singleindex_score, 1) - score_i) / gamma_derivative_n));
        qmu_d_deriv_thetaj[1] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (vector_get(y, 2) * matrix_get(x, 2, j) * kernel_deriv((vector_get(singleindex_score, 2) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                vector_get(y, 3) * matrix_get(x, 3, j) * kernel_deriv((vector_get(singleindex_score, 3) - score_i) / gamma_derivative_n));    
        mu_d_deriv_thetaj[0] = (qmu_d_deriv_thetaj[0] * h_d[0] - qmu_d[0] * h_d_deriv_thetaj[0]) / (n * pow(h_d[0], 2));
        mu_d_deriv_thetaj[1] = (qmu_d_deriv_thetaj[1] * h_d[1] - qmu_d[1] * h_d_deriv_thetaj[1]) / (n * pow(h_d[1], 2));
        // printf("TEST: h_d_deriv_thetaj[0] = %f; h_d_deriv_thetaj[1] = %f.\n", h_d_deriv_thetaj[0], h_d_deriv_thetaj[1]);
        // printf("TEST: qmu_d_deriv_thetaj[0] = %f; qmu_d_deriv_thetaj[1] = %f.\n", qmu_d_deriv_thetaj[0], qmu_d_deriv_thetaj[1]);
        // printf("TEST: mu_d_deriv_thetaj[0] = %f, mu_d_deriv_thetaj[1] = %f \n", mu_d_deriv_thetaj[0], mu_d_deriv_thetaj[1]);
        deriv_thetaj[0] = mu_d_deriv_thetaj[0];
        deriv_thetaj[1] = mu_d_deriv_thetaj[1];
        update_q0 = deriv_p[0] + deriv_thetaj[0];
        update_q1 = deriv_p[1] + deriv_thetaj[1];
        update_q0_treated = d_i * (deriv_p[0] + deriv_thetaj[0]) / p1_hat;
        update_q1_treated = d_i * (deriv_p[1] + deriv_thetaj[1]) / p1_hat;
        vector_set(q_d0, j, vector_get(q_d0, j) + update_q0);
        vector_set(q_d1, j, vector_get(q_d1, j) + update_q1);
        vector_set(q_d0_treated, j, vector_get(q_d0_treated, j) + update_q0_treated);
        vector_set(q_d1_treated, j, vector_get(q_d1_treated, j) + update_q1_treated);
        //  entry 2: intercept
        j = 2;
        x_ij = 1;
        deriv_p[0] = mu_d_deriv_p[0] * g_deriv(score_i) * x_ij;
        deriv_p[1] = mu_d_deriv_p[1] * g_deriv(score_i) * x_ij;
        h_d_deriv_thetaj[0] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (1.0 * kernel_deriv((vector_get(singleindex_score, 0) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                1.0 * kernel_deriv((vector_get(singleindex_score, 1) - score_i) / gamma_derivative_n));
        h_d_deriv_thetaj[1] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (1.0 * kernel_deriv((vector_get(singleindex_score, 2) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                1.0 * kernel_deriv((vector_get(singleindex_score, 3) - score_i) / gamma_derivative_n));
        qmu_d_deriv_thetaj[0] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (vector_get(y, 0) * 1.0 * kernel_deriv((vector_get(singleindex_score, 0) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                    vector_get(y, 1) * 1.0 * kernel_deriv((vector_get(singleindex_score, 1) - score_i) / gamma_derivative_n));
        qmu_d_deriv_thetaj[1] = (g_invderiv(propscore_i) / (n * pow(gamma_derivative_n, 2))) * (vector_get(y, 2) * 1.0 * kernel_deriv((vector_get(singleindex_score, 2) - score_i) / gamma_derivative_n) + // controls have index 0,1
                                                                                    vector_get(y, 3) * 1.0 * kernel_deriv((vector_get(singleindex_score, 3) - score_i) / gamma_derivative_n));    
        mu_d_deriv_thetaj[0] = (qmu_d_deriv_thetaj[0] * h_d[0] - qmu_d[0] * h_d_deriv_thetaj[0]) / (n * pow(h_d[0], 2));
        mu_d_deriv_thetaj[1] = (qmu_d_deriv_thetaj[1] * h_d[1] - qmu_d[1] * h_d_deriv_thetaj[1]) / (n * pow(h_d[1], 2));
        // printf("TEST: h_d_deriv_thetaj[0] = %f; h_d_deriv_thetaj[1] = %f.\n", h_d_deriv_thetaj[0], h_d_deriv_thetaj[1]);
        // printf("TEST: qmu_d_deriv_thetaj[0] = %f; qmu_d_deriv_thetaj[1] = %f.\n", qmu_d_deriv_thetaj[0], qmu_d_deriv_thetaj[1]);
        // printf("TEST: mu_d_deriv_thetaj[0] = %f, mu_d_deriv_thetaj[1] = %f \n", mu_d_deriv_thetaj[0], mu_d_deriv_thetaj[1]);
        deriv_thetaj[0] = mu_d_deriv_thetaj[0];
        deriv_thetaj[1] = mu_d_deriv_thetaj[1];
        update_q0 = deriv_p[0] + deriv_thetaj[0];
        update_q1 = deriv_p[1] + deriv_thetaj[1];
        update_q0_treated = d_i * (deriv_p[0] + deriv_thetaj[0]) / p1_hat;
        update_q1_treated = d_i * (deriv_p[1] + deriv_thetaj[1]) / p1_hat;
        vector_set(q_d0, j, vector_get(q_d0, j) + update_q0);
        vector_set(q_d1, j, vector_get(q_d1, j) + update_q1);
        vector_set(q_d0_treated, j, vector_get(q_d0_treated, j) + update_q0_treated);
        vector_set(q_d1_treated, j, vector_get(q_d1_treated, j) + update_q1_treated);
        // values from function
        cm_update_qs(i, result_q_d0, result_q_d0_treated, result_q_d1, result_q_d1_treated, cm_model, p1_hat, gamma_n, gamma_derivative_n, h_d, qmu_d);
        // for (int l=0; l<k+1; l++){
        //     printf("q_d0[%d] = %f, result_q_d0[%d] = %f \n", l, vector_get(q_d0, l), l, vector_get(result_q_d0, l));
        //     printf("q_d0_treated[%d] = %f, result_q_d0_treated[%d] = %f \n", l, vector_get(q_d0_treated, l), l, vector_get(result_q_d0_treated, l));
        //     printf("q_d1[%d] = %f, result_q_d1[%d] = %f \n", l, vector_get(q_d1, l), l, vector_get(result_q_d1, l));
        //     printf("q_d1_treated[%d] = %f, result_q_d1_treated[%d] = %f \n", l, vector_get(q_d1_treated, l), l, vector_get(result_q_d1_treated, l));
        // }
        assert(vector_isequal(q_d0, result_q_d0));
        assert(vector_isequal(q_d1, result_q_d1));
        assert(vector_isequal(q_d0_treated, result_q_d0_treated));
        assert(vector_isequal(q_d1_treated, result_q_d1_treated));
    }
    // cleanup
    vector_free(y); vector_free(propscore); vector_free(singleindex_score);
    vector_short_free(d); matrix_free(x);
    vector_free(q_d0), vector_free(q_d0_treated); vector_free(q_d1), vector_free(q_d1_treated);
    free(h_d_deriv_p); free(qmu_d_deriv_p); free(mu_d_deriv_p); free(deriv_p);
    free(h_d_deriv_thetaj); free(qmu_d_deriv_thetaj); free(mu_d_deriv_thetaj); free(deriv_thetaj);
    vector_free(result_q_d0), vector_free(result_q_d0_treated); vector_free(result_q_d1), vector_free(result_q_d1_treated);
    free(cm_model); 
    printf("Tests for `cm_update_qs` completed: no issues found.\n");
}


typedef struct CMThreadArgumentsVAR {
    long threadid;
    CMModelKnownPropscore *cm_model;
    CMResults *results;
    size_t index_low;
    size_t index_high;
    double p1_hat;
    double *var_tau;
    double *var_tau_treated;
    double *var_sigmapi;
    double *var_sigmapi_treated;
    vector *q_d0;
    vector *q_d1;
    vector *q_d0_treated;
    vector *q_d1_treated;
    matrix *fisher_info;
} CMThreadArgumentsVAR;

void *cm_variance_estimator_thread(void *thread_arguments){
    CMThreadArgumentsVAR *thread_args = (CMThreadArgumentsVAR *)thread_arguments;
    size_t n = thread_args->cm_model->d->size;
    int k = thread_args->cm_model->_called_internally ? thread_args->cm_model->_x->size2 : 0; // nr. of covarites without intercept
    double p1_hat = thread_args->p1_hat; // sample proportion of treated
    // variance component contributions
    double var_tau, var_tau_treated, var_sigmapi, var_sigmapi_treated, var_estpi, var_estpi_treated;
    var_tau = var_tau_treated = var_sigmapi = var_sigmapi_treated = var_estpi = var_estpi_treated = 0;
    vector *q_d0;
    vector *q_d1;
    vector *q_d0_treated;
    vector *q_d1_treated;
    matrix *fisher_info;
    if (thread_args->cm_model->_called_internally){ // third variance component is zero for known propscore: only allocate if necessary
        q_d0 = vector_calloc(k + 1);
        q_d1 = vector_calloc(k + 1);
        q_d0_treated = vector_calloc(k + 1);
        q_d1_treated = vector_calloc(k + 1);
        fisher_info = matrix_calloc(k + 1, k + 1);
    }
    // truncation
    double truncation_low, truncation_high;
    double pscore_min, pscore_max;
    double npower_a;
    double kappa_a;
    if (thread_args->cm_model->alpha == 0.0) { // if zero exponent is supplied, use default value
        npower_a = 1/4.000000002; // default value //1.0 / 5.5;
    } else {
        npower_a = thread_args->cm_model->alpha;
    }
    if (thread_args->cm_model->kappa_a == 0.0){ // if zero scale is supplied, use default value
        kappa_a = pow(10.0, -15); // default value
    } else {
        kappa_a = thread_args->cm_model->kappa_a;
    }
    double a_n = kappa_a * pow((double)n, -npower_a); // trunction sequence
    pscore_min = vector_get(thread_args->cm_model->propscore, thread_args->cm_model->_pscore_sortindex[0]);
    pscore_max = vector_get(thread_args->cm_model->propscore, thread_args->cm_model->_pscore_sortindex[n - 1]);
    if (CM_VERBOSE_MODE){
        printf("pscore_min = %f, pscore_max = %f.\n", pscore_min, pscore_max);
    }
    truncation_low = pscore_min + a_n;
    truncation_high = pscore_max - a_n;
    // bandwidth
    double npower_gamma;
    double kappa_gamma;
    if (thread_args->cm_model->beta == 0.0){ // if zero exponent is supplied, use default value
        npower_gamma = 1/4.000000001; // default value //1.0 / 5; 
    } else {
        npower_gamma = thread_args->cm_model->beta;
    }
    if (thread_args->cm_model->kappa_gamma == 0.0){  //if zero scale is supplied, use default value
        kappa_gamma = 0.5; // default value
    } else {
        kappa_gamma = thread_args->cm_model->kappa_gamma;
    }
    assert((npower_a < npower_gamma) && (0 < npower_a) && (0 < npower_gamma) && (npower_gamma < 1.0 / 4)); // sanity check
    double gamma_n = kappa_gamma * pow((double)n, -npower_gamma);  // bandwidth sequence
    // bandwidth for derivative estimation used in var_estpi 
    double kappa_gamma_derivative;
    double gamma_derivative_n = 0.0;
    if (thread_args->cm_model->_called_internally){  // var_estpi, hence derivatives, estimated too
        if (thread_args->cm_model->_kappa_gamma_derivative == 0.0){ // if zero scale is supplied, use default value
            kappa_gamma_derivative = 0.5; // default value
        } else {
            kappa_gamma_derivative = thread_args->cm_model->_kappa_gamma_derivative;
        }
        //scale may be different, but same power `npower_gamma` of n is used
        gamma_derivative_n = kappa_gamma_derivative * pow((double)n, -npower_gamma); // bandwidth sequence for derivatives
    }
    // process observations assigned to the thread
    for (size_t i = thread_args->index_low; i < thread_args->index_high; i++){
        // only compute nonparametric estimates for propensity score that are not truncated
        if (cm_inrange(vector_get(thread_args->cm_model->propscore, i), truncation_low, truncation_high)){
            // obtain nonparametric estimates
            double *h_d = nonpara_h_eval_filter(vector_get(thread_args->cm_model->propscore, i),
                                                thread_args->cm_model->propscore->data,
                                                thread_args->cm_model->propscore->size,
                                                thread_args->cm_model->propscore->stride,
                                                nonpara_kernel_gauss,
                                                gamma_n,
                                                0,
                                                (void *)thread_args->cm_model->d,
                                                cm_kernel_filter_func);
            double *qmu_d = nonpara_qk_eval_filter(vector_get(thread_args->cm_model->propscore, i),
                                                   thread_args->cm_model->propscore->data,
                                                   thread_args->cm_model->propscore->size,
                                                   thread_args->cm_model->propscore->stride,
                                                   thread_args->cm_model->y->data,
                                                   thread_args->cm_model->y->size,
                                                   thread_args->cm_model->y->stride,
                                                   1.0,
                                                   nonpara_kernel_gauss,
                                                   gamma_n,
                                                   0,
                                                   (void *)thread_args->cm_model->d,
                                                   cm_kernel_filter_func);
            double *qmu2_d = nonpara_qk_eval_filter(vector_get(thread_args->cm_model->propscore, i),
                                                    thread_args->cm_model->propscore->data,
                                                    thread_args->cm_model->propscore->size,
                                                    thread_args->cm_model->propscore->stride,
                                                    thread_args->cm_model->y->data,
                                                    thread_args->cm_model->y->size,
                                                    thread_args->cm_model->y->stride,
                                                    2.0,
                                                    nonpara_kernel_gauss,
                                                    gamma_n,
                                                    0,
                                                    (void *)thread_args->cm_model->d,
                                                    cm_kernel_filter_func);
            double *mu_d = malloc(2 * sizeof(double));     // first moment
            mu_d[0] = qmu_d[0] / h_d[0];                   // mu_0 control
            mu_d[1] = qmu_d[1] / h_d[1];                   // mu_1 treated
            double *mu2_d = malloc(2 * sizeof(double));    // second moment
            mu2_d[0] = qmu2_d[0] / h_d[0];                 // mu2_0 control
            mu2_d[1] = qmu2_d[1] / h_d[1];                 // mu2_1 treated
            double *sigma2_d = malloc(2 * sizeof(double)); // variance
            sigma2_d[0] = mu2_d[0] - pow(mu_d[0], 2);      // sigma2_0 control
            sigma2_d[1] = mu2_d[1] - pow(mu_d[1], 2);      // sigma2_1 treated
            //
            // printf("h_d[0] = %f; h_d[1] = %f.\n", h_d[0], h_d[1]);
            // printf("qmu_d[0] = %f; qmu_d[1] = %f.\n", qmu_d[0], qmu_d[1]);
            // printf("qmu2_d[0] = %f; qmu2_d[1] = %f.\n", qmu2_d[0], qmu2_d[1]);
            // printf("mu_d[0] = %f; mu_d[1] = %f.\n", mu_d[0], mu_d[1]);
            // printf("mu2_d[0] = %f; mu2_d[1] = %f.\n", mu2_d[0], mu2_d[1]);
            // printf("sigma2_d[0] = %f; sigma2_d[1] = %f.\n", sigma2_d[0], sigma2_d[1]);
            // update variance contributions
            // var tau
            var_tau += pow(mu_d[1] - mu_d[0], 2) / n;
            // printf("BEFORE: var_tau = %f, var_tau_treated = %f\n", var_tau, var_tau_treated);
            if (vector_short_get(thread_args->cm_model->d, i)){ // only deal with treated obs
                var_tau_treated += pow(mu_d[1] - mu_d[0], 2) / (n * pow(p1_hat, 2));
            }
            // printf("AFTER: var_tau = %f, var_tau_treated = %f\n", var_tau, var_tau_treated);
            // var sigma pi
            var_sigmapi += sigma2_d[0] / (n * (1 - vector_get(thread_args->cm_model->propscore, i))) + sigma2_d[1] / (n * vector_get(thread_args->cm_model->propscore, i));
            var_sigmapi_treated += pow(vector_get(thread_args->cm_model->propscore, i), 2) * sigma2_d[0] / (n * pow(p1_hat, 2) * (1 - vector_get(thread_args->cm_model->propscore, i))) + vector_get(thread_args->cm_model->propscore, i) * sigma2_d[1] / (n * pow(p1_hat, 2));
            // printf("var_sigmapi = %f, var_sigmapi_treated = %f\n", var_sigmapi, var_sigmapi_treated);
            if (thread_args->cm_model->_called_internally){ // compute variance contribution for estimated propscore
                cm_update_qs(i, q_d0, q_d0_treated, q_d1, q_d1_treated, thread_args->cm_model, p1_hat, gamma_n, gamma_derivative_n, h_d, qmu_d);
            }
            // cleanup
            free(qmu_d);
            free(qmu2_d);
            free(h_d);
            free(mu_d);
            free(mu2_d);
            free(sigma2_d);
        }
        if (thread_args->cm_model->_called_internally){ // compute variance contribution for estimated propscore
            cm_update_fisher_info(i, fisher_info, thread_args->cm_model);
        }
    }
    // merge in results
    pthread_mutex_lock(&mutex_lock);
    // truncation and bandwidth sequences
    thread_args->results->kappa_a = kappa_a;
    thread_args->results->kappa_gamma = kappa_gamma;
    thread_args->results->kappa_gamma_derivative = kappa_gamma_derivative;
    thread_args->results->alpha = npower_a;
    thread_args->results->beta = npower_gamma;
    thread_args->results->a_n = a_n;
    thread_args->results->gamma_n = gamma_n;
    thread_args->results->gamma_derivative_n = gamma_derivative_n;  // zero if propscore is knonw
    thread_args->results->propscore_min = pscore_min;
    thread_args->results->propscore_max = pscore_max;
    thread_args->results->truncation_low = truncation_low;
    thread_args->results->truncation_high = truncation_high;
    // variance components
    // printf("*(thread_args->var_tau) before = %f\n", *(thread_args->var_tau));
    *(thread_args->var_tau) += var_tau;
    // printf("*(thread_args->var_tau) after = %f\n", *(thread_args->var_tau));
    // printf("*(thread_args->var_tau_treated) before = %f\n", *(thread_args->var_tau_treated));
    *(thread_args->var_tau_treated) += var_tau_treated;
    // printf("*(thread_args->var_tau_treated) before = %f\n", *(thread_args->var_tau_treated));
    *(thread_args->var_sigmapi) += var_sigmapi;
    *(thread_args->var_sigmapi_treated) += var_sigmapi_treated;
    // printf("BEFORE:\n");
    // for (int j = 0; j < thread_args->cm_model->_theta->size; j++){
    //     printf("thread_args->q_d0[%d] = %f\n", j, vector_get(thread_args->q_d0, j));
    //     printf("thread_args->q_d1[%d] = %f\n", j, vector_get(thread_args->q_d1, j));
    //     printf("thread_args->q_d0_treated[%d] = %f\n", j, vector_get(thread_args->q_d0_treated, j));
    //     printf("thread_args->q_d1_treated[%d] = %f\n", j, vector_get(thread_args->q_d1_treated, j));
    // }
    if (thread_args->cm_model->_called_internally){
        matrix_add(thread_args->fisher_info, fisher_info);
        vector_add(thread_args->q_d0, q_d0);
        vector_add(thread_args->q_d1, q_d1);
        vector_add(thread_args->q_d0_treated, q_d0_treated);
        vector_add(thread_args->q_d1_treated, q_d1_treated);
        if (CM_VERBOSE_MODE){    
            printf("AFTER:\n");
            for (int j = 0; j < thread_args->cm_model->_theta->size; j++){
                printf("thread_args->q_d0[%d] = %f\n", j, vector_get(thread_args->q_d0, j));
                printf("thread_args->q_d1[%d] = %f\n", j, vector_get(thread_args->q_d1, j));
                printf("thread_args->q_d0_treated[%d] = %f\n", j, vector_get(thread_args->q_d0_treated, j));
                printf("thread_args->q_d1_treated[%d] = %f\n", j, vector_get(thread_args->q_d1_treated, j));
            }
        }
    }
    pthread_mutex_unlock(&mutex_lock);
    void *ptr;
    return ptr;
}

// Estimates asymptotic variance variance of matching estimator
void cm_variance_estimator(CMModelKnownPropscore *cm_model, CMResults *results){
    size_t n = cm_model->d->size;
    size_t k = cm_model->_called_internally ? cm_model->_x->size2 : 0;
    double var_tau, var_tau_treated, var_sigmapi, var_sigmapi_treated;
    var_tau = var_tau_treated = var_sigmapi = var_sigmapi_treated = 0;
    double p1_hat = vector_short_mean(cm_model->d);
    vector *q_d0;
    vector *q_d1;
    vector *q_d0_treated; // the factor (1/p1_hat) is assumed to be included in q_d0_treated
    vector *q_d1_treated; // the factor (1/p1_hat) is assumed to be included in q_d1_treated
    matrix *fisher_info;
    if (cm_model->_called_internally){ // third variance component is zero for known propscore (=if not _called_internally): only allocate if necessary
        q_d0 = vector_calloc(k + 1);
        q_d1 = vector_calloc(k + 1);
        q_d0_treated = vector_calloc(k + 1);
        q_d1_treated = vector_calloc(k + 1);
        fisher_info = matrix_calloc(k + 1, k + 1);
    }
    if (NUM_THREADS > 1){ // multithreading
        int result_code;
        int obs_per_thread = (int)n / NUM_THREADS;
        pthread_t threads[NUM_THREADS];
        CMThreadArgumentsVAR thread_args[NUM_THREADS];
        for (long t = 0; t < NUM_THREADS; t++){
            if (CM_VERBOSE_MODE){
                printf("`cm_variance_estimator`: creating thread %ld.\n", t);
            }
            thread_args[t].threadid = t;
            thread_args[t].cm_model = cm_model;
            thread_args[t].results = results;
            thread_args[t].index_low = t * obs_per_thread;
            thread_args[t].index_high = (t == NUM_THREADS - 1 ? n : (t + 1) * obs_per_thread);
            thread_args[t].p1_hat = p1_hat;
            thread_args[t].var_tau = &var_tau;
            thread_args[t].var_tau_treated = &var_tau_treated;
            thread_args[t].var_sigmapi = &var_sigmapi;
            thread_args[t].var_sigmapi_treated = &var_sigmapi_treated;
            thread_args[t].q_d0 = q_d0;
            thread_args[t].q_d1 = q_d1;
            thread_args[t].q_d0_treated = q_d0_treated;
            thread_args[t].q_d1_treated = q_d1_treated;
            thread_args[t].fisher_info = fisher_info;
            if (CM_VERBOSE_MODE){
                printf("`cm_variance_estimator`: thread %ld index_high is %zu.\n", t, thread_args[t].index_high);
            }
            result_code = pthread_create(&threads[t], NULL, cm_variance_estimator_thread, (void *)&thread_args[t]);
            assert(!result_code);
        }
        // wait for threads to finish
        for (long t = 0; t < NUM_THREADS; t++){
            result_code = pthread_join(threads[t], NULL);
            assert(!result_code);
            if (CM_VERBOSE_MODE){
                printf("`cm_variance_estimator`: thread %ld has finished.\n", t);
            }
        }
    } else { // single thread
        CMThreadArgumentsVAR *thread_args = malloc(sizeof(CMThreadArgumentsVAR));
        thread_args->cm_model = cm_model;
        thread_args->results = results;
        thread_args->index_low = 0;
        thread_args->index_high = n;
        thread_args->p1_hat = p1_hat;
        thread_args->var_tau = &var_tau;
        thread_args->var_tau_treated = &var_tau_treated;
        thread_args->var_sigmapi = &var_sigmapi;
        thread_args->var_sigmapi_treated = &var_sigmapi_treated;
        thread_args->q_d0 = q_d0;
        thread_args->q_d1 = q_d1;
        thread_args->q_d0_treated = q_d0_treated;
        thread_args->q_d1_treated = q_d1_treated;
        thread_args->fisher_info = fisher_info;
        cm_variance_estimator_thread((void *)thread_args);
    }
    // adjust var_tau and var_tau_treated for values not dealt with by threads
    if (CM_VERBOSE_MODE){
        printf("var_tau_treated after = %f\n", var_tau_treated);
    }
    var_tau -= pow(results->ate_hat, 2);
    var_tau_treated -= pow(results->att_hat, 2) / p1_hat;
    if (CM_VERBOSE_MODE){
        printf("var_tau_treated after after = %f\n", var_tau_treated);
        printf("results->att_hat = %f\n", results->att_hat);
        printf("pow(results->att_hat, 2) = %f\n", pow(results->att_hat, 2));
        printf("p1_hat = %f\n", p1_hat);
        printf("var_tau after after = %f\n", var_tau);
    }
    // compute third variance component coming from propscore estimation
    double var_estpi = 0;
    double var_estpi_treated = 0;
    if (cm_model->_called_internally){ // estimated propscore
        if (CM_VERBOSE_MODE){
            for (int j = 0; j < cm_model->_theta->size; j++){
                for (int l = 0; l < cm_model->_theta->size; l++){
                    printf("fisher_info[%d,%d] = %f\n", j, l, matrix_get(fisher_info, j, l));
                }
            }
        }
        double *v_estpi = cm_compute_var_estpi(q_d0, q_d0_treated, q_d1, q_d1_treated, fisher_info);  // in-place modifies q's and fisher_info!
        var_estpi = v_estpi[0];
        var_estpi_treated = v_estpi[1];
        // cleanup
        vector_free(q_d0); vector_free(q_d0_treated);
        vector_free(q_d1); vector_free(q_d1_treated);
        matrix_free(fisher_info);
    }
    // combine variance components
    double var_ate = var_tau + var_sigmapi + var_estpi;
    double var_att = var_tau_treated + var_sigmapi_treated + var_estpi_treated;
    // load to results
    results->var_hat_ate = var_ate;
    results->var_hat_att = var_att;
    results->var_hat_component_tau_ate = var_tau;
    results->var_hat_component_tau_att = var_tau_treated;
    results->var_hat_component_sigmapi_ate = var_sigmapi;
    results->var_hat_component_sigmapi_att = var_sigmapi_treated;
    results->var_hat_component_estpi_ate = var_estpi;
    results->var_hat_component_estpi_att = var_estpi_treated;
}

typedef struct CMThreadArgumentsMaxmin {
    long threadid;
    vector_short *d;
    vector *propscore;
    size_t *pscore_sortindex;
    size_t sortindex_low;
    size_t sortindex_high;
    double *maxmins;
} CMThreadArgumentsMaxmin;


void *cm_maxmin_delta_thread(void *thread_arguments){
    CMThreadArgumentsMaxmin *thread_args = (CMThreadArgumentsMaxmin *)thread_arguments;
    double maxmin = 0.0;
    double min_left, min_right, min;
    for (int i=thread_args->sortindex_low; i<thread_args->sortindex_high; i++){
        // largest possible distance between propensity scores in case there are no opposite treatment-group unit to the left or right
        min_left = min_right = 1.0;
        // find opposite treatment-group units
        for (int j=i; j>=0; j--){  // look to left
            if (vector_short_get(thread_args->d, thread_args->pscore_sortindex[j]) != vector_short_get(thread_args->d, thread_args->pscore_sortindex[i])){
                // distance between i-th and j-th smallest units' propensity scores
                min_left = cm_abs_distance(vector_get(thread_args->propscore, thread_args->pscore_sortindex[j]),
                                           vector_get(thread_args->propscore, thread_args->pscore_sortindex[i]));
                break;  // if found an opposite treatment group unit, then it is the closest one to the left, so stop looking
            }
        }
        for (int j=i; j<thread_args->d->size; j++){  // look to the right
            if (vector_short_get(thread_args->d, thread_args->pscore_sortindex[j]) != vector_short_get(thread_args->d, thread_args->pscore_sortindex[i])){
                // distance between i-th and j-th smallest units' propensity scores
                min_right = cm_abs_distance(vector_get(thread_args->propscore, thread_args->pscore_sortindex[j]),
                                            vector_get(thread_args->propscore, thread_args->pscore_sortindex[i]));
                break;  // if found an opposite treatment group unit, then it is the closest one to the right, so stop looking
            }
        }
        // take minimum of left- and right-minimums
        min = (min_left <= min_right) ? min_left : min_right;    
        // update maxmin distance if found unit's distance is larger
        if (maxmin < min){
            maxmin = min;
        }
    }
    // load local maxmin
    pthread_mutex_lock(&mutex_lock);
    thread_args->maxmins[thread_args->threadid] = maxmin;
    pthread_mutex_unlock(&mutex_lock);
    void *ptr;
    return ptr;
}

/* Computes maxmin distance max_{i\in[n]}min_{j:d_j!=d_i}|propscore_i-propscore_j|,
   where d_i is the i-th element of `d`, and similarly for propscore_i, propscore_j.
   `pscore_sortindex` is the permutation of the indices {0,...,n-1} that sorts
   `propscore` into ascending order. I.e. the first element of `pscore_sortindex` is
   the index of the smallest element in `propscore`; the last element is the index of
   the largest element in `propscore`, etc..
*/
double cm_maxmin_delta(vector_short *d, vector *propscore, size_t *pscore_sortindex){
    assert(vector_short_is_zero_one(d));
    int n_treated = vector_short_sum(d);
    assert((n_treated < d->size) && (0 < n_treated));   // ensure there's at least one unit from each treatment group
    // maxmins in each batch processed by a single thread
    double *maxmins = malloc(NUM_THREADS * sizeof(double));
    size_t n = d->size;
    if (NUM_THREADS > 1){   // multi-threaded version
        int result_code;
        int obs_per_thread = (int)n / NUM_THREADS;
        pthread_t threads[NUM_THREADS];
        CMThreadArgumentsMaxmin thread_args[NUM_THREADS];
        for (long t = 0; t < NUM_THREADS; t++){
            if (CM_VERBOSE_MODE){
                printf("`cm_maxmin_delta`: creating thread %ld.\n", t);
            }
            thread_args[t].threadid = t;
            thread_args[t].d = d;
            thread_args[t].propscore = propscore;
            thread_args[t].pscore_sortindex = pscore_sortindex;
            thread_args[t].maxmins = maxmins;
            thread_args[t].sortindex_low = t * obs_per_thread;
            thread_args[t].sortindex_high = (t == NUM_THREADS - 1 ? n : (t + 1) * obs_per_thread);
            if (CM_VERBOSE_MODE){
                printf("`cm_maxmin_delta`: thread %ld sortindex_high is %zu.\n", t, thread_args[t].sortindex_high);
            }
            result_code = pthread_create(&threads[t], NULL, cm_maxmin_delta_thread, (void *)&thread_args[t]);
            assert(!result_code);
        }
        // wait for threads to finish
        for (long t = 0; t < NUM_THREADS; t++){
            result_code = pthread_join(threads[t], NULL);
            assert(!result_code);
            if (CM_VERBOSE_MODE){
                printf("`cm_maxmin_delta`: thread %ld has finished.\n", t);
            }
        }
    } else {    // single-threaded version
        CMThreadArgumentsMaxmin *thread_args = malloc(sizeof(CMThreadArgumentsMaxmin));
        thread_args->threadid = 0;
        thread_args->d = d;
        thread_args->propscore = propscore;
        thread_args->pscore_sortindex = pscore_sortindex;
        thread_args->maxmins = maxmins;
        thread_args->sortindex_low = 0;
        thread_args->sortindex_high = n;
        cm_maxmin_delta_thread((void *)thread_args);
    }
    // merge: find largest maxmin across batches processed by threads
    double maxmin = 0.0;
    for (long t=0; t<NUM_THREADS; t++){
        if (maxmin <= maxmins[t]){
            maxmin = maxmins[t];
            // printf("maxmin = %f.\n", maxmin);
        }
    }
    return maxmin;
}

void test_cm_maxmin_delta(void){
    // setup
    double maxmin;
    size_t n = 4;
    vector_short *d = vector_short_alloc(n);
    vector *propscore = vector_alloc(n);
    size_t *pscore_sortindex = malloc(n * sizeof(size_t));
    
    // test 1
    vector_short_set(d, 0, 0);
    vector_set(propscore, 0, 0.1); // c1: first control unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.5); // c2: second control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.9); // t1: first treated unit
    vector_short_set(d, 3, 1);
    vector_set(propscore, 3, 0.2); // t2: second treated unit
    pscore_sortindex[0] = 0;  // index of smallest propscore
    pscore_sortindex[1] = 3;
    pscore_sortindex[2] = 1;
    pscore_sortindex[3] = 2;  // index of largest propscore
    // distance_matrix =
    //      c1    c2     t1      t2         min_{j:D_j!=D_i}|propscore_i-propscore_j|
    //      _      _     _       _          -
    // c1 | .      .     0.8,    0.1        0.1
    // c2 | .      .     0.4     0.3        0.3
    // t1 | 0.8,   0.4    .       .         0.4         
    // t2 | 0.1,   0.3    .       .         0.3
    // -----------------------------------------
    // max| .........................       0.4
    maxmin = 0.4;
    // printf("cm_maxmin_delta(d, propscore, pscore_sortindex) = %f\n", cm_maxmin_delta(d, propscore, pscore_sortindex));
    assert(fabs(cm_maxmin_delta(d, propscore, pscore_sortindex) - maxmin) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    
    // test 2
    vector_short_set(d, 0, 1);
    vector_set(propscore, 0, 0.9); // t1: first treated unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.1); // c1: first control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.2); // t2: second treated unit
    vector_short_set(d, 3, 0);
    vector_set(propscore, 3, 0.5); // c2: second control unit
    pscore_sortindex[0] = 1;  // index of smallest propscore
    pscore_sortindex[1] = 2;
    pscore_sortindex[2] = 3;
    pscore_sortindex[3] = 0;  // index of largest propscore
    // distance_matrix =
    //      c1    c2     t1      t2         min_{j:D_j!=D_i}|propscore_i-propscore_j|
    //      _      _     _       _          -
    // c1 | .      .     0.8,    0.1        0.1
    // c2 | .      .     0.4     0.3        0.3
    // t1 | 0.8,   0.4    .       .         0.4         
    // t2 | 0.1,   0.3    .       .         0.3
    // -----------------------------------------
    // max| .........................       0.4
    maxmin = 0.4;
    // printf("cm_maxmin_delta(d, propscore, pscore_sortindex) = %f\n", cm_maxmin_delta(d, propscore, pscore_sortindex));
    assert(fabs(cm_maxmin_delta(d, propscore, pscore_sortindex) - maxmin) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));

    // test 3
    vector_short_set(d, 0, 0);
    vector_set(propscore, 0, 0.9); // c1: first control unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.1); // c2: second control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.8); // t1: first treated unit
    vector_short_set(d, 3, 0);
    vector_set(propscore, 3, 0.5); // c3: thrid control unit
    pscore_sortindex[0] = 1;  // index of smallest propscore
    pscore_sortindex[1] = 3;
    pscore_sortindex[2] = 2;
    pscore_sortindex[3] = 0;  // index of largest propscore
    // distance_matrix =
    //      c1    c2     c3      t1         min_{j:D_j!=D_i}|propscore_i-propscore_j|
    //      _      _     _       _          -
    // c1 | .      .     .      0.1        0.1
    // c2 | .      .     .      0.7        0.7
    // c3 | .      .     .      0.3        0.3         
    // t1 | 0.1,  0.7,  0.3      .         0.1
    // -----------------------------------------
    // max| .........................      0.7
    maxmin = 0.7;
    // printf("cm_maxmin_delta(d, propscore, pscore_sortindex) = %f\n", cm_maxmin_delta(d, propscore, pscore_sortindex));
    assert(fabs(cm_maxmin_delta(d, propscore, pscore_sortindex) - maxmin) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));

    // test 4
    vector_short_set(d, 0, 1);
    vector_set(propscore, 0, 0.9); // t1: first treated unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.2); // c1: first control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.8); // t2: second treated unit
    vector_short_set(d, 3, 1);
    vector_set(propscore, 3, 0.1); // t3: third treated unit
    pscore_sortindex[0] = 3;  // index of smallest propscore
    pscore_sortindex[1] = 1;
    pscore_sortindex[2] = 2;
    pscore_sortindex[3] = 0;  // index of largest propscore
    // distance_matrix =
    //      t1    t2     t3      c1         min_{j:D_j!=D_i}|propscore_i-propscore_j|
    //      _      _     _       _          -
    // t1 | .      .     .      0.7        0.7
    // t2 | .      .     .      0.6        0.6
    // t3 | .      .     .      0.1        0.1         
    // c1 | 0.7,  0.6,  0.1      .         0.1
    // -----------------------------------------
    // max| .........................      0.7
    maxmin = 0.7;
    // printf("cm_maxmin_delta(d, propscore, pscore_sortindex) = %f\n", cm_maxmin_delta(d, propscore, pscore_sortindex));
    assert(fabs(cm_maxmin_delta(d, propscore, pscore_sortindex) - maxmin) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));

    // test 5
    vector_short_set(d, 0, 0);
    vector_set(propscore, 0, 0.6); // c1: first control unit
    vector_short_set(d, 1, 0);
    vector_set(propscore, 1, 0.2); // c2: second control unit
    vector_short_set(d, 2, 1);
    vector_set(propscore, 2, 0.7); // t1: first treated unit
    vector_short_set(d, 3, 1);
    vector_set(propscore, 3, 0.3); // t2: second treated unit
    pscore_sortindex[0] = 1;  // index of smallest propscore
    pscore_sortindex[1] = 3;
    pscore_sortindex[2] = 0;
    pscore_sortindex[3] = 2;  // index of largest propscore
    // distance_matrix =
    //      c1    c2     t1      t2         min_{j:D_j!=D_i}|propscore_i-propscore_j|
    //      _      _     _       _          -
    // c1 | .      .     0.1,    0.3        0.1
    // c2 | .      .     0.5     0.1        0.1
    // t1 | 0.1,   0.5    .       .         0.1         
    // t2 | 0.3,   0.1    .       .         0.1
    // -----------------------------------------
    // max| .........................       0.1
    maxmin = 0.1;
    // printf("cm_maxmin_delta(d, propscore, pscore_sortindex) = %f\n", cm_maxmin_delta(d, propscore, pscore_sortindex));
    assert(fabs(cm_maxmin_delta(d, propscore, pscore_sortindex) - maxmin) <= pow(10, -CM_TEST_DOUBLE_TOLERANCE_NPOWER));
    

    printf("Tests for `cm_maxmin_delta` completed: no issues found.\n");
}


//////////////////      API




// Sets caliper. Only to be called by the user, not internally.
int cm_initialise_known_propscore(CMModelKnownPropscore *cm_model){
    // validate inputs
    if (cm_cmmodelknownpropscore_init_check_values(cm_model, 0)){
        exit(EXIT_FAILURE);
    };
    // sample size
    size_t n = cm_model->d->size;
    // sort propensity score
    size_t *pscore_sortindex = vector_sort_index(cm_model->propscore);
    // set caliper
    size_t n1 = vector_short_sum(cm_model->d);
    size_t n0 = n - n1;
    double max_groups = cm_d_max(log(n1 + 0.0) / n1, log(n0 + 0.0) / n0);
    // double delta = cm_maxmin_delta(cm_model->d, propscore, pscore_sortindex);
    double delta = 0.0;
    if (cm_model->delta == 0.0){
        delta = cm_d_max(max_groups, cm_maxmin_delta(cm_model->d, cm_model->propscore, pscore_sortindex));
    } else {
        delta = cm_model->delta;
    }     
    // write to model
    cm_model->n = n;
    cm_model->_pscore_sortindex = pscore_sortindex;
    cm_model->delta = delta;
    // instantiation successfull
    cm_model->_called_internally = 0;
    cm_model->_isinstantiated = 1;
    return 0;
}

CMResults *cm_cm_known_propscore(CMModelKnownPropscore *cm_model){
    // validate inputs
    if (cm_cmmodelknownpropscore_check_values(cm_model, 0)){
        exit(EXIT_FAILURE);
    };
    // ensure model has been instantiated
    if (!cm_model->_isinstantiated){
        fprintf(stderr, "%s: line %d: CMModelKnownPropscore has not been instantiated prior to calling `cm_cm_known_propscore`. Do so with `cm_initialise_known_propscore`.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    CMResults *results = malloc(sizeof(CMResults));
    results->delta = cm_model->delta;  // write caliper to results too
    // estimates
    cm_matches(cm_model, results);
    // some sanity checks
    // printf("delta = %f", cm_model->delta);
    // for (int i=0; i<cm_model->d->size; i++){
    //     for (int j=0; j<cm_model->d->size; j++){
    //         if (dynamicarray_int_contains(results->match_indices[i], j)){
    //             printf("\nd[%d] = %hd, d[%d] = %hd", i, vector_short_get(cm_model->d, i), j, vector_short_get(cm_model->d, j));
    //             printf("; p[%d] = %.3f, p[%d] = %.3f", i, vector_get(cm_model->propscore, i), j, vector_get(cm_model->propscore, j));
    //         }
    //     }
    //     printf("\n number of matches for unit %d is %d.\n", i, vector_int_get(results->number_of_matches, i));
    // }
    cm_te_estimator(cm_model, results);
    cm_variance_estimator(cm_model, results);
    return results;
}


// Compute singleindex-score, propensity score, and sets caliper.
int cm_initialise(CMModel *cm_model){
    // validate inputs
    if (cm_cmmodel_init_check_values(cm_model, 0)){
        exit(EXIT_FAILURE);
    };
    // sample size
    size_t n = cm_model->d->size;
    // convert modeltype from character to index
    int modeltype_index = cm_find_modeltype_index(cm_model->modeltype);
    // compute propensity score
    vector *singleindex_score = cm_compute_singleindex_score(cm_model->x, cm_model->theta);
    vector *propscore = cm_compute_propscore(singleindex_score, modeltype_index);
    // sort propensity score
    size_t *pscore_sortindex = vector_sort_index(propscore);
    // set caliper
    size_t n1 = vector_short_sum(cm_model->d);
    size_t n0 = n - n1;
    double max_groups = cm_d_max(log(n1 + 0.0) / (n1 + 1.0), log(n0 + 0.0) / (n0 + 1.0));
    // double delta = cm_maxmin_delta(cm_model->d, propscore, pscore_sortindex);
    double delta = 0.0;
    if (cm_model->delta == 0.0){
        delta = cm_d_max(max_groups, cm_maxmin_delta(cm_model->d, propscore, pscore_sortindex));
    } else {
        delta = cm_model->delta;
    }   
    // double delta = 10 * log(n) / n;
    // load values to the model
    cm_model->n = n;
    cm_model->_modeltype_index = modeltype_index;
    cm_model->_singleindex_score = singleindex_score;
    cm_model->_propscore = propscore;
    cm_model->_pscore_sortindex = pscore_sortindex;
    cm_model->delta = delta;
    // instantiation successfull
    cm_model->_isinstantiated = 1;
    return 0;
}


CMResults *cm_cm(CMModel *cm_model){
    // ensure model has been instantiated
    if (!cm_model->_isinstantiated){
        fprintf(stderr, "%s: line %d: CMModel has not been instantiated prior to calling `cm_cm`. Do so with `cm_initialise`.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    // validate inputs
    if (cm_cmmodel_check_values(cm_model, 0)){
        exit(EXIT_FAILURE);
    };
    // instantiate model with computed propscore
    CMModelKnownPropscore *cm_known = malloc(sizeof(CMModelKnownPropscore));
    cm_known->n = cm_model->n;
    cm_known->y = cm_model->y;
    cm_known->d = cm_model->d;
    cm_known->alpha = cm_model->alpha;
    cm_known->beta = cm_model->beta;
    cm_known->kappa_a = cm_model->kappa_a;
    cm_known->kappa_gamma = cm_model->kappa_gamma;
    cm_known->_kappa_gamma_derivative = cm_model->kappa_gamma_derivative;
    cm_known->_modeltype = cm_model->modeltype;
    cm_known->_x = cm_model->x;
    cm_known->_theta = cm_model->theta;
    cm_known->_modeltype_index = cm_model->_modeltype_index;
    cm_known->_singleindex_score = cm_model->_singleindex_score;
    cm_known->propscore = cm_model->_propscore;
    cm_known->_pscore_sortindex = cm_model->_pscore_sortindex;
    cm_known->_called_internally = 1; // indicates that `cm_cm_known_propscore` is called by this function; used for variance estimation.
    cm_known->_isinstantiated = 1;  // CMModelKnownPropscore is instantiated when called by this function.
    cm_known->delta = cm_model->delta;
    CMResults *results = cm_cm_known_propscore(cm_known);
    // for (int i=0; i<y->size; i++){
    //     printf("propscore[%d] = %f\n", i, vector_get(propscore, i));
    // }
    // gsl_vector_view x_view = gsl_matrix_column(cm_model->x, 0);
    // for (int i=0; i<n; i++){
    //     printf("x_view[%d] = %f, x[%d,0] = %f\n", i, (&x_view.vector)->data[(&x_view.vector)->stride * i], i, matrix_get(x, i, 0));
    // }
    // printf("x_view[0] = %f, x[0,0] = %f\n", (&x_view.vector)->data[0], matrix_get(x, 0, 0));
    return results;
}


void test_cm(void){
    test_vector();
    test_vector_sort_index_ptr();
    test_vector_sort_index();
    test_linkedlist();
    test_dynamicarray();
    test_matrix();
    test_nonpara();
    // cm
    test_cm_int_min();
    test_cm_int_max();
    test_cm_d_min();
    test_cm_d_max();
    test_cm_inrange();
    test_cm_abs_distance();
    test_cm_cmmodel_check_values();
    test_cm_cmmodelknownpropscore_check_values();
    test_cm_compute_singleindex_score();
    test_cm_compute_propscore();
    test_cm_matches();
    test_cm_maxmin_delta();
    test_cm_te_estimator_thread();
    test_cm_update_qs();
    test_cm_update_fisher_info();
    test_cm_compute_var_estpi();
    printf("=== Tests: all tests for `cm` completed: no issues found.\n");
}