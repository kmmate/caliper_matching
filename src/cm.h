/* Header for cm.c */

#ifndef CM_H_
#define CM_H_

// includes
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

#include "vector.h"
#include "matrix.h"
#include "linalg.h"
#include "linkedlist.h"
#include "dynamicarray.h"
#include "propscore.h"
#include "nonpara.h"

// variables
#define NUM_THREADS  8
#define CM_TEST_DOUBLE_TOLERANCE_NPOWER 10
extern const int cm_number_of_modeltypes;

 // structs
typedef struct CMResults {
    // outputs
    double ate_hat;
    double att_hat;
    double var_hat_ate;  // estimated variance of ate_hat
    double var_hat_att;  // estimated variance of att_hat
    double var_hat_component_tau_ate;  // estimated V_tau
    double var_hat_component_tau_att;  // estimated V_tau
    double var_hat_component_sigmapi_ate;  // estimated V_sigmapi
    double var_hat_component_sigmapi_att;  // estimated V_sigmapi
    double var_hat_component_estpi_ate;  // estimated variance component deriving from the estimation of prepensity score
    double var_hat_component_estpi_att;  // estimated variance component deriving from the estimation of prepensity score
    dynamicarray_int **match_indices; // match_indices[i] is a dynamicarray_int holding the indices of `i`s matches
    vector_int *number_of_matches;  // row sums of matchmatrix
    vector *w;
    // inputs
    size_t n;
    char *modeltype; // propensity score modeltype
    vector *theta_hat;  // propensity score (estimated) parameter
    double caliper;
} CMResults;

typedef struct CMModel {
    size_t n;   // sample size 
    vector *y;  // outcome variable
    vector_short *d;  // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise
    matrix *x;  // design matrix
    char *modeltype;    // propensity score model type
    vector *theta;  // propensity score model parameter, of length  (number of columns in `x`)+1
    double delta;  // caliper
    double beta;    // negative-exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
    double alpha;   // negative-exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
    // Following fields are for internal use. Any user-provided values will be overwritten.
    int _isinstantiated;    // indicates whether the model has been instantiated.
    int _modeltype_index;  // index of propensity score model type. Optional; supplied internally
    vector *_singleindex_score;  // single-index x*theta. Optional; supplied internally.
    vector *_propscore;  // vector of propensity score values. Optional; supplied internally.
    size_t *_pscore_sortindex;  // optional; a permutation of indices {0,...,n-1} which sorts `propscore` into ascending order. If not given, it's computed internally.
} CMModel;


typedef struct CMModelKnownPropscore {
    size_t n;   // sample size
    vector *y;  // outcome variable
    vector_short *d;  // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise
    vector *propscore;  // vector of propensity score values
    char *modeltype;    // propensity score model type
    double delta;  // caliper
    double beta;    // negative exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
    double alpha;   // negative exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
    // Following fields are for internal use. Any user-provided values will be overwritten.
    matrix *_x;  // covariate matrix. Optional; supplied internally
    int _modeltype_index;  // index of propensity score model type. Optional; supplied internally
    vector *_singleindex_score;  // single-index x*theta. Optional; supplied internally.
    size_t *_pscore_sortindex;  // optional; a permutation of indices {0,...,n-1} which sorts `propscore` into ascending order. If not given, it's computed internally.
    vector *_theta;     // propensity score parameter. Optional; supplied internally
    int _called_internally;  // optional; overwritten internally if need be
    int _isinstantiated;    // indicates whether the model has been instantiated.
} CMModelKnownPropscore;

// functions
int cm_initialise(CMModel *cmmodel);
CMResults *cm_cm(CMModel *cmmodel);


// tests
void test_cm(void);


#endif // CM_H_