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
#define CM_VERBOSE_MODE 0  // printing for debugging
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
    double var_hat_component_estpi_ate;  // estimated variance component deriving from the estimation of propensity score
    double var_hat_component_estpi_att;  // estimated variance component deriving from the estimation of propensity score
    dynamicarray_int **match_indices; // match_indices[i] is a dynamicarray_int holding the indices of `i`s matches
    vector_int *number_of_matches;  // row sums of matchmatrix
    // inputs
    size_t n;
    char *modeltype; // propensity score modeltype
    vector *theta_hat;  // propensity score (estimated) parameter
    double delta; // caliper that is actually used. If zero was passed, then it is equal to the default data-driven value. If a positive value was passed, then it is equal to that instead.
    int estimate_variance; // Whether variance estimation was performed (`estimate_variance`), or not (`!estimate_variance`).
    // fields related to nonparametric variance estimation
    double a_n; // value of truncation sequence in nonparametric variance estimation. Equal to `kappa_a` * `n` ^ `alpha`.
    double gamma_n; // value of bandwidth in nonparametric variance estimation. Equal to `kappa_gamma` * `n` ^ `beta`.
    double gamma_derivative_n; // value of bandwidth in nonparametric variance estimation used in derivative estimation w.r.t. propensity score parameters; it is zero when the propensity score is known. Equal to `kappa_gamma_derivative` * `n` ^ `beta`.
    double beta;  // negative-exponent of bandwidth in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
    double alpha;  // negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
    double kappa_gamma;  // scale of bandwidth sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
    double kappa_a;  // scale of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value, otherwise it is equal to the passed value.
    double kappa_gamma_derivative;  // negative-exponent of truncation sequence in nonparametric variance estimation that is actually used. If zero was passed, it is the default value (which is zero for known propensity scores), otherwise it is equal to the passed value.
    double propscore_min; // the smallest propensity score value
    double propscore_max; // the largest propensity score value
    double truncation_low; // lower threshold for variance estimation = `propscore_min`+`a_n`
    double truncation_high; // higher threshold for variance estimation = `propscore_max` - `a_n`
} CMResults;

typedef struct CMModel {
    size_t n;   // sample size 
    vector *y;  // outcome variable
    vector_short *d;  // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise
    matrix *x;  // design matrix
    char *modeltype;    // propensity score model type
    vector *theta;  // propensity score model parameter, of length  (number of columns in `x`)+1
    double delta; // caliper. If zero is passed, then the default data-driven value is used (recommanded). If a positive value is passed, that is used instead.
    int estimate_variance; // If zero is passed, the variances are not estimated, but set to zero automatically. This gains a speed-up when variance estimates aren not required.
    double beta;    // negative-exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
    double alpha;   // negative-exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
    double kappa_a; // scale parameter of truncation sequence in nonparametric variance estimation. If zero is passed, a default value is used.
    double kappa_gamma; // scale parameter of bandwidth in nonparametric variance estimation. If zero is passed, a default value is used. 
    double kappa_gamma_derivative; // scale parameter of bandwidth in nonparametric variance estimation estimating derivatives w.r.t. theta. If zero is passed, a default value is used. 
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
    double delta; // caliper. If zero is passed, then the default data-driven value is used (recommanded). If a positive value is passed, that is used instead.
    int estimate_variance; // If zero is passed, the variances are not estimated, but set to zero automatically. This gains a speed-up when variance estimates aren not required.
    double beta;    // negative exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
    double alpha;   // negative exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
    double kappa_a; // scale parameter of truncation sequence in nonparametric variance estimation. If zero is passed, a default value is used.
    double kappa_gamma; // scale parameter of bandwidth in nonparametric variance estimation. If zero is passed, a default value is used. 
    // Following fields are for internal use. Any user-provided values will be overwritten.
    char *_modeltype;    // propensity score model type
    matrix *_x;  // covariate matrix. Optional; supplied internally
    int _modeltype_index;  // index of propensity score model type. Optional; supplied internally
    vector *_singleindex_score;  // single-index x*theta. Optional; supplied internally.
    size_t *_pscore_sortindex;  // optional; a permutation of indices {0,...,n-1} which sorts `propscore` into ascending order. If not given, it's computed internally.
    vector *_theta;     // propensity score parameter. Optional; supplied internally
    int _called_internally;  // optional; overwritten internally if need be
    int _isinstantiated;    // indicates whether the model has been instantiated.
    double _kappa_gamma_derivative; // scale parameter of bandwidth in nonparametric variance estimation estimating derivatives w.r.t. theta. If zero is passed, a default value is used. 
} CMModelKnownPropscore;

// functions
int cm_initialise(CMModel *cmmodel);
CMResults *cm_cm(CMModel *cmmodel);
CMResults *cm_cm_safe(vector *y, // outcome variable.
                      vector_short *d, // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise.
                      matrix *x, // `n`-by-`k` matrix of covariates.
                      char *modeltype, // propensity score modeltype
                      vector *theta, // the estimated propensity score parameter.
                      double delta, // caliper.
                      int estimate_variance,  // If zero is passed, the variances are not estimated, but set to zero automatically. This gains a speed-up when variance estimates are not required.
                      double beta,    // negative-exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
                      double alpha,   // negative-exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
                      double kappa_a, // scale parameter of truncation sequence in nonparametric variance estimation. If zero is passed, a default value is used.
                      double kappa_gamma, // scale parameter of bandwidth in nonparametric variance estimation. If zero is passed, a default value is used. 
                      double kappa_gamma_derivative // scale parameter of bandwidth in nonparametric variance estimation estimating derivatives w.r.t. theta. If zero is passed, a default value is used. 
);
int cm_initialise_known_propscore(CMModelKnownPropscore *cmmodel);
CMResults *cm_cm_known_propscore(CMModelKnownPropscore *cmmodel);
CMResults *cm_cm_known_propscore_safe(vector *y, // outcome variables
                                      vector_short *d, // treatment indicator: an entry `i` is one if unit `i` is treated, zero otherwise   
                                      vector *propscore, // vector of propensity score values
                                      double delta, // caliper
                                      int estimate_variance,  // If zero is passed, the variances are not estimated, but set to zero automatically. This gains a speed-up when variance estimates are not required.
                                      double beta,    // negative-exponent of bandwidth in nonparametric variance estimation. If zero is passed, a dafault value is used.
                                      double alpha,   // negative-exponent of truncation sequence in nonparemetric variance estimation. If zero is passed, a default value is used.
                                      double kappa_a, // scale parameter of truncation sequence in nonparametric variance estimation. If zero is passed, a default value is used.
                                      double kappa_gamma // scale parameter of bandwidth in nonparametric variance estimation. If zero is passed, a default value is used. 
);


// tests
void test_cm(void);


#endif // CM_H_