//////////////////////////////////////////////////////////////////////////
//
// * FILE: example_simple.c
// * DESCRIPTION:
//    Simple example of using caliper matching library.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 13/oct/2023
// * COMPILE:
//   Mac: gcc -pthread -lgsl -l gslcblas example_simple.c cm.c vector.c matrix.c dynamicarray.c linkedlist.c propscore.c nonpara.c linalg.c -o example_simple

//
//////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "cm.h"
#include "vector.h"
#include "matrix.h"


////////////////////////  Functions for Data Generating Process


// Generates `n` observations of `k`-dimensional covariates writing them into an `n`-by-`k` matrix `x`.
void generate_x(gsl_rng *rng, int k, matrix *x){
    // Note: all x_k must have zero mean to test ATE estimates
    for (int i=0; i<x->size1; i++){
        //vector_set(x, i, gsl_ran_gaussian(r, 1.0));
        for (int j=0; j<x->size2; j++){ // all covariates are Uniform(-3,3)
            matrix_set(x, i, j, gsl_ran_flat(rng, -10.0, 10.0));
        }
    }
}


double generate_d_error_term_logit(gsl_rng *rng){
    return gsl_ran_logistic(rng, 1.0);
}

double generate_d_error_term_probit(gsl_rng *rng){
    return gsl_ran_ugaussian(rng);
}

// Generates `n` observations of treatment writing them to vector `d`.
void generate_d(gsl_rng *rng, matrix *x, char *modeltype, vector *theta, vector_short *d){
    short d_i;
    double eps_i, score_i;
    // choose between logit or probit model
    double (*generate_error_term)(gsl_rng*);  // pointer to function
    if (!strcmp(modeltype, "logit")){  // logit model
        generate_error_term = generate_d_error_term_logit;
    } else if (!strcmp(modeltype, "probit")){ // probit model
        generate_error_term = generate_d_error_term_probit;
    } else {
        fprintf(stderr, "%s: line %d: `modeltype` '%s' is not supported for propensity score.\n", __FILE__, __LINE__, modeltype);
        exit(EXIT_FAILURE);
    }
    // generate treatment
    for (int i=0; i<d->size; i++){
        // latent score
        score_i = vector_get(theta, theta->size - 1);  // intercept
        for (int j=0; j<x->size2; j++){ // dot product <theta, x>.
            score_i += vector_get(theta, j) * matrix_get(x, i, j);
        }   
        eps_i = generate_error_term(rng);
        score_i += eps_i;
        d_i = (score_i > 0 ? 1 : 0);
        vector_short_set(d, i, d_i);
    }
}

// Generates `n` observations of outcome writing them to vector `y`.
void generate_y(gsl_rng *rng, matrix *x, vector_short *d, vector *y){
    double mu0,  mu1, y_i, y_0i, y_1i;
    mu0 = 0.0;
    mu1 = 0.0;  // imply ATE = mu1 - mu0 if x's are mean zero
    if (x->size2 != 2){
        fprintf(stderr, "%s: line %d: currently only 2-dimensional covariates are supported.\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    for (int i=0; i<y->size; i++){
        y_0i = mu0 + 0.8 * matrix_get(x, i, 0) + 1.5 * matrix_get(x, i, 1) + 0.9 * gsl_ran_ugaussian(rng);
        y_1i = mu1 + 3.1 * sin(matrix_get(x, i, 0)) - 1.8 * matrix_get(x, i, 1) + 1.1 * gsl_ran_ugaussian(rng);
        // y_0i = mu0 + 0.8 * matrix_get(x, i, 0) + 1.5 * matrix_get(x, i, 1) + 0.9 * gsl_ran_ugaussian(rng);
        // y_1i = mu1 + 1.0 * matrix_get(x, i, 0) + 1.7 * matrix_get(x, i, 1) + 1.1 * gsl_ran_ugaussian(rng);
        // y_0i = mu0 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        // y_1i = mu1 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        y_i = y_0i + vector_short_get(d, i) * (y_1i - y_0i);
        vector_set(y, i, y_i);
    }
}


// Computes the propensity score according to `modeltype` from covariates `x`, parameters `theta`, and write result to `propscore`.
void compute_propscore(matrix *x, char *modeltype, vector *theta, vector *propscore){
    double score_i;
    // choose between logit or probit model
    double (*propscore_g)(double);  // pointer to function
    if (!strcmp(modeltype, "logit")){  // logit model
        propscore_g = propscore_g_logit;
    } else if (!strcmp(modeltype, "probit")){ // probit model
        propscore_g = propscore_g_probit;
    } else {
        fprintf(stderr, "%s: line %d: `modeltype` '%s' is not supported for propensity score.\n", __FILE__, __LINE__, modeltype);
        exit(EXIT_FAILURE);
    }
    // compute propensity score
    for (int i=0; i<x->size1; i++){
        score_i = vector_get(theta, theta->size-1); // intercept
        for (int j=0; j<x->size2; j++){
            score_i += vector_get(theta, j) * matrix_get(x, i, j);
        }
        vector_set(propscore, i, propscore_g(score_i));
    }
}


////////////////////////  Usage of caliper matching library

// Illustrates basic usage of caliper matching library for two cases:
//  -- Case-KnownPS: when the propensity score is known; and
//  -- Case-EstimatedPS: when an estimated propensity score is supplied in the form of known model type (logit,probit) and estimated coefficients.
// For each of the two cases two possible ways of estimation is presented. It is recommended, however, that users invoke the `_safe` versions of the function.
int main(int argc, char *argv[]){
    // Test caliper matching functions
    test_cm();

    // Setup
    // Data generating process
    // note: theta0=(0.5,0.2,intercept=0.04), x1,x2 ~ Uniform[-3,3] gives propensity scores away from 0 & 1;
    //       theta0=(0.5,0.2,intercept=0.04), x1,x2 ~ Uniform[-10,10] given propensity scores very close to 0 & 1.
    size_t n = 5000;  // sample size
    size_t k = 2;  // number of covariates without constant term
    char *modeltype = "logit";  // propensity score model: "logit" or "probit"
    vector *theta0 = vector_alloc(k + 1); // propensity score parameter, last entry intercept
    vector_set(theta0, 0, 0.5);
    vector_set(theta0, 1, 0.2);
    vector_set(theta0, 2, 0.04);  // intercept
    double ate = 0.0;   // needs to be adjusted when DGP for y (generate_y) changes!

    // Generate data
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, time(NULL));  // gives different random numbers at each run of compiled binary
    matrix *x = matrix_alloc(n, k);
    vector_short *d = vector_short_alloc(n);
    vector *y = vector_alloc(n);
    generate_x(rng, k, x);
    generate_d(rng, x, modeltype, theta0, d);
    generate_y(rng, x, d, y);
    
    // Estimation
    // common settings
    double delta = 0.0; // caliper. Set to 0.0 for data-driven default choice (recommanded).
    int estimate_variance = 1;
    double beta = 0.0;  // power of n in bandwidth of variance estimation; if 0, default will be used
    double alpha = 0.0; // power of n in truncation sequence in variance estimation; if 0, default will be used
    double kappa_gamma = 0.0; // scale of bandwidth of variance estimation; if 0, default will be used
    double kappa_a = 0.0; // scale of truncation sequence in variance estimation; if 0, default will be used
    CMResults *cm_results = malloc(sizeof(CMResults)); // to store estimation results

    // -- Case-KnownPS
    CMModelKnownPropscore *cm_model_known_propscore = malloc(sizeof(CMModelKnownPropscore)); // Case-KnownPS
    cm_model_known_propscore->delta = delta;
    cm_model_known_propscore->estimate_variance = estimate_variance;
    cm_model_known_propscore->beta = beta;
    cm_model_known_propscore->alpha = alpha;
    cm_model_known_propscore->kappa_gamma = kappa_gamma;
    cm_model_known_propscore->kappa_a = kappa_a;
    vector *propscore = vector_alloc(n);
    compute_propscore(x, modeltype, theta0, propscore); // compute propensity score using known parameters theta0
    cm_model_known_propscore->y = y;
    cm_model_known_propscore->d = d;
    cm_model_known_propscore->propscore = propscore;
    cm_initialise_known_propscore(cm_model_known_propscore); // set caliper
    cm_results = cm_cm_known_propscore(cm_model_known_propscore); // estimation
    // print results
    printf("=========================================================\n");
    printf("Case-KnownPS: known propensity score\n");
    printf("=========================================================\n");
    // ATE
    printf("- ate_hat = %f (est. var. %f; s.e. %f)\n", cm_results->ate_hat, cm_results->var_hat_ate, sqrt(cm_results->var_hat_ate / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
           cm_results->var_hat_component_tau_ate, 100 * cm_results->var_hat_component_tau_ate / cm_results->var_hat_ate, 
           cm_results->var_hat_component_sigmapi_ate, 100 * cm_results->var_hat_component_sigmapi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_estpi_ate, 100 * cm_results->var_hat_component_estpi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_tau_ate + cm_results->var_hat_component_sigmapi_ate + cm_results->var_hat_component_estpi_ate, 100.0);
    // ATT
    printf("- att_hat = %f (est. var. %f; s.e. %f)\n", cm_results->att_hat, cm_results->var_hat_att, sqrt(cm_results->var_hat_att / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
           cm_results->var_hat_component_tau_att, 100 * cm_results->var_hat_component_tau_att / cm_results->var_hat_att, 
           cm_results->var_hat_component_sigmapi_att, 100 * cm_results->var_hat_component_sigmapi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_estpi_att, 100 * cm_results->var_hat_component_estpi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_tau_att + cm_results->var_hat_component_sigmapi_att + cm_results->var_hat_component_estpi_att, 100.0);
    // variance
    printf("Details of variance estimation:\n");
    printf("--- propscore_min = %.10f\n", cm_results->propscore_min);
    printf("--- propscore_max = %.10f\n", cm_results->propscore_max);
    printf("--- truncation lower threshold = %.10f\n", cm_results->truncation_low);
    printf("--- truncation higher threshold = %.10f\n", cm_results->truncation_high);
    printf("--- truncation sequence a_n = %.35f\n", cm_results->a_n);
    printf("--- truncation sequence negative-power alpha = %.10f\n", cm_results->alpha);
    printf("--- truncation sequence scale kappa_a = %.10f\n", cm_results->kappa_a);
    printf("--- bandwidth sequence gamma_n = %.35f\n", cm_results->gamma_n);
    printf("--- bandwidth sequence for derivatives gamma_derivative_n = %.35f\n", cm_results->gamma_derivative_n);
    printf("--- bandwidth sequence negative-power beta = %.35f\n", cm_results->beta);    
    printf("--- bandwidth sequence scale kappa_gamma = %.35f\n", cm_results->kappa_gamma);
    printf("--- bandwidth sequence scale for derivatives kappa_gamma_derivative = %.35f\n", cm_results->kappa_gamma_derivative);
    
    // -- Case-KnownPS (safe call)
    cm_results = cm_cm_known_propscore_safe(y, d, propscore, delta, estimate_variance, beta, alpha, kappa_a, kappa_gamma);
    // print results
    printf("=========================================================\n");
    printf("Case-KnownPS: known propensity score (safe call)\n");
    printf("=========================================================\n");
    // ATE
    printf("- ate_hat = %f (est. var. %f; s.e. %f)\n", cm_results->ate_hat, cm_results->var_hat_ate, sqrt(cm_results->var_hat_ate / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
           cm_results->var_hat_component_tau_ate, 100 * cm_results->var_hat_component_tau_ate / cm_results->var_hat_ate, 
           cm_results->var_hat_component_sigmapi_ate, 100 * cm_results->var_hat_component_sigmapi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_estpi_ate, 100 * cm_results->var_hat_component_estpi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_tau_ate + cm_results->var_hat_component_sigmapi_ate + cm_results->var_hat_component_estpi_ate, 100.0);
    // ATT
    printf("- att_hat = %f (est. var. %f; s.e. %f)\n", cm_results->att_hat, cm_results->var_hat_att, sqrt(cm_results->var_hat_att / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
           cm_results->var_hat_component_tau_att, 100 * cm_results->var_hat_component_tau_att / cm_results->var_hat_att, 
           cm_results->var_hat_component_sigmapi_att, 100 * cm_results->var_hat_component_sigmapi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_estpi_att, 100 * cm_results->var_hat_component_estpi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_tau_att + cm_results->var_hat_component_sigmapi_att + cm_results->var_hat_component_estpi_att, 100.0);
    // variance
    printf("Details of variance estimation:\n");
    printf("--- propscore_min = %.10f\n", cm_results->propscore_min);
    printf("--- propscore_max = %.10f\n", cm_results->propscore_max);
    printf("--- truncation lower threshold = %.10f\n", cm_results->truncation_low);
    printf("--- truncation higher threshold = %.10f\n", cm_results->truncation_high);
    printf("--- truncation sequence a_n = %.35f\n", cm_results->a_n);
    printf("--- truncation sequence negative-power alpha = %.10f\n", cm_results->alpha);
    printf("--- truncation sequence scale kappa_a = %.10f\n", cm_results->kappa_a);
    printf("--- bandwidth sequence gamma_n = %.35f\n", cm_results->gamma_n);
    printf("--- bandwidth sequence for derivatives gamma_derivative_n = %.35f\n", cm_results->gamma_derivative_n);
    printf("--- bandwidth sequence negative-power beta = %.35f\n", cm_results->beta);    
    printf("--- bandwidth sequence scale kappa_gamma = %.35f\n", cm_results->kappa_gamma);
    printf("--- bandwidth sequence scale for derivatives kappa_gamma_derivative = %.35f\n", cm_results->kappa_gamma_derivative);

    
    // -- Case-EstimatedPS
    double kappa_gamma_derivative = 0.0; // scale of bandwidth of variance estimation for estimating derivatives w.r.t. propscore parameters; if 0, default will be used
    CMModel *cm_model = malloc(sizeof(CMModel));
    cm_model->delta = delta;
    cm_model->modeltype = modeltype;
    cm_model->estimate_variance = estimate_variance;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma_derivative = kappa_gamma_derivative;
    cm_model->theta = theta0;  // supply theta0 pretending as if it were the estimated propensity score parameter
    cm_model->y = y;
    cm_model->d = d;
    cm_model->x = x;
    cm_initialise(cm_model);  // compute propscore and set caliper
    cm_results = cm_cm(cm_model);  // estimation
    // print results
    printf("=========================================================\n");
    printf("Case-EstimtedPS: estimated propensity score\n");
    printf("=========================================================\n");
    // ATE
    printf("Estimates: \n");
    printf("caliper delta = %f\n", cm_results->delta);
    printf("- ate_hat = %f (est. var. %f; s.e. %f)\n", cm_results->ate_hat, cm_results->var_hat_ate, sqrt(cm_results->var_hat_ate / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%).\n",
           cm_results->var_hat_component_tau_ate, 100 * cm_results->var_hat_component_tau_ate / cm_results->var_hat_ate, 
           cm_results->var_hat_component_sigmapi_ate, 100 * cm_results->var_hat_component_sigmapi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_estpi_ate, 100 * cm_results->var_hat_component_estpi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_tau_ate + cm_results->var_hat_component_sigmapi_ate + cm_results->var_hat_component_estpi_ate, 100.0);
    // ATT
    printf("- att_hat = %f (est. var. %f; s.e. %f)\n", cm_results->att_hat, cm_results->var_hat_att, sqrt(cm_results->var_hat_att / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
           cm_results->var_hat_component_tau_att, 100 * cm_results->var_hat_component_tau_att / cm_results->var_hat_att, 
           cm_results->var_hat_component_sigmapi_att, 100 * cm_results->var_hat_component_sigmapi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_estpi_att, 100 * cm_results->var_hat_component_estpi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_tau_att + cm_results->var_hat_component_sigmapi_att + cm_results->var_hat_component_estpi_att, 100.0);
    // variance
    printf("Details of variance estimation:\n");
    printf("--- propscore_min = %.10f\n", cm_results->propscore_min);
    printf("--- propscore_max = %.10f\n", cm_results->propscore_max);
    printf("--- truncation lower threshold = %.10f\n", cm_results->truncation_low);
    printf("--- truncation higher threshold = %.10f\n", cm_results->truncation_high);
    printf("--- truncation sequence a_n = %.35f\n", cm_results->a_n);
    printf("--- truncation sequence negative-power alpha = %.10f\n", cm_results->alpha);
    printf("--- truncation sequence scale kappa_a = %.10f\n", cm_results->kappa_a);
    printf("--- bandwidth sequence gamma_n = %.35f\n", cm_results->gamma_n);
    printf("--- bandwidth sequence for derivatives gamma_derivative_n = %.35f\n", cm_results->gamma_derivative_n);
    printf("--- bandwidth sequence negative-power beta = %.35f\n", cm_results->beta);    
    printf("--- bandwidth sequence scale kappa_gamma = %.35f\n", cm_results->kappa_gamma);
    printf("--- bandwidth sequence scale for derivatives kappa_gamma_derivative = %.35f\n", cm_results->kappa_gamma_derivative);
    
    // -- Case-EstimatedPS (safe call)
    cm_results = cm_cm_safe(y, d, x, modeltype, theta0, delta, estimate_variance, beta, alpha, kappa_a, kappa_gamma, kappa_gamma_derivative);  // estimation
    // print results
    printf("=========================================================\n");
    printf("Case-EstimtedPS: estimated propensity score (safe call)\n");
    printf("=========================================================\n");
    // ATE
    printf("Estimates: \n");
    printf("caliper delta = %f\n", cm_results->delta);
    printf("- ate_hat = %f (est. var. %f; s.e. %f)\n", cm_results->ate_hat, cm_results->var_hat_ate, sqrt(cm_results->var_hat_ate / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%).\n",
           cm_results->var_hat_component_tau_ate, 100 * cm_results->var_hat_component_tau_ate / cm_results->var_hat_ate, 
           cm_results->var_hat_component_sigmapi_ate, 100 * cm_results->var_hat_component_sigmapi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_estpi_ate, 100 * cm_results->var_hat_component_estpi_ate / cm_results->var_hat_ate,
           cm_results->var_hat_component_tau_ate + cm_results->var_hat_component_sigmapi_ate + cm_results->var_hat_component_estpi_ate, 100.0);
    // ATT
    printf("- att_hat = %f (est. var. %f; s.e. %f)\n", cm_results->att_hat, cm_results->var_hat_att, sqrt(cm_results->var_hat_att / n ));
    printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
           cm_results->var_hat_component_tau_att, 100 * cm_results->var_hat_component_tau_att / cm_results->var_hat_att, 
           cm_results->var_hat_component_sigmapi_att, 100 * cm_results->var_hat_component_sigmapi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_estpi_att, 100 * cm_results->var_hat_component_estpi_att / cm_results->var_hat_att,
           cm_results->var_hat_component_tau_att + cm_results->var_hat_component_sigmapi_att + cm_results->var_hat_component_estpi_att, 100.0);
    // variance
    printf("Details of variance estimation:\n");
    printf("--- propscore_min = %.10f\n", cm_results->propscore_min);
    printf("--- propscore_max = %.10f\n", cm_results->propscore_max);
    printf("--- truncation lower threshold = %.10f\n", cm_results->truncation_low);
    printf("--- truncation higher threshold = %.10f\n", cm_results->truncation_high);
    printf("--- truncation sequence a_n = %.35f\n", cm_results->a_n);
    printf("--- truncation sequence negative-power alpha = %.10f\n", cm_results->alpha);
    printf("--- truncation sequence scale kappa_a = %.10f\n", cm_results->kappa_a);
    printf("--- bandwidth sequence gamma_n = %.35f\n", cm_results->gamma_n);
    printf("--- bandwidth sequence for derivatives gamma_derivative_n = %.35f\n", cm_results->gamma_derivative_n);
    printf("--- bandwidth sequence negative-power beta = %.35f\n", cm_results->beta);    
    printf("--- bandwidth sequence scale kappa_gamma = %.35f\n", cm_results->kappa_gamma);
    printf("--- bandwidth sequence scale for derivatives kappa_gamma_derivative = %.35f\n", cm_results->kappa_gamma_derivative);
    return 0;
}