//////////////////////////////////////////////////////////////////////////
//
// * FILE: example_simulation.c
// * DESCRIPTION:
//    Example of using caliper matching library in a Monte Carlo simulation.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 13/nov/2023
// * COMPILE:
//   Mac: gcc -pthread -lgsl -l gslcblas example_simulation.c cm.c vector.c matrix.c dynamicarray.c linkedlist.c propscore.c nonpara.c linalg.c -o example_simulation

//
//////////////////////////////////////////////////////////////////////////

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

// Generate truncated-normally distributed random variable between `a` and `b` parametrized by `mu` and `sigma`.
double ran_truncgaussian(gsl_rng *rng, double mu, double sigma, double a, double b){
    double u = gsl_ran_flat(rng, 0, 1);
    double alpha = (a - mu + 0.0) / sigma;
    double beta = (b - mu + 0.0) / sigma;
    double x = gsl_cdf_ugaussian_Pinv(gsl_cdf_ugaussian_P(alpha) + u * (gsl_cdf_ugaussian_P(beta) - gsl_cdf_ugaussian_P(alpha)));
    x = x * sigma + mu;
    return x;
}


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



int main(int argc, char *argv[]){
    // Test caliper matching functions
    test_cm();
    
    // Setup
    // Data generating process
    // note: theta0=(0.5,0.2,intercept=0.04), x1,x2 ~ Uniform[-3,3] gives propensity scores away from 0 & 1;
    //       theta0=(0.5,0.2,intercept=0.04), x1,x2 ~ Uniform[-10,10] given propensity scores very close to 0 & 1.
    size_t n = 50000;  // sample size
    size_t k = 2;  // number of covariates without constant term
    char *modeltype = "logit";  // propensity score model: "logit" or "probit"
    vector *theta0 = vector_alloc(k + 1); // propensity score parameter, last entry intercept
    vector_set(theta0, 0, 0.5);
    vector_set(theta0, 1, 0.2);
    vector_set(theta0, 2, 0.04);  // intercept
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng, time(NULL));  // gives different random numbers at each run of compiled binary
    matrix *x = matrix_alloc(n, k);
    vector_short *d = vector_short_alloc(n);
    vector *y = vector_alloc(n);
    double ate = 0.0;   // needs to be adjusted when DGP for y (generate_y) changes!
    // Caliper matching estimator
    double beta = 0.0;  // power of n in bandwidth of variance estimation; if 0, default will be used
    double alpha = 0.0; // power of n in truncation sequence in variance estimation; if 0, default will be used
    double kappa_gamma = 0.0; // scale of bandwidth of variance estimation; if 0, default will be used
    double kappa_a = 0.0; // scale of truncation sequence in variance estimation; if 0, default will be used
    double kappa_gamm_derivative = 0.0; // scale of bandwidth of variance estimation for estimating derivatives w.r.t. propscore parameters; if 0, default will be used
    CMModel *cm_model = malloc(sizeof(CMModel));
    cm_model->modeltype = modeltype;
    cm_model->theta = theta0;
    cm_model->beta = beta;
    cm_model->alpha = alpha;
    cm_model->kappa_gamma = kappa_gamma;
    cm_model->kappa_a = kappa_a;
    cm_model->kappa_gamma_derivative = kappa_gamm_derivative;
    CMResults *cm_results = malloc(sizeof(CMResults));
    
    // Simulation
    int reps = 1;
    double mean_caliper = 0.0;
    int count_ate_hat_above_ate = 0;
    double percentage_obs_without_matches = 0.0;
    double mean_ate_hat = 0.0;
    time_t start_time, end_time;
    time(&start_time);
    for (int rep = 0; rep<reps; rep++){
        // generate data
        generate_x(rng, k, x);
        generate_d(rng, x, modeltype, theta0, d);
        generate_y(rng, x, d, y);
        cm_model->y = y;
        cm_model->d = d;
        cm_model->x = x;
        // perform estimation
        cm_initialise(cm_model);  // compute propscore and set caliper
        // printf("delta = %f.\n", cm_model->delta);
        // cm_model->delta *= 10;  // re-scaling
        mean_caliper += cm_model->delta / reps;
        cm_results = cm_cm(cm_model);
        // obtain metrics from estimation results
        mean_ate_hat += cm_results->ate_hat / (reps + 0.0);
        count_ate_hat_above_ate += cm_results->ate_hat >= ate;
        for (int i=0; i<n; i++){
            percentage_obs_without_matches += 100 * (vector_int_get(cm_results->number_of_matches, i) == 0) / (n + 0.0);
        }
        if (reps == 1){
            printf("ate_hat = %f (est. var. %f; s.e. %f)\n", cm_results->ate_hat, cm_results->var_hat_ate, sqrt(cm_results->var_hat_ate / n ));
            printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
                   cm_results->var_hat_component_tau_ate, 100 * cm_results->var_hat_component_tau_ate / cm_results->var_hat_ate, 
                   cm_results->var_hat_component_sigmapi_ate, 100 * cm_results->var_hat_component_sigmapi_ate / cm_results->var_hat_ate,
                   cm_results->var_hat_component_estpi_ate, 100 * cm_results->var_hat_component_estpi_ate / cm_results->var_hat_ate,
                   cm_results->var_hat_component_tau_ate + cm_results->var_hat_component_sigmapi_ate + cm_results->var_hat_component_estpi_ate, 100.0);
            printf("att_hat = %f (est. var. %f; s.e. %f)\n", cm_results->att_hat, cm_results->var_hat_att, sqrt(cm_results->var_hat_att / n ));
            printf("--- variance components: \n var_tau = %.4f (%.2f %%), \n var_sigmapi = %.4f (%.2f %%), \n var_estpi = %.4f (%.2f %%), \n total = %.4f (%.2f %%)\n",
                   cm_results->var_hat_component_tau_att, 100 * cm_results->var_hat_component_tau_att / cm_results->var_hat_att, 
                   cm_results->var_hat_component_sigmapi_att, 100 * cm_results->var_hat_component_sigmapi_att / cm_results->var_hat_att,
                   cm_results->var_hat_component_estpi_att, 100 * cm_results->var_hat_component_estpi_att / cm_results->var_hat_att,
                   cm_results->var_hat_component_tau_att + cm_results->var_hat_component_sigmapi_att + cm_results->var_hat_component_estpi_att, 100.0);
        }
    }
    time(&end_time);
    printf("Sample size n = %zu; n/NUM_THREADS = %lu; reps = %d; elapsed time: %.10lf minutes.\n", n, n/NUM_THREADS, reps, difftime(end_time, start_time)/60.0);
    printf("Mean caliper = %.10f.\n", mean_caliper);
    printf("Mean(ate_hat) = %.3f\n", mean_ate_hat);
    printf("Percentage(ate_hat is above true ate across reps) = %.1f %%\n", 100 * count_ate_hat_above_ate / (reps + 0.0));
    printf("Mean(percentage(observations without matches)) = %.10f %%\n", percentage_obs_without_matches / (reps + 0.0));
    // Benchmarks for sample size without variance estimation:
    // - maximum tested that still runs for k = 2 is 300_000_000 (15mins on 1 thread, 30mins on 8 thread)
    // - (n,k)=(200_000_000,2) 9.3mins on 1 thread, 7.7mins on 8 thread.
    // - (n,k)=(10_000_000,2) 0.36mins on 1 thread, 0.16mins on 8 thread
    /* Examples runs performed: 
    Note: at the beginning of each line, "-" = has been run at least twice, no perceived issues 
                                               (i.e. all realised values are theoretically possible);
                                         "~" = has been run at least twice, no perceived issues but poor performance
                                               (i.e. te point estimates are off, large s.e.)
                                         "x" = have not been run yet
                                         "!" = problem
        - NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        - NUM_THREADS = 8, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;        
        =========
        ~ NUM_THREADS = 1, n = 100, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        - NUM_THREADS = 1, n = 1000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        - NUM_THREADS = 1, n = 20000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        ~ NUM_THREADS = 1, n = 10000, ate=0, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        ~ NUM_THREADS = 1, n = 10000, ate=2, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        - NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-3,3], y1coeffx1 = 30, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        ~[^1] NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 3, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        -[^2] NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 300, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        ![^3] NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.001, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        -[^4] NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 0.1, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        - NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 1, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        - NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 10, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        - NUM_THREADS = 1, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 30, kappa_gamma = 100, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        -[^5] NUM_THREADS = 8, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 3, kappa_gamma = 1, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        -[^6] NUM_THREADS = 8, n = 10000, ate=20, x1,x2~U[-10,10], y1coeffx1 = 3, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;
        =========
        -[^7] NUM_THREADS = 8, n = 10000, ate=-800, x1,x2~U[-3,3], y1coeffx1 = 3, kappa_gamma = 0.01, npower_gamma = 1.0/5, kappa_a = 0.00001, npower_gamma = 1/5.5;

        ========= Notes:
        [^1]: ate estimates off, true ate not in CI95%
        [^2]: ate estimates vary a lot, but CI95% contains true ate 
        [^3]: ate estimates off, infinite var
        [^4]: ate estimates vary a lot (covarge of CI95% might be lower than 95% ?)
        [^5]: ate estimate always around 16-17. Why? S.e. estimate stable around 1; CI95% bad coverage.
        [^6]: similar to [^5], but s.e. estimate gives the good covarage.
    */

    /*
    Caliper choice Monte Carlo results for reps = 1000.
    - DGP1:
        Covariates: x_1, x_2 ~ N(0,1) i.i.d.;
        Propensity score: known (= nonestimated) logit, values close to 0 and 1, theta0_0=0.5, theta0_1=0.2, theta0_2=intercept=0.04;
        Potential outcomes:
        mu0 = 0
        mu1 = 10
        y_0i = mu0 + 0.8 * matrix_get(x, i, 0) + 1.5 * matrix_get(x, i, 1) + 0.9 * gsl_ran_ugaussian(rng)
        y_1i = mu1 + 1.0 * matrix_get(x, i, 0) + 1.7 * matrix_get(x, i, 1) + 1.1 * gsl_ran_ugaussian(rng)
        ==> ate = 10.0;
        maxmindelta = max_{i\in[n]}min_{j:D_j!=D_i}|propscore_j-propscore_i|
        maxmindelta2 := max(maxmindelta, max(log n0 / n0, log n1 / n1)), where n_d = |{i:D_i=d}|. 
        
        Caliper         n        mean(caliper)   mean(ATEhat)    mean_percentage(ATE<ATEhat)     mean(percentage(observations without matches))
        ---------------------------------------------------------------------------------------------------------------------------------------
        maxmindelta     10,000    0.0050077485      10.23            73.8 %                       0 %  
                        50,000    0.0017760549      10.06            66.6 %                       0 %
                        80,000    0.0013203142      10.04            60.6 %                       0 %   
        
        maxmindelta2    10,000    0.0049577912      10.22            71.4 %                       0 %  
                        50,000    0.0017991784      10.07            66.4 %                       0 %  
                        80,000    0.0013496377      10.04            63.1 %                       0 %  
        
        log(n)/n        10,000    0.0009210340      8.25              0.0 %                      17.61 %
                        50,000    0.0002163956      8.61              0.0 %                      13.91 %
                        80,000    0.0001411223      8.71              0.0 %                      13.0 %   
        
        5*log(n)/n      10,000    0.0046051702     10.08             60.7 %                      0.86 %
                        50,000    0.0010819778      9.95             39.5 %                      0.66 %
                        80,000    0.0007056114      9.94             33.8 %                      0.65 %   
        
        10*log(n)/n     10,000    0.0092103404      10.4             87.1 %                      0.06 %
                        50,000    0.0021639557      10.06            63.8 %                      0.12 %
                        80,000    0.0014112227      10.02            57.4 %                      0.10 %   
    
    - DGP2:
        Covariates: x_1, x_2 ~ N(0,1) i.i.d.;
        Propensity score: known (= nonestimated) logit, values close to 0 and 1, theta0_0=0.5, theta0_1=0.2, theta0_2=intercept=0.04;
        Potential outcomes:
        mu0 = 0
        mu1 = 10
        y_0i = mu0 + 0.8 * matrix_get(x, i, 0) + 1.5 * matrix_get(x, i, 1) + 0.9 * gsl_ran_ugaussian(rng);
        y_1i = mu1 + 3.1 * sin(matrix_get(x, i, 0)) - 1.8 * matrix_get(x, i, 1) + 1.1 * gsl_ran_ugaussian(rng);
        ==> ate = 10.0;
        maxmindelta = max_{i\in[n]}min_{j:D_j!=D_i}|propscore_j-propscore_i|
        maxmindelta2 := max(maxmindelta, max(log n0 / n0, log n1 / n1)), where n_d = |{i:D_i=d}|.

        Caliper         n        mean(caliper)   mean(ATEhat)    mean_percentage(ATE<ATEhat)     mean(percentage(observations without matches))
        ---------------------------------------------------------------------------------------------------------------------------------------
        maxmindelta     10,000    0.0050008412      10.006             51.9 %                          0 %
                        50,000    0.0017979590       9.99              49.0 %                          0 %
                        80,000    0.0013569337      10.00              50.5 %                          0 %    

        maxmindelta2    10,000    0.0051208164      10.01              52.8 %                          0 %
                        50,000    0.0017918744      10.00              48.5 %                          0 %
                        80,000    0.0013672854       9.99              49.1 %                          0 %
        
        log(n)/n        10,000    0.0009210340       8.22               0.3 %                          17.72 %
                        50,000    0.0002163956       8.61               0.0 %                          13.99 %
                        80,000    0.0001411223       8.70               0.0 %                          12.99 % 

        5*log(n)/n      10,000    0.0046051702       9.91               47.6 %                          0.87 %
                        50,000    0.0010819778       9.93               43.5 %                          0.72 %
                        80,000    0.0007056114       9.93               41.3 %                          0.64 %

        10*log(n)/n     10,000    0.0092103404      10.01               50.2 %                          0.05 %
                        50,000    0.0021639557      10.00               49.5 %                          0.10 %
                        80,000    0.0014112227       9.98               46.3 %                          0.12 %
    */

   /*
   Variance estimation: Monte Carlo experiments with very few repetitions to assess the performance of variance estimators.
   Aims: 
   1. See if variance estimators converge to their true values, which are known (only) for a few special cases.
      The special cases considered here are when mu1(pi(x))-mu0(pi(x))=c for a constant c. 
      These cases imply ate=att=c and var_tau=var_estpi=0 and var_sigmapi!=0.
   2. Guidance on how to choose the bandwidth and truncation sequences.
    - (1) DGP:
        Covariates: x_1, x_2 ~ N(0,1) i.i.d.;
        Propensity score: known (= nonestimated) logit, values close to 0 and 1, theta0_0=0.5, theta0_1=0.2, theta0_2=intercept=0.04;
        Potential outcomes:
        mu0 = 0
        mu1 = 0 
        y_0i = mu0 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        y_1i = mu1 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        ==> ate=att=0; mu0(pi(x))-mu1(pi(x))=0 for all x, so var_tau and var_estpi should both be zero for both ate and att.
     Settings:
        Caliper = maxmindelta2;
        Truncation:  kappa_a = 0 (=no truncation);
        Bandwidth:   npower_gamma = 1/4.000000001;
     Results from a few repetitions:
        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate   var_tau             var_sigmapi        var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        0.5           1,000      0       0.0889    1        ate         0.0038      0.0001 (0.14 %)    0.0427 (73.31 %)   0.0155 (26.55 %)  0.0583      0.0076
                                                            att        -0.0045     -0.0000 (!!!!!!)    0.1205 (!!!!!!)    0.0342 (!!!!!!)   !!!!!!      !!!!!!
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.0069      0.0000 (0.07 %)    0.0350 (87.65 %)   0.0049 (12.28 %)  0.0400      0.0063
                                                            att        -0.0001      0.0001 (0.05 %)    0.0897 (88.03 %)   0.0122 (11.92 %)  0.1019      0.0101  
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        0.0073       0.0001 (0.29 %)    0.0237 (84.11 %)   0.0044 (15.6 %)   0.0282      0.0053
                                                            att        0.0086       0.0001 (0.18 %)    0.0380 (98.46 %)   0.0005 ( 1.37 %)  0.0386      0.0062
                    ==========================================================================================================================================
                      5,000      0      0.0595     1        ate        0.0023       0.0001 (0.23 %)    0.0382 (92.31 %)   0.0031 ( 7.46 %)  0.0414      0.0029   
                                                            att        0.0035       0.0002 (0.34 %)    0.0588 (92.3 %)    0.0047 ( 7.36 %)  0.0637      0.0036 
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate        0.0002       0.0000 (0.04 %)    0.0444 (79.6 %)    0.0113 (20.36 %)  0.0557      0.0034   
                                                            att        0.0001       0.0001 (0.05 %)    0.1050 (90.77 %)   0.0106 ( 9.19 %)  0.1157      0.0048  
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        0.0017       0.0000 (0.14 %)    0.0318 (98.68 %)   0.0004 ( 1.18 %)  0.0322      0.0025
                                                            att        0.0002       0.0000 (0.05 %)    0.0587 (96.2 %)    0.0023 ( 3.76 %)  0.0610      0.0035
                    ==========================================================================================================================================
                      50,000     0      0.0334     1        ate        0.0000       0.0000 (0.08 %)    0.0354 (97.12 %)   0.0010 ( 2.8 %)   0.0365      0.0009   
                                                            att        0.0002       0.0001 (0.08 %)    0.0692 (98.35 %)   0.0011 ( 1.57 %)  0.0703      0.0012 
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate       -0.0001       0.0000 (0.06 %)    0.0356 (94.12 %)   0.0022 ( 5.83 %)  0.0378      0.0009    
                                                            att        0.0000       0.0001 (0.07 %)    0.0739 (97.17 %)   0.0021 ( 2.75 %)  0.0761      0.0012
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        0.0006       0.0001 (0.13 %)    0.0380 (96.28 %)   0.0014 ( 3.59 %)  0.0395      0.0009 
                                                            att        0.0011       0.0001 (0.15 %)    0.0819 (97.6 %)    0.0019 ( 2.24 %)  0.0839      0.0013   
                    ==========================================================================================================================================

        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate   var_tau             var_sigmapi        var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        5.0           1,000      0       0.8891    1        ate         0.0032      0.0004 (0.73 %)    0.0423 (83.71 %)   0.0079 (15.56 %)  0.0505      0.0071
                                                            att         0.0008      0.0007 (0.79 %)    0.0833 (90.09 %)   0.0084 ( 9.12 %)  0.0924      0.0100
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.0029      0.0005 (0.82 %)    0.0533 (82.09 %)   0.0111 (17.09 %)  0.0649      0.0086
                                                            att         0.0142      0.0007 (0.54 %)    0.1173 (93.24 %)   0.0078 ( 6.22 %)  0.1257      0.0112
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        -0.0002      0.0006 (0.99 %)    0.0491 (87.28 %)   0.0066 (11.74 %)  0.0563      0.0075 
                                                            att        -0.0046      0.0010 (0.93 %)    0.0990 (92.17 %)   0.0074 ( 6.9 %)   0.1074      0.0104 
                    ==========================================================================================================================================
                      5,000      0       0.5946    1        ate         0.0021      0.0002 (0.44 %)    0.0493 (94.5 %)    0.0026 (5.06 %)   0.0521      0.0032
                                                            att         0.0018      0.0004 (0.43 %)    0.0980 (96.4 %)    0.0032 (3.18 %)   0.1017      0.0045
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.0006      0.0002 (0.52 %)    0.0459 (94.25 %)   0.0026 (5.23 %)   0.0487      0.0031
                                                            att         0.0004      0.0005 (0.55 %)    0.0879 (96.67 %)   0.0025 (2.78 %)   0.0908      0.0043
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         0.0042      0.0003 (0.54 %)    0.0472 (92.76 %)   0.0034 (6.69 %)   0.0509      0.0032
                                                            att         0.0033      0.0006 (0.62 %)    0.0917 (95.04 %)   0.0042 (4.35 %)   0.0965      0.0044  
                    ==========================================================================================================================================
                     50,000      0       0.3344    1        ate         0.0009      0.0001 (0.26 %)    0.0439 (99.09 %)   0.0003 (0.65 %)   0.0443      0.0009
                                                            att         0.0007      0.0002 (0.25 %)    0.0858 (99.3 %)    0.0004 (0.45 %)   0.0864      0.0013
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate        -0.0002      0.0001 (0.25 %)    0.0444 (98.92 %)   0.0004 (0.83 %)   0.0449      0.0009
                                                            att        -0.0010      0.0002 (0.24 %)    0.0904 (99.29 %)   0.0004 (0.47 %)   0.0911      0.0014
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         0.0002      0.0001 (0.24 %)    0.0453 (99.12 %)   0.0003 (0.65 %)   0.0457      0.0010                  
                                                            att        -0.0006      0.0002 (0.22 %)    0.0895 (99.46 %)   0.0003 (0.32 %)   0.0900      0.0013
     Findings: 1. ate_hat att_hat both converge to ate as n increases.
            For ate:
               2. The (proportion of the estimated var_tau to var_total) decreases to the true value zero as n increases for both choices of kappa_gamma, but:
                    2.1 the decrease is faster for kappa_gamma=5;
                    2.2 it reaches values closer to zero percent for kappa_gamma=0.5, however.
               3. Point 2. applies also to the absolute value of var_tau, except for Point 2.2.
               4. The (proportion of the estimated var_tau to var_total) varies slightly less for kappa_gamma=5, than for kappa_gamma=0.5.
               5. The (proportion of the estimated var_estpi to var_total) decreases to the true value zero as n increases for both choices of kappa_gamma, but:
                    5.1 the decrease is faster for kappa_gamma=5 and also reaches values closer to zero percent for kappa_gamma=5.
               6. Point 5. applies to the absolute value of var_estpi too.
            For att:
               7. The (proportion of estimated var_tau to var_total) decreases to the true value zero steadily for kappa_gamma=5, but:
                    7.1 for kappa_gamma=0.5 it seems to just vary around a level close to zero percentage, which is however lower than for kappa_gamma=5. 
               8. The absolute value of estimated var_tau is higher (=worse) for kappa_gamma=5 than kappa_gamma=0.5, but varies less. 
               9. For var_estpi, kappa_gamma=5 is better than kappa_gamma=0.5 both in terms of variation and decay to zero both in proportion and absolute level.
               10. Point 5. applies to att too. 
            Both:
               11. var_total seems to converge to the same values (ate:0.03-0.04, att:0.08-0.09) for both kappa_gamma. 
    - (2) DGP:
        Same as (1) DGP:
        Covariates: x_1, x_2 ~ N(0,1) i.i.d.;
        Propensity score: known (= nonestimated) logit, values close to 0 and 1, theta0_0=0.5, theta0_1=0.2, theta0_2=intercept=0.04;
        Potential outcomes:
        mu0 = 0
        mu1 = 0
        y_0i = mu0 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        y_1i = mu1 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        ==> ate=att=0; mu0(pi(x))-mu1(pi(x))=0 for all x, so var_tau and var_estpi should both be zero for both ate and att.
     Settings:
        Caliper = maxmindelta2;
        Truncation:  kappa_a = 10^(-10);
                     npower_a = 1/4.000000002
        Bandwidth:   npower_gamma = 1/4.000000001;
     Results from a few repetitions:
        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate   var_tau             var_sigmapi        var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        0.5           1,000  2*10^(-11)  0.0889    1        ate         0.0082      0.0001 (0.29 %)     0.0314 (78.40 %)  0.0085 (21.31 %)  0.0401      0.0063
                                                            att         0.0033      0.0003 (0.23 %)     0.0715 (62.82 %)  0.0421 (36.95 %)  0.1138      0.0107
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.0046      0.0001 (0.40 %)     0.0305 (95.31 %)  0.0014 ( 4.29 %)  0.0320      0.0057
                                                            att         0.0147      0.0001 (0.14 %)     0.0763 (95.77 %)  0.0033 ( 4.09 %)  0.0796      0.0089
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         0.0085      0.0001 (0.27 %)     0.0233 (92.11 %)  0.0019 ( 7.62 %)  0.0253      0.0050  
                                                            att         0.0074      0.0002 (0.34 %)     0.0496 (83.82 %)  0.0094 (15.84 %)  0.0592      0.0077
                    ==========================================================================================================================================
                      5,000  1*10^(-11)  0.0595    1        ate        -0.0010      0.0000 (0.03 %)     0.0409 (94.08 %)  0.0026 ( 5.89 %)  0.0435      0.0030      
                                                            att         0.0011      0.0000 (0.03 %)     0.0953 (96.46 %)  0.0035 ( 3.51 %)  0.0988      0.0044
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.0018      0.0000 (0.12 %)     0.0361 (95.93 %)  0.0015 ( 3.94 %)  0.0376      0.0027 
                                                            att         0.0021      0.0001 (0.19 %)     0.0596 (98.97 %)  0.0005 ( 0.84 %)  0.0602      0.0035   
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         0.0040      0.0001 (0.33 %)     0.0275 (97.37 %)  0.0006 ( 2.30 %)  0.0282      0.0024
                                                            att         0.0007      0.0001 (0.13 %)     0.0533 (99.45 %)  0.0002 ( 0.42 %)  0.0536      0.0033  
                    ==========================================================================================================================================
                      50,000  7*10^(-12)  0.0334   1        ate         0.0003      0.0000 (0.10 %)     0.0351 (93.40 %)  0.0024 ( 6.50 %)  0.0376      0.0009  
                                                            att        -0.0004      0.0001 (0.08 %)     0.0717 (97.84 %)  0.0015 ( 2.08 %)  0.0733      0.0012 
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate        -0.0000      0.0000 (0.08 %)     0.0359 (95.83 %)  0.0015 ( 4.09 %)  0.0374      0.0009 
                                                            att         0.0001      0.0001 (0.07 %)     0.0692 (97.35 %)  0.0018 ( 2.57 %)  0.0711      0.0012  
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        -0.0005      0.0000 (0.05 %)     0.0368 (96.96 %)  0.0011 ( 2.99 %)  0.0380      0.0009
                                                            att        -0.0005      0.0000 (0.05 %)     0.0689 (98.28 %)  0.0012 ( 1.67 %)  0.0702      0.0012
                    ==========================================================================================================================================

        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate    var_tau             var_sigmapi       var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        5.0           1,000  2*10^(-11)  0.8891    1        ate         0.0046      0.0003 (0.67 %)     0.0457 (89.16 %)  0.0052 (10.17 %)  0.0513      0.0072
                                                            att         0.0044      0.0007 (0.76 %)     0.0815 (93.76 %)  0.0048 ( 5.48 %)  0.0869      0.0093
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.0057      0.0003 (0.42 %)     0.0524 (86.45 %)  0.0080 (13.13 %)  0.0606      0.0078
                                                            att         0.0033      0.0006 (0.53 %)     0.0979 (93.02 %)  0.0068 ( 6.45 %)  0.1052      0.0103
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         0.0036      0.0003 (0.62 %)     0.0482 (85.58 %)  0.0078 (13.80 %)  0.0564      0.0075
                                                            att         0.0077      0.0006 (0.51 %)     0.1133 (92.98 %)  0.0079 ( 6.51 %)  0.1218      0.0110
                    ==========================================================================================================================================
                      5,000  1*10^(-11)  0.5946    1        ate         0.0045      0.0003 (0.58 %)     0.0470 (92.43 %)  0.0036 ( 6.99 %)  0.0509      0.0032
                                                            att         0.0024      0.0006 (0.60 %)     0.0904 (94.87 %)  0.0043 ( 4.53 %)  0.0953      0.0044
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate        -0.0005      0.0002 (0.41 %)     0.0458 (94.04 %)  0.0027 ( 5.55 %)  0.0487      0.0031 
                                                            att        -0.0015      0.0004 (0.44 %)     0.0874 (97.23 %)  0.0021 ( 2.34 %)  0.0899      0.0042
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        -0.0021      0.0002 (0.49 %)     0.0465 (93.33 %)  0.0031 ( 6.19 %)  0.0498      0.0032
                                                            att        -0.0017      0.0005 (0.48 %)     0.0929 (96.35 %)  0.0031 ( 3.18 %)  0.0964      0.0044
                    ==========================================================================================================================================
                     50,000  7*10^(-12)  0.3344    1        ate        -0.0022      0.0001 (0.19 %)     0.0453 (99.09 %)  0.0003 ( 0.71 %)  0.0457      0.0010
                                                            att        -0.0028      0.0002 (0.19 %)     0.0889 (99.23 %)  0.0005 ( 0.58 %)  0.0896      0.0013
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.0004      0.0001 (0.21 %)     0.0452 (99.11 %)  0.0003 ( 0.68 %)  0.0457      0.0010
                                                            att        -0.0005      0.0002 (0.19 %)     0.0918 (99.32 %)  0.0005 ( 0.49 %)  0.0925      0.0014
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        -0.0000      0.0001 (0.24 %)     0.0448 (99.16 %)  0.0003 ( 0.60 %)  0.0452      0.0010
                                                            att         0.0001      0.0002 (0.28 %)     0.0884 (99.48 %)  0.0002 ( 0.28 %)  0.0889      0.0013
     Findings: 1. Point 1. in (1) DGP applies.
            For ate:
               2. Point 2. in (1) DGP applies, but instead of Point 2.1,
                    2.1 the speed of decrease is roughly the same for the two kappa_gamma.
               3. For the absolute value of var_tau, Point 2. applies (lower values for kappa_gamma=0.5), but the speed is faster for kappa_gamma=5.0.
               4. Point 4. in (1) DGP applies.
               5. Point 5. in (1) DGP applies.
               6. Point 6. in (1) DGP applies.
            For att:
               7. The (proportion of estimated var_tau to var_total) decreases to the true value zero steadily for kappa_gamma=5, but:
                    7.1 it reaches lower percentages for kappa_gamma=0.5.
               8. Point 7. applies to the absolute value of var_tau.
               9. Point 9. in (1) DGP applies.
               10. 
            For both:
               11. Point 11. in (1) DGP applies.  

    Conclusions regarding bandwidth (based on (1) DGP (=no truncation) and (2) DGP (=mild truncation)): 
                1. kappa_gamma=0.5 seems better for var_tau but worse for var_estpi, while it's the other way around for kappa_gamma=5.0.
                   This suggests the use of different bandwidths for var_tau and var_estpi, presumably because of the different scales (propensity score vs score) in derivative estimation.
                   var_estpi uses both scales (propensity score vs score), so when the same bandwidth is used, it needs to balance the scales.
                2. Repeating (2) DGP for kappa_gamma=50.0 and n=50,000, gives worse results for var_tau (ate:0.0005 (0.65 %), att:0.0009 (0.65 %)) and especially var_estpi (ate:0.03 (38 %, att:0.04 (31 %)), but not too bad for se (ate:0.0013, att:0.0016).
                3. Repeating (2) DGP for kappa_gamma=10.0 and n=50,000  gives better results than kappa_gamma=50.0 but worse than kappa_gamma=5.0.
    Conclusions regarding truncation:           
                1. Repeating (2) DGP for kappa_gamma=5.0, n=50,000 and kappa_a=10^(-1) gives worse results for var_estpi (ate:0.0060 (26.39 %), att:0.0077 (20.05 %)) compared to (2) DGP (kappa_a=10^(-11)), but not too bad for se (ate:0.0007, att:0.0009).
                2. Repeating (2) DGP for kappa_gamma=5.0, n=50,000 and kappa_a=10^(-30) gives the same results as (2) DGP (kappa_a=10^(-10)).
   
   - (3) DGP:
       Same as (2) DGP, except mu1=1.
       Covariates: x_1, x_2 ~ N(0,1) i.i.d.;
       Propensity score: known (= nonestimated) logit, values close to 0 and 1, theta0_0=0.5, theta0_1=0.2, theta0_2=intercept=0.04;
       Potential outcomes:
       mu0 = 0
       mu1 = 1 
       y_0i = mu0 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
       y_1i = mu1 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        ==> ate=att=1; mu0(pi(x))-mu1(pi(x))=1 for all x, so var_tau and var_estpi should both be zero for both ate and att.
     Settings:
        Caliper = maxmindelta2;
        Truncation:  kappa_a = 10^(-10);
                     npower_a = 1/4.000000002
        Bandwidth:   npower_gamma = 1/4.000000001;
     Results from a few repetitons:
        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate    var_tau             var_sigmapi       var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        5.0          50,000  7*10^(-12)  0.3344    1        ate         1.0004      0.0207 (31.5 %)     0.0448 (68.10 %)  0.0003 (0.40 %)   0.0658      0.0011  
                                                            att         1.0004      0.0395 (31.53 %)    0.0855 (68.26 %)  0.0003 (0.22 %)   0.1253      0.0016
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         1.0011      0.0194 (30.57 %)    0.0437 (68.69 %)  0.0003 (0.55 %)   0.0634      0.0011
                                                            att         1.0018      0.0362 (30.15 %)    0.0837 (69.63 %)  0.0003 (0.23 %)   0.1202      0.0016
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         1.0005      0.0195 (30.34 %)    0.0445 (69.13 %)  0.0003 (0.53 %)   0.0643      0.0011
                                                            att         1.0011      0.0368 (29.65 %)    0.0870 (70.01 %)  0.0004 (0.34 %)   0.1243      0.0016
    Comparison with (2) DGP: 1. We expect var_tau, var_sigmapi and var_estpi to be the same as in (2) DGP  (mu0(pi(x))-mu1(pi(x))=0),
                                because there's only a constant shift mu0(pi(x))-mu1(pi(x))=1. This is met for var_sigmapi and var_estpi.
                             2. Indeed, all things same, except the absolute value of var_tau is much larger here. 
                                => Maybe kappa_a=10^(-10) is too small for when the treatment effect is not zero?
                                Nonetheless, the se's are still close to (2) DGP!
                    
   - (4) DGP:
       Same as (3) DGP, except kappa_a.
       Covariates: x_1, x_2 ~ N(0,1) i.i.d.;
       Propensity score: known (= nonestimated) logit, values close to 0 and 1, theta0_0=0.5, theta0_1=0.2, theta0_2=intercept=0.04;
       Potential outcomes:
       mu0 = 0
       mu1 = 1 
       y_0i = mu0 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
       y_1i = mu1 + 0.001 * matrix_get(x, i, 0) + 0.005 * matrix_get(x, i, 1) + 0.009 * gsl_ran_ugaussian(rng);
        ==> ate=att=1; mu0(pi(x))-mu1(pi(x))=1 for all x, so var_tau and var_estpi should both be zero for both ate and att.
     Settings:
        Caliper = maxmindelta2;
        Truncation:  kappa_a = 10^(-2);
                     npower_a = 1/4.000000002
        Bandwidth:   npower_gamma = 1/4.000000001;
     Results from a few repetitons:
        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate    var_tau             var_sigmapi       var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        5.0          50,000    0.0007    0.3344    1        ate         1.0004      0.0114 (20.77 %)    0.0393 (71.28 %)  0.0044 (7.95 %)   0.0551      0.0011  
                                                            att         1.0010      0.0197 (19.45 %)    0.0759 (75.09 %)  0.0055 (5.46 %)   0.1011      0.0014
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         1.0000      0.0128 (22.85 %)    0.0394 (70.45 %)  0.0037 (6.7 %)    0.0559      0.0011
                                                            att         1.0006      0.0222 (21.80 %)    0.0754 (74.00 %)  0.0043 (4.21 %)   0.1019      0.0014
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         0.9999      0.0124 (22.33 %)    0.0392 (70.89 %)  0.0038 (6.78 %)   0.0554      0.0011
                                                            att         0.9998      0.0254 (23.11 %)    0.0793 (72.23 %)  0.0051 (4.66 %)   0.1098      0.0015
     Comparison with (3) DGP:
        1. For kappa_a=10^(-2) var_tau become smaller compared to kappa_a=10^(-10) in (3) DGP.
        2. For kappa_a=10^(-2) var_sigmapi is roughly as for kappa_a=10^(-10) in (3) DGP
        3. For kappa_a=10^(-2) var_estpi became larger compared to kappa_a=10^(-10) in (3) DGP.
        4. For kappa_a=10^(-2) var_total decreased slightly (se's are roughly same).
     
     Now with kappa_a = 10^0:
        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate    var_tau             var_sigmapi       var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        5.0          50,000    0.0669    0.3344    1        ate         1.0007     -0.4652 (!!!!!!!)    0.0030 (!!!!!!!)  0.0005 (!!!!!!)   !!!!!!!     !!!!!!  
                                                            att         1.0016     -0.9239 (!!!!!!!)    0.0043 (!!!!!!!)  0.0005 (!!!!!!)   !!!!!!!     !!!!!! 
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.9991     -0.4636 (!!!!!!!)    0.0030 (!!!!!!!)  0.0005 (!!!!!!)   !!!!!!!     !!!!!!
                                                            att         0.9991     -0.9318 (!!!!!!!)    0.0043 (!!!!!!!)  0.0005 (!!!!!!)   !!!!!!!     !!!!!!
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         0.9997     -0.4620 (!!!!!!!)    0.0031 (!!!!!!!)  0.0004 (!!!!!!)   !!!!!!!     !!!!!!
                                                            att         0.9990     -0.9226 (!!!!!!!)    0.0045 (!!!!!!!)  0.0004 (!!!!!!)   !!!!!!!     !!!!!!
     
     Comparison with (3) DGP:
        1. Truncation sequence too large => truncation too strong => too many observations zero out in the first term of estimated var_tau => 
           negative estimated var_tau because att_hat is accurate.
       
   - Conclusions based on (1)--(4) DGP:
        1. Bandwidth: suggests the use of different bandwidths for var_tau and var_estpi, presumably because of the different scales (propensity score vs score) in derivative estimation.
                      var_estpi uses both scales (propensity score vs score), so when the same bandwidth is used, it needs to balance the scales.
        2. Truncation: use small rather than large truncation sequences to avoid negative var_tau, espicially for att! This may be too conservative.
        3. Shares of variance components are sensitive to bandwidth but var_total is not so much.
        4. Too strong truncation can ruin var_tau and thus var_total for nonzero treatment effects.
   - (5) DGP:
       Nonconstant mu1(pi(x))-mu0(pi(x)) in x.
       Covariates: x_1, x_2 ~ N(0,1) i.i.d.;
       Propensity score: known (= nonestimated) logit, values close to 0 and 1, theta0_0=0.5, theta0_1=0.2, theta0_2=intercept=0.04;
       Potential outcomes:
       mu0 = 0
       mu1 = 0 
       y_0i = mu0 + 0.8 * matrix_get(x, i, 0) + 1.5 * matrix_get(x, i, 1) + 0.9 * gsl_ran_ugaussian(rng);
       y_1i = mu1 + 3.1 * sin(matrix_get(x, i, 0)) - 1.8 * matrix_get(x, i, 1) + 1.1 * gsl_ran_ugaussian(rng);
        ==> ate=0; att=?
     Settings:
        Caliper = maxmindelta2;
        Truncation:  kappa_a = 10^(-10);
                     npower_a = 1/4.000000002
        Bandwidth:   npower_gamma = 1/4.000000001;
     Results from a few repetitons:     
        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate    var_tau             var_sigmapi         var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        0.5          50,000  7*10^(-12)  0.0334    1        ate         -0.2555     77.2828 (1.94 %)    3797.2900 (97.79 %)  10.3985 (0.27 %) 3882.9744 0.2787
                                                            att         -7.5104     43.0499 (0.89 %)    4667.7443 (96.00 %) 151.2534 (3.11 %) 4862.0476 0.3118   
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate          0.1812     79.6076 (1.95 %)    3986.2228 (97.47 %)  23.8240 (0.58 %) 4089.6544 0.2860  
                                                            att         -7.0665     47.7116 (0.94 %)    4864.6109 (95.84 %) 163.4902 (3.22 %) 5075.8127 0.3186        
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate         -0.3842     74.7397 (1.93 %)    3786.6167 (97.93 %)   5.4315 (0.14 %) 3866.7879 0.2781 
                                                            att         -7.7428     38.0899 (0.76 %)    4860.3724 (97.29 %)  97.1277 (1.94 %) 4995.5900 0.3161
   
        kappa_gamma     n       a_n     gamma_n   rep.id    estimand    estimate    var_tau             var_sigmapi         var_estpi        var_total     se
        ------------------------------------------------------------------------------------------------------------------------------------------------------    
        5.0          50,000  7*10^(-12)  0.3344    1        ate        -0.0613      31.4747 (0.70 %)   4477.9339 (99.05 %)  11.3019 (0.25 %) 4520.7105  0.3007
                                                            att        -7.0785     -51.9562 (!!!!!!)   5709.0297 (!!!!!!!) 958.3225 (!!!!!!) !!!!!!!!!  !!!!!! 
                                                   -----------------------------------------------------------------------------------------------------------
                                                   2        ate         0.3141      33.2990 (0.74 %)   4459.3453 (99.01 %)  11.5028 (0.26 %) 4504.1471  0.3001 
                                                            att        -7.2144     -62.0993 (!!!!!!)   5779.2317 (!!!!!!!) 968.6504 (!!!!!!) !!!!!!!!!  !!!!!! 
                                                   -----------------------------------------------------------------------------------------------------------
                                                   3        ate        -0.1381      31.9069 (0.71 %)   4456.4076 (98.96 %)  14.9623 (0.33 %) 4556.4076  0.3001    
                                                            att        -7.5714     -68.2416 (!!!!!!)   5804.0965 (!!!!!!!) 964.3565 (!!!!!!) !!!!!!!!!  !!!!!!
     Conclusions: large kappa_gamma can lead to negative var_tau as larger bandwidth induces larger bias. Use small kappa_gamma! 
   */                                                            

    // Cleanup
    vector_short_free(d); vector_free(y);
    matrix_free(x);
    free(cm_model);
    free(cm_results);
    return 0;
}