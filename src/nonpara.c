////////////////////////////////////////////////////////////////////////////
//
// * FILE: nonpara.c
// * DESCRIPTION:
//    Nonparametric tools.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 08/nov/2022
//
//////////////////////////////////////////////////////////////////////////

#include "nonpara.h"


//////////// Kernels

double nonpara_kernel_gauss(double u){
    return (1.0 / (sqrt(M_PI) * M_SQRT2)) * exp(-pow(u, 2));
}

// Derivative of nonpara_kernel_gauss
double nonpara_kernelderiv_gauss(double u){
    return nonpara_kernel_gauss(u) * (-u);
}

////////// Kernel estimators

double nonpara_qk_eval(double x_eval, double *x, unsigned long int n_x, double *y, unsigned long int n_y, double k, double (*kernel)(double), double bandwidth){
    assert(n_x == n_y);
    unsigned long int n = n_x;
    double eval = 0.0;
    for (unsigned long int i=0; i<n; i++){
        eval += pow(y[i], k) * kernel((x_eval - x[i]) / bandwidth) / (n * bandwidth);
    }
    return eval;
}

double nonpara_h_eval(double x_eval, double *x, unsigned long int n, double (*kernel)(double), double bandwidth){
    double eval = 0.0;
    for (unsigned long int i=0; i<n; i++){
        eval += kernel((x_eval - x[i]) / bandwidth) / (n * bandwidth);
    }
    return eval;
}

// Assesses whether `x` and `y` are closer than tolarance threshold.
int _test_nonpara_dcompare(double x, double y){
    return fabs(x - y) <= pow(10, -NONPARA_TEST_DOUBLE_TOLERANCE_NPOWER);
}

////////// Kernel estimators with filters





/*  Evaluates h_b := (1 / (bandwidth * n_adjusted_b)) * sum_{i in [n]: filter(i)==b} K((x_eval-x_i)/ bandwidth) for b=0,1, where
 - n is the length of `x`;
 - x_i are elements  i=0,...,n-1 of `x`.
 - K is the kernel `kernel`;
 - filter(i) is 0 if `filter_func` evaluates to false, 1 otherwise;
 - `filter_func` takes and observation index and additional argument `filter_args` passed as void pointer.
 - n_adjusted_b is the number of observations i such that filter(i)=b if adjust_nobs is true, n_adjusted_b = n otherwise;
 
 Returns 2-long pointer `eval_value` with eval_value[0] = h_0 and eval_value[1] = h_1.
*/
double *nonpara_h_eval_filter(double x_eval, double *x, unsigned long int n, unsigned long int stride, double (*kernel)(double), double bandwidth, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *)){
    unsigned long int n_true = 0;
    unsigned long int n_false = 0;
    double eval_value_false, eval_value_true;
    eval_value_false = eval_value_true = 0;
    for (unsigned long int i=0; i<n; i++){
        if (filter_fn(i, filter_args)){  // filter true
            eval_value_true += kernel((x_eval - x[stride * i]) / bandwidth) / (n * bandwidth);
            n_true++;
        } else {    // filter false
            eval_value_false += kernel((x_eval - x[stride * i]) / bandwidth) / (n * bandwidth);
            n_false++;
        }
    }
    double adjustment_factor_nobs_true, adjustment_factor_nobs_false;
    if (adjust_nobs){
        adjustment_factor_nobs_true = (n + 0.0) / n_true;
        adjustment_factor_nobs_false = (n + 0.0) / n_false; 
    } else{
        adjustment_factor_nobs_false = adjustment_factor_nobs_true = 1;
    }
    double *eval_value = malloc(2*sizeof(double)); // entry 0 = filter false; entry 1 = filter true
    eval_value[0] = eval_value_false * adjustment_factor_nobs_false;
    eval_value[1] = eval_value_true * adjustment_factor_nobs_true;
    return eval_value;
}

// Auxiliary function for test_h_deriv_eval_filter.
int _test_nonpara_h_eval_filter_filter_func(unsigned long int i, void *filter_arg){
    int *d = (int *) filter_arg;
    return d[i];
}

void test_nonpara_h_eval_filter(void){
    // test 1: only one group d=1
    int n = 4;
    int stride = 1;
    double bandwidth = pow(n, -1.0/5);
    double (*kernel)(double);
    int (*filter_func)(unsigned long int, void *);
    kernel = nonpara_kernel_gauss;
    filter_func = _test_nonpara_h_eval_filter_filter_func;
    int adjust_obs = 0;
    int *d = malloc(n*sizeof(int));  // filter variable with two categories
    double *x = malloc(n*sizeof(double));
    d[0] = 1;   x[0] = 1.0;
    d[1] = 1;   x[1] = -2.9;
    d[2] = 1;   x[2] = 4.0;
    d[3] = 1;   x[3] = 0.0;
    double x_eval = 0.5;
    double *h_d = malloc(2*sizeof(double));
    h_d[0] = 0.0; // no observation with d=0
    h_d[1] = (1.0/(n * bandwidth)) * (kernel((x_eval- x[0])/bandwidth) + kernel((x_eval- x[1])/bandwidth) + kernel((x_eval- x[2])/bandwidth) + kernel((x_eval- x[3])/bandwidth));
    double *result = nonpara_h_eval_filter(x_eval, x, n, stride, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("h_d[0] = %f, result[0] = %f\n", h_d[0], result[0]);
    // printf("h_d[1] = %f, result[1] = %f\n", h_d[1], result[1]);
    assert( _test_nonpara_dcompare(h_d[0], result[0]) && _test_nonpara_dcompare(h_d[1],result[1]));
    
    // test 2: only one group d=0
    d[0] = 0;   x[0] = 1.0;
    d[1] = 0;   x[1] = -2.9;
    d[2] = 0;   x[2] = 4.0;
    d[3] = 0;   x[3] = 0.0;
    x_eval = 0.5;
    h_d[0] = (1.0/(n * bandwidth)) * (kernel((x_eval- x[0])/bandwidth) + kernel((x_eval- x[1])/bandwidth) + kernel((x_eval- x[2])/bandwidth) + kernel((x_eval- x[3])/bandwidth));
    h_d[1] = 0.0; // no observation with d!=0
    result = nonpara_h_eval_filter(x_eval, x, n, stride, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("h_d[0] = %f, result[0] = %f\n", h_d[0], result[0]);
    // printf("h_d[1] = %f, result[1] = %f\n", h_d[1], result[1]);
    assert( _test_nonpara_dcompare(h_d[0], result[0]) && _test_nonpara_dcompare(h_d[1],result[1]));

    // test 3: 2-2 observations in two groups
    d[0] = 0;   x[0] = 1.0;
    d[1] = 0;   x[1] = -2.9;
    d[2] = 1;   x[2] = 4.0;
    d[3] = 1;   x[3] = 0.0;
    x_eval = 0.7;
    h_d[0] = (1.0/(n * bandwidth)) * (kernel((x_eval- x[0])/bandwidth) + kernel((x_eval- x[1])/bandwidth));
    h_d[1] = (1.0/(n * bandwidth)) * (kernel((x_eval- x[2])/bandwidth) + kernel((x_eval- x[3])/bandwidth));
    result = nonpara_h_eval_filter(x_eval, x, n, stride, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("h_d[0] = %f, result[0] = %f\n", h_d[0], result[0]);
    // printf("h_d[1] = %f, result[1] = %f\n", h_d[1], result[1]);
    assert( _test_nonpara_dcompare(h_d[0], result[0]) && _test_nonpara_dcompare(h_d[1],result[1]));

    // test 4: 2-2 observations in two groups, adjust number of observations
    adjust_obs = 1;
    d[0] = 0;   x[0] = 1.0;
    d[1] = 1;   x[1] = -2.9;
    d[2] = 0;   x[2] = 4.0;
    d[3] = 1;   x[3] = 0.0;
    x_eval = 0.7;
    h_d[0] = (1.0/((n / 2.0) * bandwidth)) * (kernel((x_eval- x[0])/bandwidth) + kernel((x_eval- x[2])/bandwidth));
    h_d[1] = (1.0/((n / 2.0) * bandwidth)) * (kernel((x_eval- x[1])/bandwidth) + kernel((x_eval- x[3])/bandwidth));
    result = nonpara_h_eval_filter(x_eval, x, n, 1, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("h_d[0] = %f, result[0] = %f\n", h_d[0], result[0]);
    // printf("h_d[1] = %f, result[1] = %f\n", h_d[1], result[1]);
    assert( _test_nonpara_dcompare(h_d[0], result[0]) && _test_nonpara_dcompare(h_d[1],result[1]));

    // test 4: 2-2 observations in two groups, adjust number of obs, stride_x = 2
    adjust_obs = 1;
    int stride_x = 2;
    double blank_value = 1000000;
    double *x_s = malloc(stride_x * n * sizeof(double));
    d[0] = 0;   x_s[0] = 1.0;
                x_s[1] = blank_value;
    d[1] = 1;   x_s[2] = -2.9;
                x_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;  
                x_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 0.0;
                x_s[7] = blank_value;
    x_eval = 0.5;
    h_d[0] = (1.0/((n / 2.0) * bandwidth)) * (kernel((x_eval - x_s[0])/bandwidth) +
                                      kernel((x_eval - x_s[4])/bandwidth));
    h_d[1] = (1.0/((n / 2.0) * bandwidth)) * (kernel((x_eval - x_s[2])/bandwidth) + 
                                      kernel((x_eval - x_s[6])/bandwidth));
    result = nonpara_h_eval_filter(x_eval, x_s, n, stride_x, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("h_d[0] = %f, result[0] = %f\n", h_d[0], result[0]);
    // printf("h_d[1] = %f, result[1] = %f\n", h_d[1], result[1]);
    assert( _test_nonpara_dcompare(h_d[0], result[0]) && _test_nonpara_dcompare(h_d[1],result[1]));
    
    // cleanup
    free(x);
    free(x_s);
    free(d);
    free(result);
    free(h_d);
    printf("Tests for `nonpara_h_eval_filter` completed: no issues found.\n");
}


/*  Evaluates q_b := (1 / (bandwidth * n_adjusted_b)) * sum_{i in [n]: filter(i)==b} (y_i^k) * K((x_eval-x_i)/ bandwidth) for b=0,1, where
 - n is the length of `x`;
 - x_i, y_i are elements  i=0,...,n-1 of `x` and `y`.
 - K is the kernel `kernel`;
 - filter(i) is 0 if `filter_func` evaluates to false, 1 otherwise;
 - `filter_func` takes and observation index and additional argument `filter_args` passed as void pointer.
 - n_adjusted_b is the number of observations i such that filter(i)=b if adjust_nobs is true, n_adjusted_b = n otherwise;
 
 Returns 2-long pointer `eval_value` with eval_value[0] = q_0 and eval_value[1] = q_1.
*/
double *nonpara_qk_eval_filter(double x_eval, double *x, unsigned long int n_x, unsigned long int stride_x, double *y, unsigned long int n_y, unsigned long int stride_y, double k, double (*kernel)(double), double bandwidth, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *)){
    assert(n_x == n_y);
    unsigned long int n = n_x;
    unsigned long int n_true = 0;
    unsigned long int n_false = 0;
    double eval_value_false, eval_value_true;
    eval_value_false = eval_value_true = 0;
    for (unsigned long int i=0; i<n; i++){
        if (filter_fn(i, filter_args)){  // filter true
            eval_value_true += pow(y[stride_y * i], k) * kernel((x_eval - x[stride_x * i]) / bandwidth) / (n * bandwidth);
            n_true++;
        } else {    // filter false
            eval_value_false += pow(y[stride_y * i], k) * kernel((x_eval - x[stride_x * i]) / bandwidth) / (n * bandwidth);
            n_false++;
        }
    }
    double adjustment_factor_nobs_true, adjustment_factor_nobs_false;
    if (adjust_nobs){
        adjustment_factor_nobs_true = (n + 0.0) / n_true;
        adjustment_factor_nobs_false = (n + 0.0) / n_false; 
    } else{
        adjustment_factor_nobs_false = adjustment_factor_nobs_true = 1;
    }
    double *eval_value = malloc(2*sizeof(double)); // entry 0 = filter false; entry 1 = filter true
    eval_value[0] = eval_value_false * adjustment_factor_nobs_false;
    eval_value[1] = eval_value_true * adjustment_factor_nobs_true;
    return eval_value;
}



// Auxiliary function for test_nonpara_qk_eval_filter.
int _test_nonpara_qk_eval_filter_filter_func(unsigned long int i, void *filter_arg){
    int *d = (int *) filter_arg;
    return d[i];
}

void test_nonpara_qk_eval_filter(void){
    // test 1: only one group d=1
    int n = 4;
    int stride_x, stride_y;    
    stride_x = stride_y = 1;
    double bandwidth = pow(n, -1.0/5);
    double (*kernel)(double);
    int (*filter_func)(unsigned long int, void *);
    kernel = nonpara_kernel_gauss;
    filter_func = _test_nonpara_qk_eval_filter_filter_func;
    int adjust_obs = 0;
    double power = 1;
    int *d = malloc(n*sizeof(int));  // filter variable with two categories
    double *x = malloc(n*sizeof(double));
    double *y = malloc(n*sizeof(double));
    d[0] = 1;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 1;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 1;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;
    double x_eval = 0.5;
    double *q_d = malloc(2*sizeof(double));
    q_d[0] = 0.0; // no observation with d=0
    q_d[1] = (1.0/(n * bandwidth)) * (pow(y[0], power) * kernel((x_eval- x[0])/bandwidth) +
                                      pow(y[1], power) * kernel((x_eval- x[1])/bandwidth) + 
                                      pow(y[2], power) * kernel((x_eval- x[2])/bandwidth) + 
                                      pow(y[3], power) * kernel((x_eval- x[3])/bandwidth));
    double *result = nonpara_qk_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("q_d[0] = %f, result[0] = %f\n", q_d[0], result[0]);
    // printf("q_d[1] = %f, result[1] = %f\n", q_d[1], result[1]);
    assert( _test_nonpara_dcompare(q_d[0], result[0]) && _test_nonpara_dcompare(q_d[1],result[1]));
    
    // test 2: only one group d=0
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 0;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 0;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 0;   x[3] = 0.0;     y[3] = 50.4;
    x_eval = 0.5;
    q_d[0] = (1.0/(n * bandwidth)) * (pow(y[0], power) * kernel((x_eval- x[0])/bandwidth) +
                                      pow(y[1], power) * kernel((x_eval- x[1])/bandwidth) + 
                                      pow(y[2], power) * kernel((x_eval- x[2])/bandwidth) + 
                                      pow(y[3], power) * kernel((x_eval- x[3])/bandwidth));
    q_d[1] = 0.0; // no observation with d=0
    result = nonpara_qk_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("q_d[0] = %f, result[0] = %f\n", q_d[0], result[0]);
    // printf("q_d[1] = %f, result[1] = %f\n", q_d[1], result[1]);
    assert( _test_nonpara_dcompare(q_d[0], result[0]) && _test_nonpara_dcompare(q_d[1],result[1]));

    // test 3: 2-2 observations in two groups
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 0;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 1;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;
    x_eval = 0.5;
    q_d[0] = (1.0/(n * bandwidth)) * (pow(y[0], power) * kernel((x_eval- x[0])/bandwidth) +
                                      pow(y[1], power) * kernel((x_eval- x[1])/bandwidth));
    q_d[1] = (1.0/(n * bandwidth)) * (pow(y[2], power) * kernel((x_eval- x[2])/bandwidth) + 
                                      pow(y[3], power) * kernel((x_eval- x[3])/bandwidth));
    result = nonpara_qk_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("q_d[0] = %f, result[0] = %f\n", q_d[0], result[0]);
    // printf("q_d[1] = %f, result[1] = %f\n", q_d[1], result[1]);
    assert( _test_nonpara_dcompare(q_d[0], result[0]) && _test_nonpara_dcompare(q_d[1],result[1]));
    
    // test 3: 2-2 observations in two groups, adjust number of obs
    adjust_obs = 1;
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 1;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 0;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;
    x_eval = 0.5;
    q_d[0] = (1.0/((n / 2.0) * bandwidth)) * (pow(y[0], power) * kernel((x_eval- x[0])/bandwidth) +
                                      pow(y[2], power) * kernel((x_eval- x[2])/bandwidth));
    q_d[1] = (1.0/((n / 2.0) * bandwidth)) * (pow(y[1], power) * kernel((x_eval- x[1])/bandwidth) + 
                                      pow(y[3], power) * kernel((x_eval- x[3])/bandwidth));
    result = nonpara_qk_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("q_d[0] = %f, result[0] = %f\n", q_d[0], result[0]);
    // printf("q_d[1] = %f, result[1] = %f\n", q_d[1], result[1]);
    assert( _test_nonpara_dcompare(q_d[0], result[0]) && _test_nonpara_dcompare(q_d[1],result[1]));

    // test 4: 2-2 observations in two groups, adjust number of obs, stride_x = 2
    adjust_obs = 1;
    stride_x = 2;
    double blank_value = 1000000;
    double *x_s = malloc(stride_x * n * sizeof(double));
    d[0] = 0;   x_s[0] = 1.0;     y[0] = 10.1;
                x_s[1] = blank_value;
    d[1] = 1;   x_s[2] = -2.9;    y[1] = -30.5;
                x_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;     y[2] = 100.45;
                x_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 0.0;     y[3] = 50.4;
                x_s[7] = blank_value;
    x_eval = 0.5;
    q_d[0] = (1.0/((n / 2.0) * bandwidth)) * (pow(y[0], power) * kernel((x_eval - x_s[0])/bandwidth) +
                                      pow(y[2], power) * kernel((x_eval - x_s[4])/bandwidth));
    q_d[1] = (1.0/((n / 2.0) * bandwidth)) * (pow(y[1], power) * kernel((x_eval - x_s[2])/bandwidth) + 
                                      pow(y[3], power) * kernel((x_eval - x_s[6])/bandwidth));
    result = nonpara_qk_eval_filter(x_eval, x_s, n, stride_x, y, n, stride_y, power, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("q_d[0] = %f, result[0] = %f\n", q_d[0], result[0]);
    // printf("q_d[1] = %f, result[1] = %f\n", q_d[1], result[1]);
    assert( _test_nonpara_dcompare(q_d[0], result[0]) && _test_nonpara_dcompare(q_d[1], result[1]));

    // test 5: 2-2 observations in two groups, adjust number of obs, stride_y = 2
    stride_x = 1;
    stride_y = 2;
    double *y_s = malloc(stride_y * n * sizeof(double));
    d[0] = 0;   x[0] = 1.0;     y_s[0] = 10.1;
                                y_s[1] = blank_value;
    d[1] = 1;   x[1] = -2.9;    y_s[2] = -30.5;
                                y_s[3] = blank_value;
    d[2] = 0;   x[2] = 4.0;     y_s[4] = 100.45;
                                y_s[5] = blank_value;
    d[3] = 1;   x[3] = 0.0;     y_s[6] = 50.4;
                                y_s[7] = blank_value;
    x_eval = 0.5;
    q_d[0] = (1.0/((n / 2.0) * bandwidth)) * (pow(y_s[0], power) * kernel((x_eval - x[0])/bandwidth) +
                                      pow(y_s[4], power) * kernel((x_eval- x[2])/bandwidth));
    q_d[1] = (1.0/((n / 2.0) * bandwidth)) * (pow(y_s[2], power) * kernel((x_eval - x[1])/bandwidth) + 
                                      pow(y_s[6], power) * kernel((x_eval- x[3])/bandwidth));
    result = nonpara_qk_eval_filter(x_eval, x, n, stride_x, y_s, n, stride_y, power, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("q_d[0] = %f, result[0] = %f\n", q_d[0], result[0]);
    // printf("q_d[1] = %f, result[1] = %f\n", q_d[1], result[1]);
    assert( _test_nonpara_dcompare(q_d[0], result[0]) && _test_nonpara_dcompare(q_d[1], result[1]));

    // test 6: 2-2 observations in two groups, adjust number of obs, stride_x = 2, stride_y = 2
    stride_x = 2;
    stride_y = 2;
    d[0] = 0;   x_s[0] = 1.0;          y_s[0] = 10.1;
                x_s[1] = blank_value;  y_s[1] = blank_value;
    d[1] = 1;   x_s[2] = -2.9;         y_s[2] = -30.5;
                x_s[3] = blank_value;  y_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;          y_s[4] = 100.45;
                x_s[5] = blank_value;  y_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 0.0;          y_s[6] = 50.4;
                x_s[7] = blank_value;  y_s[7] = blank_value;
    x_eval = 0.5;
    q_d[0] = (1.0/((n / 2.0) * bandwidth)) * (pow(y_s[0], power) * kernel((x_eval - x_s[0])/bandwidth) +
                                      pow(y_s[4], power) * kernel((x_eval - x_s[4])/bandwidth));
    q_d[1] = (1.0/((n / 2.0) * bandwidth)) * (pow(y_s[2], power) * kernel((x_eval - x_s[2])/bandwidth) + 
                                      pow(y_s[6], power) * kernel((x_eval - x_s[6])/bandwidth));
    result = nonpara_qk_eval_filter(x_eval, x_s, n, stride_x, y_s, n, stride_y, power, kernel, bandwidth, adjust_obs, (void *) d, filter_func);
    // printf("q_d[0] = %f, result[0] = %f\n", q_d[0], result[0]);
    // printf("q_d[1] = %f, result[1] = %f\n", q_d[1], result[1]);
    assert( _test_nonpara_dcompare(q_d[0], result[0]) && _test_nonpara_dcompare(q_d[1], result[1]));
    // cleanup
    free(x);
    free(x_s);
    free(d);
    free(y);
    free(y_s);
    free(result);
    free(q_d);
    printf("Tests for `nonpara_qk_eval_filter` completed: no issues found.\n");
}



////////// Kernel derivative estimators with filters


/*  Evaluates h_deriv_b := (- scale / (bandwidth^2 * n_adjusted_b)) * sum_{i in [n]: filter(i)==b} K'((x_i-x_eval)/ bandwidth) for b=0,1, where
 - `scale` is a real scalar;
 - n is the length of `x`;
 - x_i are elements  i=0,...,n-1 of `x`.
 - K' is the derivative `kernel_deriv` of a kernel;
 - filter(i) is 0 if `filter_func` evaluates to false, 1 otherwise;
 - `filter_func` takes and observation index and additional argument `filter_args` passed as void pointer.
 - n_adjusted_b is the number of observations i such that filter(i)=b if adjust_nobs is true, n_adjusted_b = n otherwise;
 
 Returns 2-long pointer `eval_value` with eval_value[0] = h_deriv_0 and eval_value[1] = h_deriv_1.
*/
double *nonpara_h_deriv_eval_filter(double x_eval, double *x, unsigned long int n, unsigned long int stride, double (*kernel_deriv)(double), double bandwidth, double scale, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *)){
    unsigned long int n_true = 0;
    unsigned long int n_false = 0;
    double eval_value_false, eval_value_true;
    eval_value_false = eval_value_true = 0;
    for (unsigned long int i=0; i<n; i++){
        if (filter_fn(i, filter_args)){  // filter true
            eval_value_true += (-scale) * kernel_deriv((x[stride * i] - x_eval) / bandwidth) / (n * pow(bandwidth, 2));
            n_true++;
        } else {    // filter false
            eval_value_false += (-scale) * kernel_deriv((x[stride * i] - x_eval) / bandwidth) / (n * pow(bandwidth, 2));
            n_false++;
        }
    }
    double adjustment_factor_nobs_true, adjustment_factor_nobs_false;
    if (adjust_nobs){
        adjustment_factor_nobs_true = (n + 0.0) / n_true;
        adjustment_factor_nobs_false = (n + 0.0) / n_false; 
    } else{
        adjustment_factor_nobs_false = adjustment_factor_nobs_true = 1;
    }
    double *eval_value = malloc(2*sizeof(double)); // entry 0 = filter false; entry 1 = filter true
    eval_value[0] = eval_value_false * adjustment_factor_nobs_false;
    eval_value[1] = eval_value_true * adjustment_factor_nobs_true;
    return eval_value;
}

// Auxiliary function for test_h_deriv_eval_filter.
int _test_nonpara_h_deriv_eval_filter_filter_func(unsigned long int i, void *filter_arg){
    int *d = (int *) filter_arg;
    return d[i];
}

void test_nonpara_h_deriv_eval_filter(void){
    // test 1: only one group d=1, stride = 1
    int n = 4;
    int stride = 1;
    double scale = -3.9;
    double bandwidth = pow(n, -1.0/5);
    double (*kernel_deriv)(double);
    int (*filter_func)(unsigned long int, void *);
    kernel_deriv = nonpara_kernelderiv_gauss;
    filter_func = _test_nonpara_h_deriv_eval_filter_filter_func;
    int adjust_obs = 0;
    int *d = malloc(n*sizeof(int));  // filter variable with two categories
    double *x = malloc(n*sizeof(double));
    d[0] = 1;   x[0] = 1.0;
    d[1] = 1;   x[1] = -2.9;
    d[2] = 1;   x[2] = 4.0;
    d[3] = 1;   x[3] = 1.2;
    double x_eval = 0.5;
    double *h_d_deriv = malloc(2*sizeof(double));
    h_d_deriv[0] = 0.0; // no observation with d=0
    h_d_deriv[1] = (-scale/(n * pow(bandwidth, 2))) * (kernel_deriv((x[0] - x_eval)/bandwidth) + kernel_deriv((x[1] - x_eval)/bandwidth) + kernel_deriv((x[2] - x_eval)/bandwidth) + kernel_deriv((x[3] - x_eval)/bandwidth));
    double *result = nonpara_h_deriv_eval_filter(x_eval, x, n, stride, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("h_d_deriv[0] = %f, result[0] = %f\n", h_d_deriv[0], result[0]);
    // printf("h_d_deriv[1] = %.10f, result[1] = %.10f\n", h_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(h_d_deriv[0], result[0]) && _test_nonpara_dcompare(h_d_deriv[1],result[1]));

    // test 2: only one group d=0, stride = 1
    stride = 1;
    scale = -3.9;
    adjust_obs = 0;
    d[0] = 0;   x[0] = 1.0;
    d[1] = 0;   x[1] = -2.9;
    d[2] = 0;   x[2] = 4.0;
    d[3] = 0;   x[3] = 1.2;
    x_eval = 0.5;
    h_d_deriv[0] = (-scale/(n * pow(bandwidth, 2))) * (kernel_deriv((x[0] - x_eval)/bandwidth) + kernel_deriv((x[1] - x_eval)/bandwidth) + kernel_deriv((x[2] - x_eval)/bandwidth) + kernel_deriv((x[3] - x_eval)/bandwidth));
    h_d_deriv[1] = 0.0; // no observation with d=1
    result = nonpara_h_deriv_eval_filter(x_eval, x, n, stride, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("h_d_deriv[0] = %.10f, result[0] = %.10f\n", h_d_deriv[0], result[0]);
    // printf("h_d_deriv[1] = %.10f, result[1] = %.10f\n", h_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(h_d_deriv[0], result[0]) && _test_nonpara_dcompare(h_d_deriv[1],result[1]));
    
    // test 3: 2-2 observations in two groups
    d[0] = 0;   x[0] = 1.0;
    d[1] = 0;   x[1] = -2.9;
    d[2] = 1;   x[2] = 4.0;
    d[3] = 1;   x[3] = 1.2;
    x_eval = 0.5;
    h_d_deriv[0] = (-scale/(n * pow(bandwidth, 2))) * (kernel_deriv((x[0] - x_eval)/bandwidth) + kernel_deriv((x[1] - x_eval)/bandwidth));
    h_d_deriv[1] = (-scale/(n * pow(bandwidth, 2))) * (kernel_deriv((x[2] - x_eval)/bandwidth) + kernel_deriv((x[3] - x_eval)/bandwidth));; // no observation with d=1
    result = nonpara_h_deriv_eval_filter(x_eval, x, n, stride, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("h_d_deriv[0] = %.10f, result[0] = %.10f\n", h_d_deriv[0], result[0]);
    // printf("h_d_deriv[1] = %.10f, result[1] = %.10f\n", h_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(h_d_deriv[0], result[0]) && _test_nonpara_dcompare(h_d_deriv[1],result[1]));
    
    // test 3: 2-2 observations in two groups, adjust obs
    adjust_obs = 1;
    d[0] = 0;   x[0] = 1.0;
    d[1] = 0;   x[1] = -2.9;
    d[2] = 1;   x[2] = 4.0;
    d[3] = 1;   x[3] = 1.2;
    x_eval = 0.5;
    h_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (kernel_deriv((x[0] - x_eval)/bandwidth) + kernel_deriv((x[1] - x_eval)/bandwidth));
    h_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (kernel_deriv((x[2] - x_eval)/bandwidth) + kernel_deriv((x[3] - x_eval)/bandwidth));; // no observation with d=1
    result = nonpara_h_deriv_eval_filter(x_eval, x, n, stride, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("h_d_deriv[0] = %.10f, result[0] = %.10f\n", h_d_deriv[0], result[0]);
    // printf("h_d_deriv[1] = %.10f, result[1] = %.10f\n", h_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(h_d_deriv[0], result[0]) && _test_nonpara_dcompare(h_d_deriv[1],result[1]));
    
    // test 4: 2-2 observations in two groups, adjust obs, stride = 2
    adjust_obs = 1;
    stride = 2;
    double blank_value = 0.0;
    double *x_s = malloc(stride * n * sizeof(double));
    d[0] = 0;   x_s[0] = 1.0;
                x_s[1] = blank_value;
    d[1] = 1;   x_s[2] = -2.9;
                x_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;
                x_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 1.2;
                x_s[7] = blank_value;
    x_eval = 0.5;
    h_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (kernel_deriv((x_s[0] - x_eval)/bandwidth) + kernel_deriv((x_s[4] - x_eval)/bandwidth));
    h_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (kernel_deriv((x_s[2] - x_eval)/bandwidth) + kernel_deriv((x_s[6] - x_eval)/bandwidth));; // no observation with d=1
    result = nonpara_h_deriv_eval_filter(x_eval, x_s, n, stride, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("h_d_deriv[0] = %.10f, result[0] = %.10f\n", h_d_deriv[0], result[0]);
    // printf("h_d_deriv[1] = %.10f, result[1] = %.10f\n", h_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(h_d_deriv[0], result[0]) && _test_nonpara_dcompare(h_d_deriv[1],result[1]));

    // cleanup
    free(x);
    free(x_s);
    free(d);
    free(result);
    free(h_d_deriv);
    printf("Tests for `nonpara_h_deriv_eval_filter` completed: no issues found.\n");
}



/*  Evaluates q_deriv_b := (- scale / (bandwidth^2 * n_adjusted_b)) * sum_{i in [n]: filter(i)==b} (y_i^k) * K'((x_i-x_eval)/ bandwidth) for b=0,1, where
 - `scale` is real scalar;
 - n is the length of `x`;
 - x_i, y_i are elements  i=0,...,n-1 of `x` and `y`;
 - K' is the derivative `kernel_deriv` of a kernel;
 - filter(i) is 0 if `filter_func` evaluates to false, 1 otherwise;
 - `filter_func` takes and observation index and additional argument `filter_args` passed as void pointer.
 - n_adjusted_b is the number of observations i such that filter(i)=b if adjust_nobs is true, n_adjusted_b = n otherwise.
 
 Returns 2-long pointer `eval_value` with eval_value[0] = q_deriv_0 and eval_value[1] = q_deriv_1.
*/
double *nonpara_qk_deriv_eval_filter(double x_eval, double *x, unsigned long int n_x, unsigned long int stride_x, double *y, unsigned long int n_y, unsigned long int stride_y, double k, double (*kernel_deriv)(double), double bandwidth, double scale, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *)){
    assert(n_x == n_y);
    unsigned long int n = n_x;
    unsigned long int n_true = 0;
    unsigned long int n_false = 0;
    double eval_value_false, eval_value_true;
    eval_value_false = eval_value_true = 0;
    for (unsigned long int i=0; i<n; i++){
        if (filter_fn(i, filter_args)){  // filter true
            eval_value_true += (-scale) * pow(y[stride_y * i], k) * kernel_deriv((x[stride_x * i] - x_eval) / bandwidth) / (n * pow(bandwidth, 2));
            n_true++;
        } else {    // filter false
            eval_value_false += (-scale) * pow(y[stride_y * i], k) * kernel_deriv((x[stride_x * i] - x_eval) / bandwidth) / (n * pow(bandwidth, 2));
            n_false++;
        }
    }
    double adjustment_factor_nobs_true, adjustment_factor_nobs_false;
    if (adjust_nobs){
        adjustment_factor_nobs_true = (n + 0.0) / n_true;
        adjustment_factor_nobs_false = (n + 0.0) / n_false; 
    } else{
        adjustment_factor_nobs_false = adjustment_factor_nobs_true = 1;
    }
    double *eval_value = malloc(2*sizeof(double)); // entry 0 = filter false; entry 1 = filter true
    eval_value[0] = eval_value_false * adjustment_factor_nobs_false;
    eval_value[1] = eval_value_true * adjustment_factor_nobs_true;
    return eval_value;
}


// Auxiliary function for test_qk_deriv_eval_filter.
int _test_nonpara_qk_deriv_eval_filter_filter_func(unsigned long int i, void *filter_arg){
    int *d = (int *) filter_arg;
    return d[i];
}


void test_nonpara_qk_deriv_eval_filter(void){
    // test 1: only one group d=1
    int n = 4;
    int stride_x, stride_y;    
    stride_x = stride_y = 1;
    double bandwidth = pow(n, -1.0/5);
    double scale = -3.9;
    double (*kernel_deriv)(double);
    int (*filter_func)(unsigned long int, void *);
    kernel_deriv = nonpara_kernel_gauss;
    filter_func = _test_nonpara_qk_deriv_eval_filter_filter_func;
    int adjust_obs = 0;
    double power = 1;
    int *d = malloc(n*sizeof(int));  // filter variable with two categories
    double *x = malloc(n*sizeof(double));
    double *y = malloc(n*sizeof(double));
    d[0] = 1;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 1;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 1;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;
    double x_eval = 0.5;
    double *q_d_deriv = malloc(2*sizeof(double));
    q_d_deriv[0] = 0.0; // no observation with d=0
    q_d_deriv[1] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[0], power) * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y[1], power) * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y[2], power) * kernel_deriv((x[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * kernel_deriv((x[3] - x_eval)/bandwidth));
    double *result = nonpara_qk_deriv_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));
    
    // test 2: only one group d=0
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 0;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 0;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 0;   x[3] = 0.0;     y[3] = 50.4;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[0], power) * kernel_deriv((x[0]- x_eval)/bandwidth) +
                                      pow(y[1], power) * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y[2], power) * kernel_deriv((x[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * kernel_deriv((x[3] - x_eval)/bandwidth));
    q_d_deriv[1] = 0.0; // no observation with d=0
    result = nonpara_qk_deriv_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));

    // test 3: 2-2 observations in two groups
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 0;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 1;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[0], power) * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y[1], power) * kernel_deriv((x[1] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[2], power) * kernel_deriv((x[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * kernel_deriv((x[3] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));
    
    // test 3: 2-2 observations in two groups, adjust number of obs
    adjust_obs = 1;
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;
    d[1] = 1;   x[1] = -2.9;    y[1] = -30.5;
    d[2] = 0;   x[2] = 4.0;     y[2] = 100.45;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[0], power) * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y[2], power) * kernel_deriv((x[2] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[1], power) * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * kernel_deriv((x[3] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));

    // test 4: 2-2 observations in two groups, adjust number of obs, stride_x = 2
    adjust_obs = 1;
    stride_x = 2;
    double blank_value = 1000000;
    double *x_s = malloc(stride_x * n * sizeof(double));
    d[0] = 0;   x_s[0] = 1.0;     y[0] = 10.1;
                x_s[1] = blank_value;
    d[1] = 1;   x_s[2] = -2.9;    y[1] = -30.5;
                x_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;     y[2] = 100.45;
                x_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 0.0;     y[3] = 50.4;
                x_s[7] = blank_value;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[0], power) * kernel_deriv((x_s[0] - x_eval)/bandwidth) +
                                      pow(y[2], power) * kernel_deriv((x_s[4] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[1], power) * kernel_deriv((x_s[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * kernel_deriv((x_s[6] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_eval_filter(x_eval, x_s, n, stride_x, y, n, stride_y, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1], result[1]));

    // test 5: 2-2 observations in two groups, adjust number of obs, stride_y = 2
    stride_x = 1;
    stride_y = 2;
    double *y_s = malloc(stride_y * n * sizeof(double));
    d[0] = 0;   x[0] = 1.0;     y_s[0] = 10.1;
                                y_s[1] = blank_value;
    d[1] = 1;   x[1] = -2.9;    y_s[2] = -30.5;
                                y_s[3] = blank_value;
    d[2] = 0;   x[2] = 4.0;     y_s[4] = 100.45;
                                y_s[5] = blank_value;
    d[3] = 1;   x[3] = 0.0;     y_s[6] = 50.4;
                                y_s[7] = blank_value;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[0], power) * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y_s[4], power) * kernel_deriv((x[2] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[2], power) * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y_s[6], power) * kernel_deriv((x[3] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_eval_filter(x_eval, x, n, stride_x, y_s, n, stride_y, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1], result[1]));

    // test 6: 2-2 observations in two groups, adjust number of obs, stride_x = 2, stride_y = 2
    stride_x = 2;
    stride_y = 2;
    d[0] = 0;   x_s[0] = 1.0;          y_s[0] = 10.1;
                x_s[1] = blank_value;  y_s[1] = blank_value;
    d[1] = 1;   x_s[2] = -2.9;         y_s[2] = -30.5;
                x_s[3] = blank_value;  y_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;          y_s[4] = 100.45;
                x_s[5] = blank_value;  y_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 0.0;          y_s[6] = 50.4;
                x_s[7] = blank_value;  y_s[7] = blank_value;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[0], power) * kernel_deriv((x_s[0] - x_eval)/bandwidth) +
                                      pow(y_s[4], power) * kernel_deriv((x_s[4] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[2], power) * kernel_deriv((x_s[2] - x_eval)/bandwidth) + 
                                      pow(y_s[6], power) * kernel_deriv((x_s[6] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_eval_filter(x_eval, x_s, n, stride_x, y_s, n, stride_y, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1], result[1]));
    
    // cleanup
    free(x);
    free(x_s);
    free(d);
    free(y);
    free(y_s);
    free(result);
    free(q_d_deriv);    
    printf("Tests for `nonpara_qk_deriv_eval_filter` completed: no issues found.\n");
}


////////// Kernel derivative estimators with weights and filters


/*  Evaluates q_deriv_b := (- scale / (bandwidth^2 * n_adjusted_b)) * sum_{i in [n]: filter(i)==b} weight_i * (y_i^k) * K'((x_eval-x_i)/ bandwidth) for b=0,1, where
 - `scale` is real scalar;
 - n is the length of `x`;
 - x_i, y_i, weight_i are elements  i=0,...,n-1 of `x`, `y` and `weight`;
 - K' is the derivative `kernel_deriv` of a kernel;
 - filter(i) is 0 if `filter_func` evaluates to false, 1 otherwise;
 - `filter_func` takes and observation index and additional argument `filter_args` passed as void pointer.
 - n_adjusted_b is the number of observations i such that filter(i)=b if adjust_nobs is true, n_adjusted_b = n otherwise.
 
 Returns 2-long pointer `eval_value` with eval_value[0] = q_deriv_0 and eval_value[1] = q_deriv_1.
*/
double *nonpara_qk_deriv_weighted_eval_filter(double x_eval, double *x, unsigned long int n_x, unsigned long int stride_x, double *y, unsigned long int n_y, unsigned long int stride_y, double *weight, unsigned long int n_weight, unsigned long int stride_weight, double k, double (*kernel_deriv)(double), double bandwidth, double scale, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *)){
    assert( (n_x == n_y) && (n_y == n_weight));
    unsigned long int n = n_x;
    unsigned long int n_true = 0;
    unsigned long int n_false = 0;
    double eval_value_false, eval_value_true;
    eval_value_false = eval_value_true = 0;
    for (unsigned long int i=0; i<n; i++){
        if (filter_fn(i, filter_args)){  // filter true
            eval_value_true += (-scale) * pow(y[stride_y * i], k) * weight[stride_weight * i] * kernel_deriv((x[stride_x * i] - x_eval) / bandwidth) / (n * pow(bandwidth, 2));
            n_true++;
        } else {    // filter false
            eval_value_false += (-scale) * pow(y[stride_y * i], k) * weight[stride_weight * i] * kernel_deriv((x[stride_x * i] - x_eval) / bandwidth) / (n * pow(bandwidth, 2));
            n_false++;
        }
    }
    double adjustment_factor_nobs_true, adjustment_factor_nobs_false;
    if (adjust_nobs){
        adjustment_factor_nobs_true = (n + 0.0) / n_true;
        adjustment_factor_nobs_false = (n + 0.0) / n_false; 
    } else{
        adjustment_factor_nobs_false = adjustment_factor_nobs_true = 1;
    }
    double *eval_value = malloc(2*sizeof(double)); // entry 0 = filter false; entry 1 = filter true
    eval_value[0] = eval_value_false * adjustment_factor_nobs_false;
    eval_value[1] = eval_value_true * adjustment_factor_nobs_true;
    return eval_value;   
}



// Auxiliary function for test_qk_deriv_eval_filter.
int _test_nonpara_qk_deriv_weighted_eval_filter_filter_func(unsigned long int i, void *filter_arg){
    int *d = (int *) filter_arg;
    return d[i];
}


void test_nonpara_qk_deriv_weighted_eval_filter(void){
    // test 1: only one group d=1
    int n = 4;
    int stride_x, stride_y, stride_w;    
    stride_x = stride_y = stride_w = 1;
    double bandwidth = pow(n, -1.0/5);
    double scale = -3.9;
    double (*kernel_deriv)(double);
    int (*filter_func)(unsigned long int, void *);
    kernel_deriv = nonpara_kernel_gauss;
    filter_func = _test_nonpara_qk_deriv_weighted_eval_filter_filter_func;
    int adjust_obs = 0;
    double power = 1;
    int *d = malloc(n*sizeof(int));  // filter variable with two categories
    double *x = malloc(n*sizeof(double));
    double *y = malloc(n*sizeof(double));
    double *w = malloc(n*sizeof(double));  // weights
    d[0] = 1;   x[0] = 1.0;     y[0] = 10.1;    w[0] = -10.5;
    d[1] = 1;   x[1] = -2.9;    y[1] = -30.5;   w[1] = 30.2;
    d[2] = 1;   x[2] = 4.0;     y[2] = 100.45;  w[2] = 100.15;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;    w[3] = -2.42;
    double x_eval = 0.5;
    double *q_d_deriv = malloc(2*sizeof(double));
    q_d_deriv[0] = 0.0; // no observation with d=0
    q_d_deriv[1] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[0], power) * w[0] * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y[1], power) * w[1] * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y[2], power) * w[2] * kernel_deriv((x[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * w[3] * kernel_deriv((x[3] - x_eval)/bandwidth));
    double *result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, w, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));
    
    // test 2: only one group d=0
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;    w[0] = -10.5;
    d[1] = 0;   x[1] = -2.9;    y[1] = -30.5;   w[1] = 30.2;
    d[2] = 0;   x[2] = 4.0;     y[2] = 100.45;  w[2] = 100.15;
    d[3] = 0;   x[3] = 0.0;     y[3] = 50.4;    w[3] = -2.42;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[0], power) * w[0] * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y[1], power) * w[1] * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y[2], power) * w[2] * kernel_deriv((x[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * w[3] * kernel_deriv((x[3] - x_eval)/bandwidth));
    q_d_deriv[1] = 0.0; // no observation with d=0
    result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, w, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));

    // test 3: 2-2 observations in two groups
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;    w[0] = -10.5;
    d[1] = 0;   x[1] = -2.9;    y[1] = -30.5;   w[1] = 30.2;
    d[2] = 1;   x[2] = 4.0;     y[2] = 100.45;  w[2] = 100.15;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;    w[3] = -2.42;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[0], power) * w[0] * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y[1], power) * w[1] * kernel_deriv((x[1] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/(n * pow(bandwidth, 2))) * (pow(y[2], power) * w[2] * kernel_deriv((x[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * w[3] * kernel_deriv((x[3] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, w, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));
    
    // test 3: 2-2 observations in two groups, adjust number of obs
    adjust_obs = 1;
    d[0] = 0;   x[0] = 1.0;     y[0] = 10.1;    w[0] = -10.5;
    d[1] = 1;   x[1] = -2.9;    y[1] = -30.5;   w[1] = 30.2;
    d[2] = 0;   x[2] = 4.0;     y[2] = 100.45;  w[2] = 100.15;
    d[3] = 1;   x[3] = 0.0;     y[3] = 50.4;    w[3] = -2.42;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[0], power) * w[0] * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y[2], power) * w[2] * kernel_deriv((x[2] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[1], power) * w[1] * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * w[3] * kernel_deriv((x[3] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x, n, stride_x, y, n, stride_y, w, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1],result[1]));

    // test 4: 2-2 observations in two groups, adjust number of obs, stride_x = 2
    adjust_obs = 1;
    stride_x = 2;
    double blank_value = 1000000;
    double *x_s = malloc(stride_x * n * sizeof(double));
    d[0] = 0;   x_s[0] = 1.0;     y[0] = 10.1;   w[0] = -10.5;
                x_s[1] = blank_value;
    d[1] = 1;   x_s[2] = -2.9;    y[1] = -30.5;  w[1] = 30.2;
                x_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;     y[2] = 100.45; w[2] = 100.15;
                x_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 0.0;     y[3] = 50.4;   w[3] = -2.42;
                x_s[7] = blank_value;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[0], power) * w[0] * kernel_deriv((x_s[0] - x_eval)/bandwidth) +
                                      pow(y[2], power) * w[2] * kernel_deriv((x_s[4] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y[1], power) * w[1] * kernel_deriv((x_s[2] - x_eval)/bandwidth) + 
                                      pow(y[3], power) * w[3] * kernel_deriv((x_s[6] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x_s, n, stride_x, y, n, stride_y, w, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1], result[1]));

    // test 5: 2-2 observations in two groups, adjust number of obs, stride_y = 2
    stride_x = 1;
    stride_y = 2;
    double *y_s = malloc(stride_y * n * sizeof(double));
    d[0] = 0;   x[0] = 1.0;     y_s[0] = 10.1;          w[0] = -10.5;
                                y_s[1] = blank_value;   
    d[1] = 1;   x[1] = -2.9;    y_s[2] = -30.5;         w[1] = 30.2;
                                y_s[3] = blank_value;   
    d[2] = 0;   x[2] = 4.0;     y_s[4] = 100.45;        w[2] = 100.15;
                                y_s[5] = blank_value;
    d[3] = 1;   x[3] = 0.0;     y_s[6] = 50.4;          w[3] = -2.42;
                                y_s[7] = blank_value;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[0], power) * w[0] * kernel_deriv((x[0] - x_eval)/bandwidth) +
                                      pow(y_s[4], power) * w[2] * kernel_deriv((x[2] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[2], power) * w[1] * kernel_deriv((x[1] - x_eval)/bandwidth) + 
                                      pow(y_s[6], power) * w[3] * kernel_deriv((x[3] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x, n, stride_x, y_s, n, stride_y, w, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1], result[1]));

    // test 6: 2-2 observations in two groups, adjust number of obs, stride_x = 2, stride_y = 2
    stride_x = 2;
    stride_y = 2;
    d[0] = 0;   x_s[0] = 1.0;          y_s[0] = 10.1;           w[0] = -10.5;
                x_s[1] = blank_value;  y_s[1] = blank_value;   
    d[1] = 1;   x_s[2] = -2.9;         y_s[2] = -30.5;          w[1] = 30.2;
                x_s[3] = blank_value;  y_s[3] = blank_value;
    d[2] = 0;   x_s[4] = 4.0;          y_s[4] = 100.45;         w[2] = 100.5;
                x_s[5] = blank_value;  y_s[5] = blank_value;
    d[3] = 1;   x_s[6] = 0.0;          y_s[6] = 50.4;           w[3] = -2.42;
                x_s[7] = blank_value;  y_s[7] = blank_value;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[0], power) * w[0] * kernel_deriv((x_s[0] - x_eval)/bandwidth) +
                                      pow(y_s[4], power) * w[2] * kernel_deriv((x_s[4] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[2], power) * w[1] * kernel_deriv((x_s[2] - x_eval)/bandwidth) + 
                                      pow(y_s[6], power) * w[3] * kernel_deriv((x_s[6] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x_s, n, stride_x, y_s, n, stride_y, w, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1], result[1]));
    
    // test 7: 2-2 observations in two groups, adjust number of obs, stride_x = 2, stride_y = 2, stride_w = 2
    stride_w = 2;
    double blank_value_w = 9000000;
    double *w_s = malloc(stride_w * n * sizeof(double));
    d[0] = 0;   x_s[0] = 1.0;          y_s[0] = 10.1;           w_s[0] = -10.5;
                x_s[1] = blank_value;  y_s[1] = blank_value;    w_s[1] = blank_value_w;
    d[1] = 1;   x_s[2] = -2.9;         y_s[2] = -30.5;          w_s[2] = 30.2;
                x_s[3] = blank_value;  y_s[3] = blank_value;    w_s[3] = blank_value_w; 
    d[2] = 0;   x_s[4] = 4.0;          y_s[4] = 100.45;         w_s[4] = 100.5;
                x_s[5] = blank_value;  y_s[5] = blank_value;    w_s[5] = blank_value_w;
    d[3] = 1;   x_s[6] = 0.0;          y_s[6] = 50.4;           w_s[6] = -2.42;
                x_s[7] = blank_value;  y_s[7] = blank_value;    w_s[7] = blank_value_w;
    x_eval = 0.5;
    q_d_deriv[0] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[0], power) * w_s[0] * kernel_deriv((x_s[0] - x_eval)/bandwidth) +
                                      pow(y_s[4], power) * w_s[4] * kernel_deriv((x_s[4] - x_eval)/bandwidth));
    q_d_deriv[1] = (-scale/((n / 2.0) * pow(bandwidth, 2))) * (pow(y_s[2], power) * w_s[2] * kernel_deriv((x_s[2] - x_eval)/bandwidth) + 
                                      pow(y_s[6], power) * w_s[6] * kernel_deriv((x_s[6] - x_eval)/bandwidth));
    result = nonpara_qk_deriv_weighted_eval_filter(x_eval, x_s, n, stride_x, y_s, n, stride_y, w_s, n, stride_w, power, kernel_deriv, bandwidth, scale, adjust_obs, (void *) d, filter_func);
    // printf("q_d_deriv[0] = %f, result[0] = %f\n", q_d_deriv[0], result[0]);
    // printf("q_d_deriv[1] = %f, result[1] = %f\n", q_d_deriv[1], result[1]);
    assert( _test_nonpara_dcompare(q_d_deriv[0], result[0]) && _test_nonpara_dcompare(q_d_deriv[1], result[1]));
    
    // cleanup
    free(x);
    free(x_s);
    free(d);
    free(y);
    free(y_s);
    free(result);
    free(q_d_deriv);
    printf("Tests for `nonpara_qk_deriv_weighted_eval_filter` completed: no issues found.\n");
}



void test_nonpara(void){
    test_nonpara_h_eval_filter();
    test_nonpara_qk_eval_filter();
    test_nonpara_h_deriv_eval_filter();
    test_nonpara_qk_deriv_eval_filter();
    test_nonpara_qk_deriv_weighted_eval_filter();
    printf("=== Tests: all tests for `nonpara` completed: no issues found.\n");
}