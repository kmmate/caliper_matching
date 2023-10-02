/* Header file for nonpara.c */

#ifndef NONPARA_H_
#define NONPARA_H_

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#define NONPARA_TEST_DOUBLE_TOLERANCE_NPOWER 9

// functions
double nonpara_kernel_gauss(double u);
double nonpara_kernelderiv_gauss(double u);
double nonpara_h_eval(double x_eval, double *x, unsigned long int n, double (*kernel)(double), double bandwidth);
double nonpara_qk_eval(double x_eval, double *x, unsigned long int n_x, double *y, unsigned long int n_y, double k, double (*kernel)(double), double bandwidth);
// filters
double *nonpara_h_eval_filter(double x_eval, double *x, unsigned long int n, unsigned long int stride, double (*kernel)(double), double bandwidth, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *));
double *nonpara_qk_eval_filter(double x_eval, double *x, unsigned long int n_x, unsigned long int stride_x, double *y, unsigned long int n_y,unsigned long int stride_y, double k, double (*kernel)(double), double bandwidth, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *));
double *nonpara_h_deriv_eval_filter(double x_eval, double *x, unsigned long int n, unsigned long int stride, double (*kernel_deriv)(double), double bandwidth, double scale, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *));
double *nonpara_qk_deriv_eval_filter(double x_eval, double *x, unsigned long int n_x, unsigned long int stride_x, double *y, unsigned long int n_y, unsigned long int stride_y, double k, double (*kernel_deriv)(double), double bandwidth, double scale, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *));
// filters ans weights
double *nonpara_qk_deriv_weighted_eval_filter(double x_eval, double *x, unsigned long int n_x, unsigned long int stride_x, double *y, unsigned long int n_y, unsigned long int stride_y, double *weight, unsigned long int n_weight, unsigned long int stride_weight, double k, double (*kernel_deriv)(double), double bandwidth, double scale, int adjust_nobs, void *filter_args, int (*filter_fn)(unsigned long int, void *));

// tests
void test_nonpara(void);

#endif