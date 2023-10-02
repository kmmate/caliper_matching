/* Header for propscore.c */

#ifndef PROPSCORE_H_
#define PROPSCORE_H_

// includes
#include <math.h>
#include <gsl/gsl_cdf.h>


// functions
// logit
double propscore_g_logit(double x);
double propscore_ginv_logit(double p);
double propscore_gderiv_logit(double x);
double propscore_ginvderiv_logit(double p);
// probit
double propscore_g_probit(double x);
double propscore_ginv_probit(double p);
double propscore_gderiv_probit(double x);
double propscore_ginvderiv_probit(double p);



// variables




#endif // PROPSCORE_H_