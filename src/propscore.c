//////////////////////////////////////////////////////////////////////////
//
// * FILE: cm.c
// * DESCRIPTION:
//    Propensity score model.
// * AUTHOR: Mate Kormos
// * LAST REVISED: 11/oct/2022
// * NOTE: To add a new model implement its functions similarly to logit
//          and probit and add them to the header file. Also modify the
//          cm.c file as described therein.
//
//////////////////////////////////////////////////////////////////////////

#include "propscore.h"


// logit g
double propscore_g_logit(double x){
    return 1.0 / (1.0 + exp(-x));
}
// logit g inverse
double propscore_ginv_logit(double p){
    return log(p / (1.0 - p));
}
// logit derivative of g
double propscore_gderiv_logit(double x){
    return exp(-x) / pow(1.0 + exp(-x), 2);
}

// logit (g^{-1})'(p)
double propscore_ginvderiv_logit(double p){
    //return 1.0 / propscore_gderiv_logit(propscore_ginv_logit(p));
    return 1.0 / (p * (1.0 - p));  // analytically computed
}

// probit g
double propscore_g_probit(double x){
    return 0.5 * (1 + erf(x / M_SQRT2 ));
}
// probit g inverse
double propscore_ginv_probit(double p){
    return gsl_cdf_ugaussian_Pinv(p);
}
// probit g derivative
double propscore_gderiv_probit(double x){
    return (1.0 / (sqrt(M_PI) * M_SQRT2)) * exp( - pow(x, 2));
}
// probit (g^{-1})'(p)
double propscore_ginvderiv_probit(double p){
    return 1.0 / propscore_gderiv_probit(propscore_ginv_probit(p));
}