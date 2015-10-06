#include <math.h>
#include <gsl/gsl_cdf.h>
#include "stats.h"

double lrt(double log10lik_null, double log10lik_alt, int nparam_null, int nparam_alt)
{
    double teststat = - 2 * log10lik_null / log10(exp(1)) + 2 * log10lik_alt / log10(exp(1));
    int df = nparam_alt - nparam_null;
    return 1 - gsl_cdf_chisq_P(teststat, df);
}

double bic(double log10lik, int nparam, int ndata)
{
    return -2 * log10lik / log10(exp(1)) + nparam * log(ndata);
}

double aic(double log10lik, int nparam)
{
    return -2 * log10lik / log10(exp(1)) + 2 * nparam;
}
