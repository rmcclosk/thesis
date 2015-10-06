#ifndef STATS_H
#define STATS_H

typedef enum {
    AIC,
    BIC,
    LRT
} model_selector;

/** Calculate a likelihood ratio test.
 *
 * \param[in] log10lik_null log10 likelihood of null model
 * \param[in] log10lik_alt log10 likelihood of alternative model
 * \param[in] nparam_null number of parameters of null model
 * \param[in] nparam_alt number of parameters of alternative model
 * \return p-value for support of the alternative model
 */
double lrt(double log10lik_null, double log10lik_alt, int nparam_null, int nparam_alt);

/** Calculate the Bayesian information criterion for a model fit.
 * 
 * \param[in] log10lik log10 likelihood of the fitted model
 * \param[in] nparam number of parameters of the model
 * \param[in] ndata number of data points of the model
 * \return the BIC
 */
double bic(double log10lik, int nparam, int ndata);

/** Calculate the Akaike information criterion for a model fit.
 *
 * \param[in] log10lik log 10 likelihood of the fitted model
 * \param[in] nparam number of parameters of the model
 */
double aic(double log10lik, int nparam);

#endif
