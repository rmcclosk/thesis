#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "smc.h"

#define NT 200
#define M 15
#define POPSIZE 100
#define END_EPSILON 0.01
#define alpha 0.9

double ess(const double *W);
double calc_new_epsilon(int nparam, double *W, void **X, metric d, 
        double old_epsilon, const void *data);

// see Del Moral et al. 2012: An adaptive sequential Monte Carlo method for
// approximate Bayesian computation
void abc_smc(int nparam, real_rv *prior_rvs, sampler data_sampler, metric d,
             const void *data)
{
    int i, j;
    double *theta = malloc(nparam * POPSIZE * sizeof(double));
    void **X = malloc(M * POPSIZE * sizeof(void*));
    double *W = malloc(POPSIZE * sizeof(double));
    double prev_ess = POPSIZE, prev_epsilon, cur_epsilon = INFINITY;
    double num, denom;

    // step 0: sample initial parameters from prior distributions, and initial
    // data from parameters
    for (i = 0; i < POPSIZE; ++i)
    {
        W[i] = 1.0 / POPSIZE;
        for (j = 0; j < nparam; ++j)
        {
            theta[i * nparam + j] = prior_rvs[j]();
        }
        for (j = 0; j < M; ++j)
        {
            X[i * M + j] = data_sampler(&theta[i * nparam]);
        }
    }

    // step 1: find new epsilon
    if (cur_epsilon == END_EPSILON)
        return;
    prev_epsilon = cur_epsilon;
    cur_epsilon = calc_new_epsilon(nparam, W, X, d, prev_epsilon, data);

    // step 2: resample the particles
    free(theta);
    free(W);
    free(X);
}

/* Private. */

#define epsilon_tol 0.001

double calc_new_epsilon(int nparam, double *W, void **X, metric d, 
                        double old_epsilon, const void *data)
{
    int i;
    double sum, num = 0, denom = 0;
    double prev_ess = ess(W);
    double new_epsilon, lbound = 0, ubound = old_epsilon;
    double *new_W = malloc(POPSIZE * sizeof(double));

    // ESS(W, epsilon) is an increasing function of epsilon
    while (ubound - lbound > epsilon_tol)
    {
        new_epsilon = (lbound + ubound) / 2.;
        for (i = 0; i < M; ++i)
        {
            if (d(X[i], data) < old_epsilon)
                denom += 1;
            if (d(X[i], data) < new_epsilon)
                num += 1;
        }
    
        sum = 0;
        for (i = 0; i < POPSIZE; ++i)
        {
            new_W[i] = W[i] * num / denom;
            sum += new_W[i];
        }
        for (i = 0; i < POPSIZE; ++i)
        {
            new_W[i] /= sum;
        }
    
        if (fabs(ess(new_W) / prev_ess) > alpha)
            ubound = new_epsilon;
        else if (fabs(ess(new_W) / prev_ess) < alpha)
            lbound = new_epsilon;
        else
            break;
    }
    memcpy(W, new_W, POPSIZE * sizeof(double));
    free(new_W);
    return new_epsilon;
}

double ess(const double *W)
{
    int i;
    double sum = 0;
    for (i = 0; i < POPSIZE; ++i) {
        sum += W[i] * W[i];
    }
    return 1.0 / sum;
}
