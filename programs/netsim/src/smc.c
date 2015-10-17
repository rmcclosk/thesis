#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_roots.h>
#include "smc.h"
#include "util.h"

#define RESIZE_AMOUNT 100
#define NTHREAD 1

struct objfun_params {
    double *X;
    double *W;
    double *new_W;
    double epsilon;
    const smc_config *config;
};

struct perturb_workspace {
    gsl_rng *rng;
    const smc_config *config;
    const smc_functions *functions;
    double *W;
    double *theta;
    double *X;
    double epsilon;
    const void *data;
    char *fdbk;
    double *new_theta;
    double *new_X;
    char *z;
    int nparticle;
};

void resample(const smc_config config, double *W, double *theta);
double next_epsilon(struct objfun_params *objpar);
double epsilon_objfun(double epsilon, void *params);
double ess(const double *W, int n);
double perturb(struct perturb_workspace *ws);

// see Del Moral et al. 2012: An adaptive sequential Monte Carlo method for
// approximate Bayesian computation
smc_result *abc_smc(const smc_config config, const smc_functions functions,
                    int seed, const void *data)
{
    int i, j;
    double *X = malloc(config.nsample * config.nparticle * sizeof(double));
    double *W = malloc(config.nparticle * sizeof(double));
    char *z = malloc(config.dataset_size);
    double epsilon = 1;
    size_t new_size;
    gsl_rng *rng = set_seed(seed);

    struct objfun_params objpar = { 
        .X = X, 
        .W = W,
        .new_W = malloc(config.nparticle * sizeof(double)),
        .epsilon = epsilon, 
        .config = &config 
    };

    smc_result *res = malloc(sizeof(smc_result));
    res->epsilon = malloc(RESIZE_AMOUNT * sizeof(double));
    res->acceptance_rate = malloc(RESIZE_AMOUNT * sizeof(double));
    res->theta = malloc(config.nparticle * config.nparam * sizeof(double));
    res->epsilon[0] = epsilon;

    struct perturb_workspace ws = {
        .config = &config,
        .functions = &functions,
        .epsilon = epsilon,
        .data = data,
        .nparticle = config.nparticle
    };

    struct perturb_workspace *ws2 = malloc(NTHREAD * sizeof(struct perturb_workspace));
    for (i = 0; i < NTHREAD; ++i)
    {
        memcpy(&ws2[i], &ws, sizeof(struct perturb_workspace));
        ws2[i].rng = set_seed(seed+i+1);
        ws2[i].W = &W[config.nparticle * i / NTHREAD];
        ws2[i].theta = &res->theta[config.nparticle * config.nparam * i / NTHREAD];
        ws2[i].X = &X[config.nparticle * config.nsample * config.dataset_size * i / NTHREAD];
        ws2[i].new_theta = malloc(config.nparticle * config.nparam * sizeof(double));
        ws2[i].fdbk = malloc(config.feedback_size);
        ws2[i].new_X = malloc(config.nsample * sizeof(double));
        ws2[i].z = malloc(config.dataset_size);
    }

    // step 0: sample particles from prior
    for (i = 0; i < config.nparticle; ++i)
    {
        functions.sample_from_prior(rng, &res->theta[i * config.nparam]);
        W[i] = 1. / config.nparticle;
        for (j = 0; j < config.nsample; ++j)
        {
            functions.sample_dataset(rng, &res->theta[i * config.nparam], z);
            X[i * config.nsample + j] = functions.distance(z, data);
        }
    }

    res->niter = 1;
    while (epsilon != config.final_epsilon)
    {
        fprintf(stderr, "%d\t%f\n", res->niter, epsilon);

        // step 1: update epsilon
        epsilon = fmax(next_epsilon(&objpar), config.final_epsilon);
        objpar.epsilon = epsilon;
        for (i = 0; i < NTHREAD; ++i)
            ws2[i].epsilon = epsilon;
        res->epsilon[res->niter] = epsilon;
    
        // step 2: resample particles according to their weights
        if (ess(W, config.nparticle) < config.ess_tolerance) {
            resample(config, W, res->theta);
        }
    
        // step 3: perturb particles
        functions.feedback(res->theta, config.nparticle, ws2[0].fdbk);
        for (i = 1; i < NTHREAD; ++i)
            memcpy(ws2[i].fdbk, ws2[0].fdbk, config.feedback_size);
        res->acceptance_rate[res->niter] = perturb(&ws2[0]);
        ++res->niter;

        // allocate more space
        if (res->niter % RESIZE_AMOUNT == 0)
        {
            new_size = RESIZE_AMOUNT * (res->niter / RESIZE_AMOUNT + 1) * sizeof(double);
            res->epsilon = safe_realloc(res->epsilon, new_size);
            res->acceptance_rate = safe_realloc(res->acceptance_rate, new_size);
        }
    }

    // finally, sample from the estitmated posterior
    resample(config, W, res->theta);

    // clean up
    free(W);
    free(X);
    free(z);

    for (i = 0; i < NTHREAD; ++i)
    {
        free(ws2[i].new_theta);
        free(ws2[i].fdbk);
        free(ws2[i].new_X);
        free(ws2[i].z);
    }
    free(ws2);

    return res;
}

void smc_result_free(smc_result *r)
{
    free(r->epsilon);
    free(r->acceptance_rate);
    free(r->theta);
}

/* Private. */

double next_epsilon(struct objfun_params *objpar)
{
    int status, iter = 0, max_iter = 10000;
    double r, x_lo = 0, x_hi = objpar->epsilon;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (gsl_root_fsolver_bisection);
    gsl_function F = { .function = &epsilon_objfun, .params = objpar };
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);

    do
    {
        status = gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (x_lo, x_hi, 0, objpar->config->step_tolerance);
        ++iter;
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    // TODO: handle running out of iterations?
    if (iter == max_iter) {
        fprintf(stderr, "Warning: hit max iterations solving for next epsilon\n");
        epsilon_objfun(objpar->epsilon * objpar->config->quality, objpar);
    }
    memcpy(objpar->W, objpar->new_W, objpar->config->nparticle * sizeof(double));

    gsl_root_fsolver_free(s);
    return r;
}

double epsilon_objfun(double epsilon, void *params)
{
    int i, j; 
    double num, denom;
    struct objfun_params *par = (struct objfun_params *) params;
    const smc_config *config = par->config;
    double obj, sum = 0;

    for (i = 0; i < config->nparticle; ++i)
    {
        num = 0; denom = 0;
        for (j = 0; j < config->nsample; ++j)
        {
            num += par->X[i * config->nsample + j] < epsilon;
            denom += par->X[i * config->nsample + j] < par->epsilon;
        }
        if (num == denom) {
            par->new_W[i] = par->W[i];
        }
        else {
            par->new_W[i] = par->W[i] * num / denom;
        }
        sum += par->new_W[i];
    }
    for (i = 0; i < config->nparticle; ++i) 
    {
        if (sum == 0) {
            num = 0; denom = 0;
            for (j = 0; j < config->nsample; ++j)
            {
                num += par->X[i * config->nsample + j] < epsilon;
                denom += par->X[i * config->nsample + j] < par->epsilon;
            }
        }
        par->new_W[i] /= sum;
    }

    if (epsilon == 0 || sum == 0)
        obj = -1;
    else
        obj = ess(par->new_W, config->nparticle) - par->config->quality * ess(par->W, config->nparticle);
    return obj;
}

double ess(const double *W, int n)
{
    int i;
    double sum = 0;
    for (i = 0; i < n; ++i) {
        sum += W[i] * W[i];
    }
    return 1.0 / sum;
}

void resample(const smc_config config, double *W, double *theta)
{
    int i, wcur;
    double r, wsum;
    double *new_theta = malloc(config.nparticle * config.nparam * sizeof(double));

    for (i = 0; i < config.nparticle; ++i)
    {
        wcur = 0;
        wsum = W[wcur];
        r = (double) rand() / (double) RAND_MAX;
        while (r > wsum && wcur < config.nparticle) {
            wsum += W[++wcur];
        }
        memcpy(&new_theta[i * config.nparam], 
               &theta[wcur * config.nparam], 
               config.nparam * sizeof(double));
    }

    for (i = 0; i < config.nparticle; ++i) {
        W[i] = 1.0 / config.nparticle;
    }

    memcpy(theta, new_theta, config.nparticle * config.nparam * sizeof(double));
    free(new_theta);
}

double perturb(struct perturb_workspace *ws)
{
    int i, j;
    double mh_ratio, old_nbhd, new_nbhd;
    double *cur_theta, *prev_theta;
    double accept = 0, alive = 0;

    // gather feedback from current particle population
    memcpy(ws->new_theta, ws->theta, ws->nparticle * ws->config->nparam * sizeof(double));

    for (i = 0; i < ws->nparticle; ++i)
    {
        if (ws->W[i] == 0)
            continue;
        ++alive;

        cur_theta = &ws->new_theta[i * ws->config->nparam];
        prev_theta = &ws->theta[i * ws->config->nparam];

        // perturb the particle
        ws->functions->propose(ws->rng, cur_theta, ws->fdbk);

        // proposal ratio
        mh_ratio = ws->functions->proposal_density(cur_theta, prev_theta, ws->fdbk) /
                   ws->functions->proposal_density(prev_theta, cur_theta, ws->fdbk);

        // prior ratio
        mh_ratio *= ws->functions->prior_density(cur_theta) / 
                    ws->functions->prior_density(prev_theta);

        if (mh_ratio == 0)
            continue;

        // sample new datasets
        for (j = 0; j < ws->config->nsample; ++j)
        {
            ws->functions->sample_dataset(ws->rng, cur_theta, ws->z);
            ws->new_X[j] = ws->functions->distance(ws->z, ws->data);
        }

        // SMC approximation to likelihood ratio
        old_nbhd = 0; new_nbhd = 0;
        for (j = 0; j < ws->config->nsample; ++j) {
            old_nbhd += ws->X[i * ws->config->nsample + j] < ws->epsilon;
            new_nbhd += ws->new_X[j] < ws->epsilon;
        }
        mh_ratio *= new_nbhd / old_nbhd;

        // accept or reject the proposal
        if ((double) rand() / (double) RAND_MAX < mh_ratio)
        {
            ++accept;
            memcpy(prev_theta, cur_theta, ws->config->nparam * sizeof(double));
            memcpy(&ws->X[i * ws->config->nsample], ws->new_X, 
                   ws->config->nsample * sizeof(double));
        }
    }
    return accept / alive;
}
