#ifndef SMC_H
#define SMC_H

#include <gsl/gsl_rng.h>

typedef struct {
    int nparam; /**< number of parameters in the model */
    int nparticle; /**< number of particles used to approximate the posterior */
    int nsample; /**< number of sampled data points per particle */
    int ess_tolerance; /**< ESS below this value triggers resampling */
    double final_epsilon; /**< tolerance level to end at */
    double quality; /**< between 0 and 1, where 0 is fast and coarse, 1 is slow and accurate */
    double step_tolerance; /**< tolerance for bisection solution of next epsilon */
    size_t dataset_size; /**< size of data objects */
    size_t feedback_size; /**< size of feedback objects */
} smc_config;

typedef struct {
    void   (*sample_from_prior) (gsl_rng *rng, double *);
    double (*prior_density)     (const double *);
    void   (*propose)           (gsl_rng *rng, double *, const void *);
    double (*proposal_density)  (const double *, const double *, const void *);
    void   (*sample_dataset)    (gsl_rng *rng, const double *, void *);
    double (*distance)          (const void *, const void *);
    void   (*feedback)          (const double *, int, void *);
    void   (*destroy_dataset)   (void *);
} smc_functions;

typedef struct {
    int niter;
    double *epsilon;
    double *acceptance_rate;
    double *theta;
} smc_result;

/** Perform ABC-SMC.
 *
 * This implements the adaptive tolerance algorithm from DelMoral et al. 2012.
 * The object returned must be passed to smc_result_free() when you are done
 * with it.
 *
 * \param[in] config control parameters for the algorithm
 * \param[in] functions functions to use in the algorithm
 * \param[in] seed random seed (if negative, use time)
 * \panam[in] nthread number of threads to use
 * \param[in] data true data
 * \return an smc_result object with the posterior distribution (theta) and
 * other information about the run
 */
smc_result *abc_smc(const smc_config config, const smc_functions functions,
                    int seed, int nthread, const void *data);

/** Free an smc_result object.
 *
 * \param[in] r object to free
 */
void smc_result_free(smc_result *r);

#endif
