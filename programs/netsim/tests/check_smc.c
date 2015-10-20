#include <check.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include "../src/smc.h"
#include "../src/util.h"

Suite *smc_suite(void);

void toy_sample_from_prior (gsl_rng *rng, double *theta)
{
    *theta = gsl_ran_flat(rng, -10, 10);
}

double toy_prior_density(const double *theta)
{
    if (*theta < -10 || *theta > 10)
        return 0.;
    return 1./20.;
}

void toy_propose(gsl_rng *rng, double *theta, const void *params)
{
    double var = *((double *) params);
    *theta += gsl_ran_gaussian(rng, sqrt(2*var));
}

double toy_proposal_density(const double *from, const double *to, const void *params)
{
    double var = * ((double *) params);
    return gsl_ran_gaussian_pdf(*to - *from, sqrt(2*var));
}

void toy_sample_dataset(gsl_rng *rng, const double *theta, void *X)
{
    double x;
    if (rand() < RAND_MAX / 2) {
        x = (*theta + gsl_ran_gaussian(rng, 1)); 
    }
    else {
        x = (*theta + gsl_ran_gaussian(rng, 0.1));
    }
    memcpy(X, &x, sizeof(double));
}

double toy_distance(const void *x, const void *y)
{
    return fabs(*(double *) x - *(double *) y);
}

void toy_feedback(const double *theta, int nparticle, void *params)
{
    double var = gsl_stats_variance(theta, 1, nparticle);
    memcpy(params, &var, sizeof(double));
}

void toy_destroy_dataset(void *z)
{
    return;
}

smc_functions toy_functions = {
    .sample_from_prior = toy_sample_from_prior,
    .prior_density = toy_prior_density,
    .propose = toy_propose,
    .proposal_density = toy_proposal_density,
    .sample_dataset = toy_sample_dataset,
    .distance = toy_distance,
    .feedback = toy_feedback,
    .destroy_dataset = toy_destroy_dataset
};

void plot(const double *x, const double *y, int nx, const char *pdf, 
          const char *plot_cmd)
{
    int i, fd;
    char fn[BUFSIZ], cmd[BUFSIZ];

    sprintf(fn, "%s/smcXXXXXX", P_tmpdir);
    fd = mkstemp(fn);
    for (i = 0; i < nx; ++i) {
        if (y) {
            dprintf(fd, "%e\t%e\n", x[i], y[i]);
        } 
        else {
            dprintf(fd, "%e\n", x[i]);
        }
    }
    close(fd);
    
    sprintf(cmd, "R --vanilla --silent -e 'd <- read.table(\"%s\")' ", fn);
    sprintf(cmd, "%s -e 'pdf(\"%s\")' ", cmd, pdf);
    sprintf(cmd, "%s -e '%s' ", cmd, plot_cmd);
    sprintf(cmd, "%s -e 'dev.off()'\n", cmd);
    i = system(cmd);

    unlink(fn);
}

START_TEST (test_smc_toy_steps)
{
    int i;
    double y = 0;
    smc_config toy_config = {
        .nparam = 1,
        .nparticle = 10000,
        .nsample = 1,
        .ess_tolerance = 5000,
        .final_epsilon = 0.01,
        .quality = 0.95,
        .step_tolerance = 1e-5,
        .dataset_size = sizeof(double),
        .feedback_size = sizeof(double)
    };

    smc_result *res = abc_smc(toy_config, toy_functions, 0, (void *) &y);

    fprintf(stderr, "%d steps\n", res->niter);
    plot(&res->epsilon[1], NULL, res->niter-1, "check_smc_epsilon.pdf", 
         "plot(d[,1], xlab=\"time index\", ylab=\"epsilon\", main=NA, type=\"l\")");
    plot(&res->acceptance_rate[1], NULL, res->niter-1, "check_smc_accept.pdf", 
         "plot(d[,1], xlab=\"time index\", ylab=\"acceptance rate\", main=NA, type=\"l\")");
    smc_result_free(res);
}
END_TEST

START_TEST (test_smc_toy)
{
    double y = 0;
    int i, fd;

    smc_config config = {
        .nparam = 1,
        .nparticle = 10000,
        .nsample = 2,
        .ess_tolerance = 5000,
        .final_epsilon = 0.01,
        .quality = 0.95,
        .step_tolerance = 1e-9,
        .dataset_size = sizeof(double),
        .feedback_size = sizeof(double)
    };

    smc_result *res = abc_smc(config, toy_functions, 0, (void *) &y);

    for (i = 0; i < config.nparticle; ++i) {
        ck_assert(res->theta[i] > -10 && res->theta[i] < 10);
    }

    plot(res->theta, NULL, config.nparticle, "check_smc_hist.pdf",
         "hist(d[,1], xlim=c(-3, 3), xlab=\"theta\", ylab=\"density\", main=NA, breaks=200)");
    smc_result_free(res);
}
END_TEST

Suite *smc_suite(void)
{
    Suite *s;
    TCase *tc_io, *tc_smc;

    s = suite_create("smc");

    tc_smc = tcase_create("Core");
    tcase_add_test(tc_smc, test_smc_toy_steps);
    tcase_add_test(tc_smc, test_smc_toy);
    tcase_set_timeout(tc_smc, 30);
    suite_add_tcase(s, tc_smc);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = smc_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}