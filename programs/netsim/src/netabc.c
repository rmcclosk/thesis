#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <igraph/igraph.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>

#include "tree.h"
#include "treestats.h"
#include "simulate.h"
#include "smc.h"
#include "util.h"

#define NNODE 5000
#define NSIMNODE 1000
#define NTIP 200
#define MEAN_DEGREE 8
#define PA_POWER_MIN 0.1
#define PA_POWER_MAX 1.0
#define DECAY_FACTOR 0.4
#define RBF_VARIANCE 2

struct netabc_options {
    FILE *tree_file;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: netabc [options] [tree]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
}

struct netabc_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct netabc_options opts = {
        .tree_file = stdin
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "hl:g:s:c:ndb:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case '?':
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }

    // TODO: safety
    if (optind < argc) {
        opts.tree_file = fopen(argv[optind++], "r");
    }

    return opts;
}

void ba_sample_from_prior(gsl_rng *rng, double *theta)
{
    *theta = gsl_ran_flat(rng, 0.1, 1.0);
}

double ba_prior_density(const double *theta)
{
    if (*theta < PA_POWER_MIN || *theta > PA_POWER_MAX) {
        return 0;
    }
    return 1. / (PA_POWER_MAX - PA_POWER_MIN);
}

void ba_propose(gsl_rng *rng, double *theta, const void *params)
{
    double var = * (double *) params;
    *theta += gsl_ran_gaussian(rng, sqrt(2*var));
}

double ba_proposal_density(const double *from, const double *to, const void *params)
{
    double var = * (double *) params;
    return gsl_ran_gaussian_pdf(*to - *from, sqrt(2*var));
}

// TODO: I think this leaks memory
void ba_sample_dataset(gsl_rng *rng, const double *theta, void *X)
{
    int i;
    igraph_t net;
    igraph_t *tree = (igraph_t *) X;
    igraph_vector_t v;

    igraph_vector_init(&v, NNODE);
    igraph_barabasi_game(&net, NNODE, *theta, MEAN_DEGREE/2, NULL, 0, 1, 0, 
                         IGRAPH_BARABASI_PSUMTREE, NULL);
    igraph_to_directed(&net, IGRAPH_TO_DIRECTED_MUTUAL);

    igraph_vector_fill(&v, 0);
    SETVANV(&net, "remove", &v);
    for (i = 0; i < NNODE; ++i) {
        VECTOR(v)[i] = i;
    }
    SETVANV(&net, "id", &v);

    igraph_vector_resize(&v, igraph_ecount(&net));
    igraph_vector_fill(&v, 1);
    SETEANV(&net, "transmit", &v);

    simulate_phylogeny(tree, &net, rng, INFINITY, NSIMNODE, 1);
    subsample_tips(tree, NTIP, rng);
    ladderize(tree);
    scale_branches(tree, MEAN);
    SETGAN(tree, "kernel", kernel(tree, tree, DECAY_FACTOR, RBF_VARIANCE, 1, INFINITY));

    igraph_vector_destroy(&v);
    igraph_destroy(&net);
}

double ba_distance(const void *x, const void *y)
{
    igraph_t *gx = (igraph_t *) x;
    igraph_t *gy = (igraph_t *) y;
    double kx = GAN(gx, "kernel");
    double ky = GAN(gy, "kernel");
    double kxy = kernel(gx, gy, DECAY_FACTOR, RBF_VARIANCE, 1, INFINITY);
    return 1 - kxy / sqrt(kx * ky);
}

void ba_feedback(const double *theta, int nparticle, void *params)
{
    double var = gsl_stats_variance(theta, 1, nparticle);
    memcpy(params, &var, sizeof(double));
}

smc_functions ba_functions = {
    .sample_from_prior = ba_sample_from_prior,
    .prior_density = ba_prior_density,
    .propose = ba_propose,
    .proposal_density = ba_proposal_density,
    .sample_dataset = ba_sample_dataset,
    .distance = ba_distance,
    .feedback = ba_feedback
};

smc_config ba_config = {
    .nparam = 1,
    .nparticle = 1000,
    .nsample = 1,
    .ess_tolerance = 500,
    .final_epsilon = 0.01,
    .quality = 0.9,
    .step_tolerance = 1e-5,
    .dataset_size = sizeof(igraph_t),
    .feedback_size = sizeof(double)
};

int main (int argc, char **argv)
{
    int i;
    struct netabc_options opts = get_options(argc, argv);
    igraph_t *tree;
    smc_result *ba_result;

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    tree = parse_newick(opts.tree_file);
    ladderize(tree);
    scale_branches(tree, MEAN);
    SETGAN(tree, "kernel", 
           kernel(tree, tree, DECAY_FACTOR, RBF_VARIANCE, 1, INFINITY));
    
    ba_result = abc_smc(ba_config, ba_functions, 0, tree);
    for (i = 0; i < ba_config.nparticle; ++i)
        printf("%f\n", ba_result->theta[i]);

    igraph_destroy(tree);
    smc_result_free(ba_result);
    if (opts.tree_file != stdin) {
        fclose(opts.tree_file);
    }
    return EXIT_SUCCESS;
}
