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
#define MEAN_DEGREE 8
#define PA_POWER_MIN 0.0
#define PA_POWER_MAX 1.0

struct netabc_options {
    FILE *tree_file;
    int nthread;
    int seed;
    int nparticle;
    int nsample;
    int use_nltt;
    double decay_factor;
    double rbf_variance;
    double quality;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"num-threads", required_argument, 0, 't'},
    {"seed", required_argument, 0, 's'},
    {"decay-factor", required_argument, 0, 'l'},
    {"rbf-variance", required_argument, 0, 'g'},
    {"num-particles", required_argument, 0, 'n'},
    {"num-samples", required_argument, 0, 'p'},
    {"quality", required_argument, 0, 'q'},
    {"nltt", no_argument, 0, 'c'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: netabc [options] [tree]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -t, --num-threads         number of threads\n");
    fprintf(stderr, "  -s, --seed                random seed\n");
    fprintf(stderr, "  -l, --decay-factor        decay factor for tree kernel\n");
    fprintf(stderr, "  -g, --rbf-variance        variance for tree kernel radial basis function\n");
    fprintf(stderr, "  -n, --num-particles       number of particles for SMC\n");
    fprintf(stderr, "  -p, --num-samples         number of sampled datasets per particle\n");
    fprintf(stderr, "  -q, --quality             tradeoff between speed and accuracy (0.9=fast, 0.99=accurate)\n");
    fprintf(stderr, "  -c, --nltt                multiply the kernel by the nLTT statistic\n");
}

struct netabc_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct netabc_options opts = {
        .tree_file = stdin,
        .nthread = 1,
        .seed = -1,
        .decay_factor = 0.2,
        .rbf_variance = 2,
        .nparticle = 1000,
        .nsample = 5,
        .quality = 0.95,
        .use_nltt = 0
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "hcg:l:n:p:q:s:t:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'c':
                opts.use_nltt = 1;
                break;
            case 'g':
                opts.rbf_variance = atof(optarg);
                break;
            case 'l':
                opts.decay_factor = atof(optarg);
                break;
            case 'n':
                opts.nparticle = atoi(optarg);
                break;
            case 'p':
                opts.nsample = atoi(optarg);
                break;
            case 'q':
                opts.quality = atof(optarg);
                break;
            case 's':
                opts.seed = atoi(optarg);
                break;
            case 't':
                opts.nthread = atoi(optarg);
                break;
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

void ba_propose(gsl_rng *rng, double *theta, const void *params)
{
    double var = * (double *) params;
    *theta += gsl_ran_gaussian(rng, sqrt(2*var));
}

double ba_proposal_density(const double *from, const double *to, const void *params)
{
    double var = * (double *) params;
    if (*to < PA_POWER_MIN || *to > PA_POWER_MAX)
        return 0.;
    return gsl_ran_gaussian_pdf(*to - *from, sqrt(2*var));
}

void ba_sample_dataset(gsl_rng *rng, const double *theta, const void *arg, void *X)
{
    int i;
    igraph_t net;
    igraph_t *tree = (igraph_t *) X;
    igraph_vector_t v;
    igraph_rng_t igraph_rng;
    unsigned long int igraph_seed = gsl_rng_get(rng);
    int ntip = (int) ((double *) arg)[0];
    double decay_factor = ((double *) arg)[1];
    double rbf_variance = ((double *) arg)[2];

    // because of thread-local storage this might work?
    igraph_rng_init(&igraph_rng, &igraph_rngtype_mt19937);
    igraph_rng_seed(&igraph_rng, igraph_seed);
    igraph_rng_set_default(&igraph_rng);

    igraph_vector_init(&v, NNODE);
    igraph_barabasi_game(&net, NNODE, *theta, MEAN_DEGREE/2, NULL, 0,
            1, 0, IGRAPH_BARABASI_PSUMTREE, NULL);
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
    subsample_tips(tree, ntip, rng);
    ladderize(tree);
    scale_branches(tree, MEAN);
    SETGAN(tree, "kernel", kernel(tree, tree, decay_factor, rbf_variance, 1));

    igraph_rng_set_default(igraph_rng_default());
    igraph_rng_destroy(&igraph_rng);
    igraph_vector_destroy(&v);
    igraph_destroy(&net);
}

double ba_distance(const void *x, const void *y, const void *arg)
{
    igraph_t *gx = (igraph_t *) x;
    igraph_t *gy = (igraph_t *) y;
    double decay_factor = ((double *) arg)[0];                                   
    double rbf_variance = ((double *) arg)[1]; 
    double nx = NTIP(gx), ny = NTIP(gy);
    double kx = GAN(gx, "kernel");
    double ky = GAN(gy, "kernel");
    double d[14];
    double w[14] = {9.33802949226047e-01, 2.18403741961014e+05,
    1.13819995020497e+01, 7.13025882178519e-02, 5.477396722075e-01,
    6.00418896265072e-02, 4.41223417569416e+05, 6.51247258113592e+00,
    1.24508698036331e+01, 9.17636517221493e-01, 6.89583117577029e-04,
    8.78845237640112e-02, 1, 1};
    double dist = 0;
    int i;

    /*
    d[0] = ((double) colless(gx) - COLLESS_YULE(nx)) / nx -
           ((double) colless(gy) - COLLESS_YULE(ny)) / ny;
    d[1] = (double) il_nodes(gx) / nx - (double) il_nodes(gy) / ny;
    d[2] = (double) ladder_length(gx) / nx - (double) ladder_length(gy) / ny;
    d[3] = (double) width(gx) / (double) ladder_length(gx) - 
           (double) width(gy) / (double) ladder_length(gy);
    d[4] = (double) max_delta_width(gx) / nx - (double) max_delta_width(gy) / ny;
    d[5] = (double) width(gx) / nx - (double) width(gy) / ny;
    d[6] = (double) cherries(gx) / nx - (double) cherries(gy) / ny;
    d[7] = prop_unbalanced(gx) - prop_unbalanced(gy);
    d[8] = avg_unbalance(gx) - avg_unbalance(gy);
    d[9] = (double) (sackin(gx, 0) - SACKIN_YULE(nx)) / nx -
           (double) (sackin(gy, 0) - SACKIN_YULE(ny)) / ny;
    d[10] = (double) (cophenetic(gx, 0) - COPHENETIC_YULE(nx)) / nx -
            (double) (cophenetic(gy, 0) - COPHENETIC_YULE(ny)) / ny;
    d[11] = pybus_gamma(gx) - pybus_gamma(gy);
    */
    d[12] = 1.0 - kernel(gx, gy, decay_factor, rbf_variance, 1) / sqrt(kx) / sqrt(ky);
    //d[13] = 1.0 - nLTT(gx, gy);

    /*
    for (i = 0; i < 14; ++i) {
        dist += pow(w[i] * d[i], 2);
    }
    */
    //dist = pow(d[2], 2) + pow(d[8], 2);
    //return sqrt(dist);
    return d[12];
}

void ba_feedback(const double *theta, int nparticle, void *params)
{
    double var = gsl_stats_variance(theta, 1, nparticle);
    memcpy(params, &var, sizeof(double));
}

void ba_destroy_dataset(void *z)
{
    igraph_destroy((igraph_t *) z);
}

smc_functions ba_functions = {
    .propose = ba_propose,
    .proposal_density = ba_proposal_density,
    .sample_dataset = ba_sample_dataset,
    .distance = ba_distance,
    .feedback = ba_feedback,
    .destroy_dataset = ba_destroy_dataset
};

smc_config ba_config = {
    .nparam = 1,
    .final_epsilon = 0.005,
    .step_tolerance = 1e-5,
    .dataset_size = sizeof(igraph_t),
    .feedback_size = sizeof(double)
};

int main (int argc, char **argv)
{
    int i, ntip;
    struct netabc_options opts = get_options(argc, argv);
    igraph_t *tree;
    smc_result *ba_result;
    double distance_arg[3];
    double sample_dataset_arg[3];

#if IGRAPH_THREAD_SAFE == 0
    if (opts.nthread > 1)
    {
        fprintf(stderr, "Warning: igraph is not thread-safe\n");
        fprintf(stderr, "Disabling multithreading\n");
        fprintf(stderr, "To fix this, install a recent igraph compiled with thread-local storage\n");
        opts.nthread = 1;
    }
#endif

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    tree = parse_newick(opts.tree_file);
    ladderize(tree);
    scale_branches(tree, MEAN);
    SETGAN(tree, "kernel", 
           kernel(tree, tree, opts.decay_factor, opts.rbf_variance, 1));

    ntip = (igraph_vcount(tree) + 1) / 2;
    
    // set up SMC configuration
    ba_config.nparticle = opts.nparticle;
    ba_config.ess_tolerance = opts.nparticle / 2;
    ba_config.nsample = opts.nsample;
    ba_config.quality = opts.quality;

    ba_config.priors = malloc(sizeof(smc_distribution));
    ba_config.priors[0] = UNIFORM;
    ba_config.prior_params = malloc(2 * sizeof(double));
    ba_config.prior_params[0] = PA_POWER_MIN;
    ba_config.prior_params[1] = PA_POWER_MAX;

    sample_dataset_arg[0] = ntip;
    sample_dataset_arg[1] = opts.decay_factor;
    sample_dataset_arg[2] = opts.rbf_variance;
    ba_config.sample_dataset_arg = &sample_dataset_arg;

    distance_arg[0] = opts.decay_factor;
    distance_arg[1] = opts.rbf_variance;
    distance_arg[2] = opts.use_nltt;
    ba_config.distance_arg = &distance_arg;

    ba_result = abc_smc(ba_config, ba_functions, opts.seed, opts.nthread, tree);
    for (i = 0; i < ba_config.nparticle; ++i)
        printf("%f\n", ba_result->theta[i]);

    igraph_destroy(tree);
    smc_result_free(ba_result);
    if (opts.tree_file != stdin) {
        fclose(opts.tree_file);
    }
    free(ba_config.priors);
    free(ba_config.prior_params);
    return EXIT_SUCCESS;
}
