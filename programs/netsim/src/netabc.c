#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <yaml.h>
#include <igraph/igraph.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>

#include "tree.h"
#include "treestats.h"
#include "simulate.h"
#include "smc.h"
#include "util.h"

#define MAX_PARAMS 7
#define SIMULATE_MAX_TRIES 100

typedef enum {
    PREF_ATTACH,
    RANDOM_GNP,
    SMALL_WORLD,
    PARETO
} net_type;

struct netabc_options {
    FILE *tree_file;
    FILE *yaml_file;
    FILE *trace_file;
    int nthread;
    int seed;
    int nparticle;
    int nsample;
    net_type net;
    double decay_factor;
    double rbf_variance;
    double quality;
    double final_epsilon;
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
    {"trace", required_argument, 0, 'd'},
    {"net-type", required_argument, 0, 'm'},
    {"final-epsilon", required_argument, 0, 'e'},
    {0, 0, 0, 0}
};

typedef enum {
    NNODE = 0,
    NSIMNODE = 1,
    SIM_TIME = 2,
    TRANSMIT_RATE = 3,
    REMOVE_RATE = 4,

    // for random networks
    PR_EDGE = 5,

    // for preferential attachment networks
    EDGES_PER_VERTEX = 5,
    ATTACH_POWER = 6,

    // for smallworld networks
    NBHD_SIZE = 5,
    REWIRE_PROB = 6,

    // for Pareto networks
    PARETO_SHAPE = 5
} net_parameter;

void usage(void)
{
    fprintf(stderr, "Usage: netabc [options] [tree] [yaml_file]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -t, --num-threads         number of threads\n");
    fprintf(stderr, "  -s, --seed                random seed\n");
    fprintf(stderr, "  -l, --decay-factor        decay factor for tree kernel\n");
    fprintf(stderr, "  -g, --rbf-variance        variance for tree kernel radial basis function\n");
    fprintf(stderr, "  -n, --num-particles       number of particles for SMC\n");
    fprintf(stderr, "  -p, --num-samples         number of sampled datasets per particle\n");
    fprintf(stderr, "  -q, --quality             tradeoff between speed and accuracy (0.9=fast, 0.99=accurate)\n");
    fprintf(stderr, "  -d, --trace               write population and weights at each iteration to this file\n");
    fprintf(stderr, "  -m, --net-type            type of network (pa/gnp/sw/pareto)\n");
    fprintf(stderr, "  -e, --final-epsilon       last epsilon for SMC\n");
}

struct netabc_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct netabc_options opts = {
        .tree_file = stdin,
        .yaml_file = NULL,
        .trace_file = NULL,
        .nthread = 1,
        .seed = -1,
        .nparticle = 1000,
        .nsample = 5,
        .net = RANDOM_GNP,
        .decay_factor = 0.2,
        .rbf_variance = 2,
        .quality = 0.95,
        .final_epsilon = 0.01
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "hd:e:g:l:m:n:p:q:s:t:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 'd':
                opts.trace_file = fopen(optarg, "w");
                break;
            case 'e':
                opts.final_epsilon = atof(optarg);
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
            case 'm':
                if (strcmp(optarg, "pa") == 0) {
                    opts.net = PREF_ATTACH;
                }
                else if (strcmp(optarg, "gnp") == 0) {
                    opts.net = RANDOM_GNP;
                }
                else if (strcmp(optarg, "sw") == 0) {
                    opts.net = SMALL_WORLD;
                }
                else if (strcmp(optarg, "pareto") == 0) {
                    opts.net = PARETO;
                }
                else {
                    fprintf(stderr, "Error: unrecognized network type \"%s\"\n", optarg);
                    exit(EXIT_FAILURE);
                }
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
    if (optind < argc) {
        opts.yaml_file = fopen(argv[optind++], "r");
    }

    return opts;
}

void set_parameter_defaults(smc_distribution *priors, double *prior_params, net_type net)
{
    priors[NNODE] = DELTA;
    prior_params[NNODE * MAX_DIST_PARAMS] = 5000;

    priors[NSIMNODE] = DELTA;
    prior_params[NSIMNODE * MAX_DIST_PARAMS] = 1000;

    priors[SIM_TIME] = DELTA;
    prior_params[SIM_TIME * MAX_DIST_PARAMS] = DBL_MAX;

    priors[TRANSMIT_RATE] = DELTA;
    prior_params[TRANSMIT_RATE * MAX_DIST_PARAMS] = 1;

    priors[REMOVE_RATE] = DELTA;
    prior_params[REMOVE_RATE * MAX_DIST_PARAMS] = 0;

    if (net == PREF_ATTACH) {
        priors[EDGES_PER_VERTEX] = DELTA;
        prior_params[EDGES_PER_VERTEX * MAX_DIST_PARAMS] = 2;

        priors[ATTACH_POWER] = UNIFORM;
        prior_params[ATTACH_POWER * MAX_DIST_PARAMS] = 0;
        prior_params[ATTACH_POWER * MAX_DIST_PARAMS + 1] = 1;
    }
    else if (net == RANDOM_GNP) {
        priors[PR_EDGE] = UNIFORM;
        prior_params[PR_EDGE * MAX_DIST_PARAMS] = 0;
        prior_params[PR_EDGE * MAX_DIST_PARAMS + 1] = 0.1;
    }
    else if (net == SMALL_WORLD) {
        priors[NBHD_SIZE] = DELTA;
        prior_params[NBHD_SIZE * MAX_DIST_PARAMS] = 2;

        priors[REWIRE_PROB] = UNIFORM;
        prior_params[REWIRE_PROB * MAX_DIST_PARAMS] = 0;
        prior_params[REWIRE_PROB * MAX_DIST_PARAMS + 1] = 0.1;
    }
    else if (net == PARETO) {
        priors[PARETO_SHAPE] = UNIFORM;
        prior_params[PARETO_SHAPE * MAX_DIST_PARAMS] = 1;
        prior_params[PARETO_SHAPE * MAX_DIST_PARAMS + 1] = 2;
    }
}

void get_parameters(FILE *f, smc_distribution *priors, double *prior_params, net_type net)
{
    int key = 1, seq = 0, seq_cur = 0, param;
    yaml_parser_t parser;
    yaml_event_t event;

    set_parameter_defaults(priors, prior_params, net);
    if (f == NULL) {
        return;
    }

    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, f);
    do
    {
        yaml_parser_parse(&parser, &event);
        switch (event.type)
        {
            case YAML_SEQUENCE_START_EVENT:
                seq = 1;
                seq_cur = 0;
                break;
            case YAML_SEQUENCE_END_EVENT:
                seq = 0;
                break;
            case YAML_SCALAR_EVENT:
                if (key) {
                    if (strcmp(event.data.scalar.value, "nodes") == 0) {
                        param = NNODE;
                    }
                    else if (strcmp(event.data.scalar.value, "sim_nodes") == 0) {
                        param = NSIMNODE;
                    }
                    else if (strcmp(event.data.scalar.value, "sim_time") == 0) {
                        param = SIM_TIME;
                    }
                    else if (strcmp(event.data.scalar.value, "transmit") == 0) {
                        param = TRANSMIT_RATE;
                    }
                    else if (strcmp(event.data.scalar.value, "remove") == 0) {
                        param = REMOVE_RATE;
                    }
                    else if (strcmp(event.data.scalar.value, "pr_edge") == 0) {
                        param = PR_EDGE;
                    }
                    else if (strcmp(event.data.scalar.value, "edges_per_vertex") == 0) {
                        param = EDGES_PER_VERTEX;
                    }
                    else if (strcmp(event.data.scalar.value, "nbhd_size") == 0) {
                        param = NBHD_SIZE;
                    }
                    else if (strcmp(event.data.scalar.value, "rewire_prob") == 0) {
                        param = REWIRE_PROB;
                    }
                    else if (strcmp(event.data.scalar.value, "pareto_shape") == 0) {
                        param = PARETO_SHAPE;
                    }
                    else {
                        fprintf(stderr, "Error: unrecognized parameter \"%s\" in YAML file\n",
                                event.data.scalar.value);
                        exit(EXIT_FAILURE);
                    }
                }
                else if (seq) {
                    if (seq_cur == 0) {
                        if (strcmp(event.data.scalar.value, "uniform") == 0) {
                            priors[param] = UNIFORM;
                        }
                        else if (strcmp(event.data.scalar.value, "gaussian") == 0) {
                            priors[param] = GAUSSIAN;
                        }
                        else if (strcmp(event.data.scalar.value, "delta") == 0) {
                            priors[param] = GAUSSIAN;
                        }
                        else {
                            fprintf(stderr, "Error: unrecognized distribution \"%s\"\n",
                                    event.data.scalar.value);
                            exit(EXIT_FAILURE);
                        }
                    }
                    else if (seq_cur <= MAX_DIST_PARAMS) {
                        prior_params[MAX_DIST_PARAMS * param + seq_cur - 1] = atof(event.data.scalar.value);
                    } 
                    else {
                        fprintf(stderr, "Error: too many distribution parameters\n");
                        exit(EXIT_FAILURE);
                    }
                    ++seq_cur;
                }
                else {
                    priors[param] = DELTA;
                    prior_params[MAX_DIST_PARAMS * param] = atof(event.data.scalar.value);
                }

                if (!seq) {
                    key = !key;
                }
                break;
            default:
                break;
        }
        if (event.type != YAML_STREAM_END_EVENT)
            yaml_event_delete(&event);
    }
    while (event.type != YAML_STREAM_END_EVENT);

    yaml_event_delete(&event);
    yaml_parser_delete(&parser);
}

void propose(gsl_rng *rng, double *theta, const void *params, int nparam)
{
    int i;
    double *var = (double *) params;
    for (i = 0; i < nparam; ++i) {
        theta[i] += gsl_ran_gaussian(rng, sqrt(2*var[i]));
    }
}

void propose_gnp(gsl_rng *rng, double *theta, const void *params)
{
    propose(rng, theta, params, 6);
}

void propose_pa(gsl_rng *rng, double *theta, const void *params)
{
    propose(rng, theta, params, 7);
}

void propose_sw(gsl_rng *rng, double *theta, const void *params) {
    propose(rng, theta, params, 7);
}

void propose_pareto(gsl_rng *rng, double *theta, const void *params)
{
    propose(rng, theta, params, 6);
}

double proposal_density(const double *from, const double *to, const void *params, int nparam)
{
    int i;
    double *var = (double *) params;
    double p = 1;

    for (i = 0; i < nparam; ++i) {
        if (fabs(var[i]) > FLT_EPSILON || fabs(to[i] - from[i]) > FLT_EPSILON) {
            p *= gsl_ran_gaussian_pdf(to[i] - from[i], sqrt(2*var[i]));
        }
    }
    return p;
}

double proposal_density_gnp(const double *from, const double *to, const void *params)
{
    return proposal_density(from, to, params, 6);
}

double proposal_density_pa(const double *from, const double *to, const void *params)
{
    return proposal_density(from, to, params, 7);
}

double proposal_density_sw(const double *from, const double *to, const void *params)
{
    return proposal_density(from, to, params, 7);
}

double proposal_density_pareto(const double *from, const double *to, const void *params)
{
    return proposal_density(from, to, params, 6);
}

struct sample_dataset_arg {
    int ntip;
    double decay_factor;
    double rbf_variance;
    int (*sample_network) (igraph_t *, gsl_rng *, const double *);
};

int sample_network_gnp(igraph_t *net, gsl_rng *rng, const double *theta)
{
    igraph_erdos_renyi_game(net, IGRAPH_ERDOS_RENYI_GNP, theta[NNODE], 
                            theta[PR_EDGE], 0, 0);
    return 0;
}

int sample_network_pa(igraph_t *net, gsl_rng *rng, const double *theta)
{
    igraph_barabasi_game(net, theta[NNODE], theta[ATTACH_POWER], theta[EDGES_PER_VERTEX],
                         NULL, 0, 1, 0, IGRAPH_BARABASI_PSUMTREE, NULL);
    return 0;
}

int sample_network_sw(igraph_t *net, gsl_rng *rng, const double *theta)
{
    igraph_watts_strogatz_game(net, 1, theta[NNODE], theta[NBHD_SIZE], theta[REWIRE_PROB],
                               0, 0);
    return 1;
}

int sample_network_pareto(igraph_t *net, gsl_rng *rng, const double *theta)
{
    int done = 0, i = 0, j, is_graphical;
    igraph_vector_t s;

    igraph_vector_init(&s, theta[NNODE]);

    igraph_set_error_handler(igraph_error_handler_ignore);
    while (!done) {
        if (i == SIMULATE_MAX_TRIES) {
            fprintf(stderr, "Too many tries to generate a network with parameter %f\n",
                            theta[PARETO_SHAPE]);
            return 1;
        }

        for (j = 0; j < theta[NNODE]; ++j) {
            VECTOR(s)[j] = floor(gsl_ran_pareto(rng, theta[PARETO_SHAPE], 1));
        }
        igraph_is_graphical_degree_sequence(&s, NULL, &is_graphical);
        if (!is_graphical) {
            VECTOR(s)[0] += 1;
        }
        igraph_is_graphical_degree_sequence(&s, NULL, &is_graphical);
        if (is_graphical)
        {
            done = !igraph_degree_sequence_game(net, &s, NULL, IGRAPH_DEGSEQ_VL);
        }
        ++i;
    }
    igraph_set_error_handler(igraph_error_handler_abort);
    igraph_vector_destroy(&s);
    return 0;
}

void sample_dataset(gsl_rng *rng, const double *theta, const void *arg, void *X)
{
    int i, failed = 0;
    igraph_t net;
    igraph_t *tree = (igraph_t *) X;
    igraph_vector_t v;
    igraph_rng_t igraph_rng;
    unsigned long int igraph_seed = gsl_rng_get(rng);

    struct sample_dataset_arg *sarg = (struct sample_dataset_arg *) arg;
    int ntip = sarg->ntip;
    double decay_factor = sarg->decay_factor;
    double rbf_variance = sarg->rbf_variance;

    // seems to work because of thread-local storage
    igraph_rng_init(&igraph_rng, &igraph_rngtype_mt19937);
    igraph_rng_seed(&igraph_rng, igraph_seed);
    igraph_rng_set_default(&igraph_rng);

    igraph_vector_init(&v, theta[NNODE]);

    failed = sarg->sample_network(&net, rng, theta);
    if (!failed) {
        igraph_to_directed(&net, IGRAPH_TO_DIRECTED_MUTUAL);
        
        igraph_vector_fill(&v, 0);
        SETVANV(&net, "remove", &v);
        for (i = 0; i < theta[NNODE]; ++i) {
            VECTOR(v)[i] = i;
        }
        SETVANV(&net, "id", &v);
    
        igraph_vector_resize(&v, igraph_ecount(&net));
        igraph_vector_fill(&v, theta[TRANSMIT_RATE]);
        SETEANV(&net, "transmit", &v);
        
        simulate_phylogeny(tree, &net, rng, theta[SIM_TIME], theta[NSIMNODE], 1);
        i = 0;
        while (igraph_vcount(tree) < (ntip - 1) / 2) {
            if (i == 20) {
                fprintf(stderr, "Too many tries to simulate a tree\n");
                failed = 1;
                break;
            }
            igraph_destroy(tree);
            simulate_phylogeny(tree, &net, rng, theta[SIM_TIME], theta[NSIMNODE], 1);
            ++i;
        }
    }

    if (failed) {
        memset(tree, 0, sizeof(igraph_t));
    }
    else {
        subsample_tips(tree, ntip, rng);
        ladderize(tree);
        scale_branches(tree, MEAN);
        SETGAN(tree, "kernel", kernel(tree, tree, decay_factor, rbf_variance, 1));
        igraph_destroy(&net);
    }

    igraph_rng_set_default(igraph_rng_default());
    igraph_rng_destroy(&igraph_rng);
    igraph_vector_destroy(&v);
}

double distance(const void *x, const void *y, const void *arg)
{
    igraph_t *gx = (igraph_t *) x;
    igraph_t *gy = (igraph_t *) y;
    double decay_factor = ((double *) arg)[0];                                   
    double rbf_variance = ((double *) arg)[1]; 
    char *zeroes = calloc(sizeof(igraph_t), 1);
    double kx, ky, dist;

    if (memcmp(x, zeroes, sizeof(igraph_t)) == 0 ||
        memcmp(y, zeroes, sizeof(igraph_t)) == 0) {
        dist = INFINITY;
    }
    else {
        kx = GAN(gx, "kernel");
        ky = GAN(gy, "kernel");
        dist = 1.0 - kernel(gx, gy, decay_factor, rbf_variance, 1) / sqrt(kx) / sqrt(ky);
    }

    free(zeroes);
    return dist;
}

void feedback(const double *theta, int nparticle, void *params, int nparam)
{
    int i;
    double *var = (double *) params;
    for (i = 0; i < nparam; ++i) {
        var[i] = gsl_stats_variance(&theta[i], nparam, nparticle);
    }
}

void feedback_gnp(const double *theta, int nparticle, void *params)
{
    feedback(theta, nparticle, params, 6);
}

void feedback_pa(const double *theta, int nparticle, void *params)
{
    feedback(theta, nparticle, params, 7);
}

void feedback_sw(const double *theta, int nparticle, void *params)
{
    feedback(theta, nparticle, params, 7);
}

void feedback_pareto(const double *theta, int nparticle, void *params)
{
    feedback(theta, nparticle, params, 6);
}

void destroy_dataset(void *z)
{
    igraph_destroy((igraph_t *) z);
}

smc_functions functions = {
    .sample_dataset = sample_dataset,
    .distance = distance,
    .destroy_dataset = destroy_dataset
};

smc_config config = {
    .step_tolerance = 1e-5,
    .dataset_size = sizeof(igraph_t),
};

int main (int argc, char **argv)
{
    int i, j, k;
    struct netabc_options opts = get_options(argc, argv);
    igraph_t *tree;
    smc_result *result;
    double distance_arg[2];
    struct sample_dataset_arg sarg;
    smc_distribution *priors = malloc(MAX_PARAMS * sizeof(smc_distribution));
    double *prior_params = malloc(MAX_PARAMS * MAX_DIST_PARAMS * sizeof(double));

#if IGRAPH_THREAD_SAFE == 0
    if (opts.nthread > 1)
    {
        fprintf(stderr, "Warning: igraph is not thread-safe\n");
        fprintf(stderr, "Disabling multithreading\n");
        fprintf(stderr, "To fix this, install a recent igraph compiled with thread-local storage\n");
        opts.nthread = 1;
    }
#endif

    get_parameters(opts.yaml_file, priors, prior_params, opts.net);
    if (opts.yaml_file != NULL) {
        fclose(opts.yaml_file);
    }

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    tree = parse_newick(opts.tree_file);
    if (opts.tree_file != stdin) {
        fclose(opts.tree_file);
    }

    ladderize(tree);
    scale_branches(tree, MEAN);
    SETGAN(tree, "kernel", 
           kernel(tree, tree, opts.decay_factor, opts.rbf_variance, 1));

    // set up SMC configuration
    config.nparticle = opts.nparticle;
    config.ess_tolerance = opts.nparticle / 2;
    config.nsample = opts.nsample;
    config.quality = opts.quality;
    config.final_epsilon = opts.final_epsilon;

    config.priors = priors;
    config.prior_params = prior_params;

    sarg.ntip = NTIP(tree);
    sarg.decay_factor = opts.decay_factor;
    sarg.rbf_variance = opts.rbf_variance;
    switch (opts.net) {
        case PREF_ATTACH:
            sarg.sample_network = sample_network_pa;
            functions.propose = propose_pa;
            functions.proposal_density = proposal_density_pa;
            functions.feedback = feedback_pa;
            config.nparam = 7;
            config.feedback_size = 7 * sizeof(double);
            break;
        case SMALL_WORLD:
            sarg.sample_network = sample_network_sw;
            functions.propose = propose_sw;
            functions.proposal_density = proposal_density_sw;
            functions.feedback = feedback_sw;
            config.nparam = 7;
            config.feedback_size = 7 * sizeof(double);
            break;
        case PARETO:
            sarg.sample_network = sample_network_pareto;
            functions.propose = propose_pareto;
            functions.proposal_density = proposal_density_pareto;
            functions.feedback = feedback_pareto;
            config.nparam = 6;
            config.feedback_size = 6 * sizeof(double);
            break;
        default:
            sarg.sample_network = sample_network_gnp;
            functions.propose = propose_gnp;
            functions.proposal_density = proposal_density_gnp;
            functions.feedback = feedback_gnp;
            config.nparam = 6;
            config.feedback_size = 6 * sizeof(double);
            break;
    }
    config.sample_dataset_arg = &sarg;

    distance_arg[0] = opts.decay_factor;
    distance_arg[1] = opts.rbf_variance;
    config.distance_arg = &distance_arg;

    if (opts.trace_file != NULL) {
        fprintf(opts.trace_file, "iter\tweight\tnodes\tsim_nodes\tsim_time\ttransmit\tremove\t");
        switch (opts.net) {
            case PREF_ATTACH:
                fprintf(opts.trace_file, "edges_per_vertex\tattach_power\n");
                break;
            case SMALL_WORLD:
                fprintf(opts.trace_file, "nbhd_size\trewire_prob\n");
                break;
            case PARETO:
                fprintf(opts.trace_file, "pareto_shape\n");
                break;
            default:
                fprintf(opts.trace_file, "pr_edge\n");
                break;
        }
    }

    result = abc_smc(config, functions, opts.seed, opts.nthread, tree, opts.trace_file);

    if (opts.trace_file != NULL) {
        fclose(opts.trace_file);
    }

    igraph_destroy(tree);
    smc_result_free(result);
    free(priors);
    free(prior_params);
    return EXIT_SUCCESS;
}
