#include <stdio.h>
#include <getopt.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "simulate.h"
#include "util.h"
#include "tree.h"

struct nettree_options {
    double sim_time;
    int sim_nodes;
    double tree_height;
    int extant_only;
    int ntip;
    FILE *net_file;
    FILE *tree_file;
    int seed;
};

struct option long_options[] =
{
    {"sim-time", required_argument, 0, 's'},
    {"sim-nodes", required_argument, 0, 'n'},
    {"tree-height", required_argument, 0, 'h'},
    {"tree-tips", required_argument, 0, 't'},
    {"extant-only", no_argument, 0, 'e'},
    {"seed", required_argument, 0, 'd'},
    {0, 0, 0, 0}
};

struct nettree_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct nettree_options opts = {
        .sim_time = DBL_MAX, 
        .sim_nodes = INT_MAX,
        .tree_height = DBL_MAX,
        .extant_only = 0,
        .ntip = INT_MAX,
        .net_file = stdin,
        .tree_file = stdout,
        .seed = -1
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "s:n:h:t:ed:", long_options, &i);

        switch (c)
        {
            case 's':
                opts.sim_time = atof(optarg);
                break;
            case 'n':
                opts.sim_nodes = atoi(optarg);
                break;
            case 'h':
                opts.tree_height = atof(optarg);
                break;
            case 't':
                opts.ntip = atoi(optarg);
                break;
            case 'e':
                opts.extant_only = 1;
                break;
            case 'd':
                opts.seed = atoi(optarg);
                break;
            default:
                break;
        }
    }

    // TODO: safety
    if (optind < argc)
    {
        opts.net_file = fopen(argv[optind++], "r");
        if (optind < argc)
            opts.tree_file = fopen(argv[optind++], "w");
    }

    return opts;
}

int net_ok(igraph_t *net)
{
    if (!igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_VERTEX, "remove"))
    {
        fprintf(stderr, "Vertices must have the 'remove' attribute\n");
        return 0;
    }
    if (!igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_EDGE, "transmit"))
    {
        fprintf(stderr, "Edges must have the 'transmit' attribute\n");
        return 0;
    }
    return 1;
}

int main (int argc, char **argv)
{
    struct nettree_options opts = get_options(argc, argv);
    gsl_rng *rng = set_seed(opts.seed < 0 ? time(NULL) : opts.seed);
    igraph_t net, *tree;

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_gml(&net, opts.net_file);

    if (net_ok(&net))
    {
        tree = simulate_phylogeny(&net, rng, opts.sim_time, opts.sim_nodes);
        fprintf(stderr, "Simulated a tree of height %.2f with %d tips\n", 
                height(tree), (igraph_vcount(tree) + 1)/2);
        cut_at_time(tree, opts.tree_height, opts.extant_only);
        subsample_tips(tree, opts.ntip, rng);
        ladderize(tree);
        write_tree_newick(tree, opts.tree_file);
        igraph_destroy(tree);
    }

    if (opts.net_file != stdin)
        fclose(opts.net_file);
    if (opts.tree_file != stdout)
        fclose(opts.tree_file);
    gsl_rng_free(rng);
    igraph_destroy(&net);
    return EXIT_SUCCESS;
}
