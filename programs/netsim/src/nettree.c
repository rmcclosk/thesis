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
    {"help", no_argument, 0, 'h'},
    {"sim-time", required_argument, 0, 's'},
    {"sim-nodes", required_argument, 0, 'n'},
    {"tree-height", required_argument, 0, 'r'},
    {"tree-tips", required_argument, 0, 't'},
    {"extant-only", no_argument, 0, 'e'},
    {"seed", required_argument, 0, 'd'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: nettree [options] [graph]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -s, --sim-time            amount of epidemic time to simulate (default infinity)\n");
    fprintf(stderr, "  -n, --sim-nodes           stop simulation after this many nodes are infected (default infinity)\n");
    fprintf(stderr, "  -r, --tree-height         cut the tree at this height (default infinity)\n");
    fprintf(stderr, "  -t, --tree-tips           sample this many tips from the tree (default infinity)\n");
    fprintf(stderr, "  -e, --extant-only         sample only extant tips (default false)\n");
    fprintf(stderr, "  -d, --seed                random seed\n");
}

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
        c = getopt_long(argc, argv, "hs:n:r:t:ed:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 's':
                opts.sim_time = atof(optarg);
                break;
            case 'n':
                opts.sim_nodes = atoi(optarg);
                break;
            case 'r':
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
            case '?':
            case 0:
                break;
            default:
                usage();
                exit(EXIT_FAILURE);
        }
    }

    // TODO: safety
    if (optind < argc)
    {
        opts.net_file = fopen(argv[optind++], "r");
        if (optind < argc)
        {
            opts.tree_file = fopen(argv[optind++], "w");
        }
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
    if (!igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_VERTEX, "id"))
    {
        fprintf(stderr, "Vertices must have the 'id' attribute\n");
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
    int i, numeric_ids = 0;
    char buf[128];
    struct nettree_options opts = get_options(argc, argv);
    gsl_rng *rng = set_seed(opts.seed < 0 ? time(NULL) : opts.seed);
    igraph_t net, tree;

    igraph_strvector_t gnames, vnames, enames;
    igraph_vector_t gtypes, vtypes, etypes;

    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);
    igraph_vector_init(&gtypes, 0);
    igraph_vector_init(&vtypes, 0);
    igraph_vector_init(&etypes, 0);

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_read_graph_gml(&net, opts.net_file);
    igraph_cattribute_list(&net, &gnames, &gtypes, &vnames, &vtypes, &enames, &etypes);

    for (i = 0; i < igraph_strvector_size(&vnames); ++i) {
        if (strcmp(STR(vnames, i), "id") == 0 && VECTOR(vtypes)[i] == 1) {
            numeric_ids = 1;
            break;
        }
    }

    if (net_ok(&net))
    {
        simulate_phylogeny(&tree, &net, rng, opts.sim_time, opts.sim_nodes, numeric_ids);
        fprintf(stderr, "Simulated a tree of height %.2f with %d tips\n", 
                height(&tree), (igraph_vcount(&tree) + 1)/2);
        cut_at_time(&tree, opts.tree_height, opts.extant_only);
        subsample_tips(&tree, opts.ntip, rng);
        ladderize(&tree);
        write_tree_newick(&tree, opts.tree_file);
        igraph_destroy(&tree);
    }

    if (opts.net_file != stdin)
        fclose(opts.net_file);
    if (opts.tree_file != stdout)
        fclose(opts.tree_file);
    gsl_rng_free(rng);
    igraph_destroy(&net);
    igraph_strvector_destroy(&gnames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&enames);
    igraph_vector_destroy(&gtypes);
    igraph_vector_destroy(&vtypes);
    igraph_vector_destroy(&etypes);
    return EXIT_SUCCESS;
}
