#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include "tree.h"
#include "treestats.h"

struct treekernel_options {
    double decay_factor;
    double gauss_factor;
    double sst_control;
    double coal_power;
    int normalize;
    int ladderize;
    scaling scale_branches;
    FILE *tree1_file;
    FILE *tree2_file;
};

struct option long_options[] =
{
    {"decay-factor", required_argument, 0, 'l'},
    {"gauss-factor", required_argument, 0, 'g'},
    {"sst-control", required_argument, 0, 's'},
    {"coal-power", required_argument, 0, 'c'},
    {"normalize", no_argument, 0, 'n'},
    {"ladderize", no_argument, 0, 'd'},
    {"scale-branches", required_argument, 0, 'b'},
    {0, 0, 0, 0}
};

struct treekernel_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct treekernel_options opts = {
        .decay_factor = 0.2,
        .gauss_factor = 2,
        .sst_control = 1.0,
        .coal_power = 0.0,
        .normalize = 0,
        .ladderize = 0,
        .scale_branches = NONE,
        .tree1_file = stdin,
        .tree2_file = stdin
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "l:g:s:c:ndb:", long_options, &i);

        switch (c)
        {
            case 'l':
                opts.decay_factor = atof(optarg);
                break;
            case 'g':
                opts.gauss_factor = atof(optarg);
                break;
            case 's':
                opts.sst_control = atof(optarg);
                break;
            case 'c':
                opts.coal_power = atof(optarg);
                break;
            case 'n':
                opts.normalize = 1;
                break;
            case 'd':
                opts.ladderize = 1;
                break;
            case 'b':
                if (strcmp(optarg, "mean") == 0) {
                    opts.scale_branches = MEAN;
                }
                else if (strcmp(optarg, "median") == 0) {
                    opts.scale_branches = MEDIAN;
                }
            default:
                break;
        }
    }

    // TODO: safety
    if (optind < argc) {
        opts.tree1_file = fopen(argv[optind++], "r");
    }
    if (optind < argc) {
        opts.tree2_file = fopen(argv[optind++], "r");
    }

    return opts;
}

int main (int argc, char **argv)
{
    struct treekernel_options opts = get_options(argc, argv);
    double kdenom = 1;

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_t *t1 = parse_newick(opts.tree1_file);
    igraph_t *t2 = parse_newick(opts.tree2_file);

    if (opts.ladderize) {
        ladderize(t1);
        ladderize(t2);
    }

    scale_branches(t1, opts.scale_branches);
    scale_branches(t2, opts.scale_branches);

    if (opts.normalize) {
        kdenom = sqrt(kernel(t1, t1, opts.decay_factor, opts.gauss_factor, opts.sst_control, opts.coal_power)) *
                 sqrt(kernel(t2, t2, opts.decay_factor, opts.gauss_factor, opts.sst_control, opts.coal_power));
    }

    printf("%f\n", kernel(t1, t2, opts.decay_factor, opts.gauss_factor,
                opts.sst_control, opts.coal_power) / kdenom);

    igraph_destroy(t1);
    igraph_destroy(t2);

    if (opts.tree1_file != stdin) {
        fclose(opts.tree1_file);
    }
    if (opts.tree2_file != stdin) {
        fclose(opts.tree2_file);
    }
    return EXIT_SUCCESS;
}
