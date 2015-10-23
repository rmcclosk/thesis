#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <math.h>

#include "tree.h"
#include "treestats.h"

typedef enum {
    TREESTAT_HEIGHT,
    TREESTAT_NTIP,
    TREESTAT_SACKIN,
    TREESTAT_COLLESS
} tree_statistic;

struct treestat_options {
    tree_statistic stat;
    treeshape_norm norm_type;
    int ladderize;
    int use_branch_lengths;
    scaling scale_branches;
    FILE *tree_file;
};

struct option long_options[] =
{
    {"help", no_argument, 0, 'h'},
    {"statistic", required_argument, 0, 's'},
    {"treeshape-norm", required_argument, 0, 'k'},
    {"ignore-branch-lengths", no_argument, 0, 'i'},
    {"ladderize", no_argument, 0, 'd'},
    {"scale-branches", required_argument, 0, 'b'},
    {0, 0, 0, 0}
};

void usage(void)
{
    fprintf(stderr, "Usage: treestat [options] [tree1] [tree2]\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                display this message\n");
    fprintf(stderr, "  -s, --statistic           statistic to compute (see options below)\n");
    fprintf(stderr, "  -k, --treeshape-norm      null model for index statistics (yule/pda/none)\n");
    fprintf(stderr, "  -d, --ladderize           ladderize the tree before computing the statistic\n");
    fprintf(stderr, "  -b, --scale-branches      type of branch scaling to apply (mean/median/max/none, default none)\n");
    fprintf(stderr, "  -i, --ignore-branch-lengths  treat all branches as if they have unit length\n\n");
    fprintf(stderr, "Available statistics:\n");
    fprintf(stderr, "  height                    tree height\n");
    fprintf(stderr, "  ntip                      number of tips\n");
    fprintf(stderr, "  sackin                    Sackin's index\n");
    fprintf(stderr, "  colless                   Colless' index\n");
}

struct treestat_options get_options(int argc, char **argv)
{
    int i, c = 0;
    struct treestat_options opts = {
        .stat = TREESTAT_NTIP,
        .norm_type = TREESHAPE_NORM_NONE,
        .ladderize = 0,
        .use_branch_lengths = 1,
        .scale_branches = NONE,
        .tree_file = stdin
    };

    while (c != -1)
    {
        c = getopt_long(argc, argv, "hs:k:idb:", long_options, &i);
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                break;
            case 'h':
                usage();
                exit(EXIT_SUCCESS);
            case 's':
                if (strcmp(optarg, "height") == 0) {
                    opts.stat = TREESTAT_HEIGHT;
                }
                else if (strcmp(optarg, "sackin") == 0) {
                    opts.stat = TREESTAT_SACKIN;
                }
                else if (strcmp(optarg, "colless") == 0) {
                    opts.stat = TREESTAT_COLLESS;
                }
                else if (strcmp(optarg, "ntip") == 0) {
                    opts.stat = TREESTAT_NTIP;
                }
                else {
                    fprintf(stderr, "Unrecognized tree statistic \"%s\"\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'k':
                if (strcmp(optarg, "yule") == 0) {
                    opts.norm_type = TREESHAPE_NORM_YULE;
                }
                else if (strcmp(optarg, "pda") == 0) {
                    opts.norm_type = TREESHAPE_NORM_PDA;
                }
                break;
            case 'i':
                opts.use_branch_lengths = 0;
                break;
            case 'b':
                if (strcmp(optarg, "mean") == 0) {
                    opts.scale_branches = MEAN;
                }
                else if (strcmp(optarg, "median") == 0) {
                    opts.scale_branches = MEDIAN;
                }
                else if (strcmp(optarg, "max") == 0) {
                    opts.scale_branches = MAX;
                }
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

int main (int argc, char **argv)
{
    int i;
    double s = 0;
    struct treestat_options opts = get_options(argc, argv);

    igraph_i_set_attribute_table(&igraph_cattribute_table);
    igraph_t *t = parse_newick(opts.tree_file);

    if (opts.ladderize) {
        ladderize(t);
    }

    if (!opts.use_branch_lengths) {
        for (i = 0; i < igraph_ecount(t); ++i) {
            SETEAN(t, "length", i, 1.0);
        }
    }
    else {
        scale_branches(t, opts.scale_branches);
    }

    switch (opts.stat) 
    {
        case TREESTAT_NTIP:
            s = (igraph_vcount(t) + 1) / 2;
            break;
        case TREESTAT_HEIGHT:
            s = height(t);
            break;
        case TREESTAT_SACKIN:
            s = sackin(t, 1, opts.norm_type);
            break;
        case TREESTAT_COLLESS:
            s = colless(t, opts.norm_type);
            break;
        default:
            fprintf(stderr, "Unrecognized tree statistic\n");
            return 1;
    }

    if (s == (int) s) {
        printf("%d\n", (int) s);
    }
    else {
        printf("%f\n", s);
    }

    igraph_destroy(t);

    if (opts.tree_file != stdin) {
        fclose(opts.tree_file);
    }
    return EXIT_SUCCESS;
}
