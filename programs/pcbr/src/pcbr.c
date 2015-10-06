#define NDEBUG

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <pll/pll.h>
#include <pthread.h>
#include <signal.h>
#include "likelihood.h"
#include "util.h"
#include "tree.h"
#include "thread.h"
#include "../c-cmaes/cmaes_interface.h"
#include "../c-cmaes/boundary_transformation.h"
#include "../C-Thread-Pool/thpool.h"

#define MAX_NRATES 6
#define NUM_THREADS 1
#define CMAES_POP_SIZE 100

double init_sd[MAX_NRATES*MAX_NRATES+1];

typedef struct {
    int nrate;
    int nthread;
    use_tips ut;
    trans_per_branch tpb;
    int tip_pdf;
    int cheating;
    FILE *newick_outfile;
    FILE *annot_outfile;
    FILE *cluster_outfile;
} command_args;

command_args parse_args(int argc, char *argv[]);
double fit(pllNewickTree *tree, double *result, int nrates, threadpool thpool, int nthread, command_args *cmd_args);
void print_rates(double *result, int nrates, double branch_scale);

void usage(void)
#if HAVE_FUNC_ATTRIBUTE_NORETURN
__attribute ((noreturn))
#endif
;

void print_rates(double *result, int nrates, double branch_scale)
{
    int i, j;
    double sorted_rates[MAX_NRATES];
    int rate_order[MAX_NRATES];

    order(result, sorted_rates, rate_order, nrates);

    fprintf(stderr, "rates: ");
    for (i = 0; i < nrates; ++i)
	    fprintf(stderr, "%f ", sorted_rates[i] / branch_scale);
    fprintf(stderr, "\n");

    fprintf(stderr, "Q: ");
    for (i = 0; i < nrates; ++i) {
        fprintf(stderr, "[ ");
        for (j = 0; j < nrates; ++j) {
            if (i == j)
                fprintf(stderr, "   *   ");
            else if (i > j)
	            fprintf(stderr, "%f", result[nrates + rate_order[i]*(nrates-1) + rate_order[j]] / branch_scale);
            else
	            fprintf(stderr, "%f", result[nrates + rate_order[i]*(nrates-1) + rate_order[j] - 1] / branch_scale);
        }
        if (i < nrates - 1)
            fprintf(stderr, " ]\n   ");
        else
            fprintf(stderr, " ]\n");
    }
}

void usage(void)
{
    fprintf(stderr, "Usage: %s [OPTIONS] TREE_FILE\n", PACKAGE_NAME);
    fprintf(stderr, "Model-based phylogenetic clustering by branching rates\n");
    fprintf(stderr, "Example: %s -n 3 -a rates.tsv -n out.nwk -c clusters.tsv tree.nwk\n\n", PACKAGE_NAME);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h, --help                Show this message\n");
    fprintf(stderr, "  -r, --rates               Number of branching rates to fit\n");
    fprintf(stderr, "                            If not supplied, decide by likelihood ratio test\n");
    fprintf(stderr, "  -a, --annot-outfile       File to write fitted branching rates for each node\n");
    fprintf(stderr, "                            Default: no output\n");
    fprintf(stderr, "  -f, --tip-pdf             Treat sampling as branching\n");
    fprintf(stderr, "  -x, --cheating            Always assign tips same state as parents\n");
    fprintf(stderr, "  -n, --newick-outfile      File to write annotated Newick tree\n");
    fprintf(stderr, "                            Default: stdout\n");
    fprintf(stderr, "  -u, --use-tips            Use tips for inference (yes/trans/no)\n");
    fprintf(stderr, "  -p, --trans-per-branch    Allowed transitions per branch (one/any)\n");
    fprintf(stderr, "  -c, --cluster-outfile     File to write cluster assignments for each node\n");
    fprintf(stderr, "                            Default: no output\n");
    fprintf(stderr, "  -t, --num-threads         Number of threads to use\n");
    fprintf(stderr, "                            Default: 1\n");
    exit(EXIT_FAILURE);
}

command_args parse_args(int argc, char *argv[])
{
    int i = 0, c = 0;
    command_args args = {
        .nrate = 0, 
        .nthread = 1, 
        .ut = YES,
        .tpb = ANY,
        .tip_pdf = 0,
        .cheating = 0,
        .newick_outfile = stdout, 
        .annot_outfile = NULL, 
        .cluster_outfile = NULL};

	// read arguments
    static struct option long_options[] = {
        {"rates", required_argument, 0, 'r'},
        {"use-tips", required_argument, 0, 'u'},
        {"trans-per-branch", required_argument, 0, 'p'},
        {"tip-pdf", no_argument, 0, 'f'},
        {"cheating", no_argument, 0, 'x'},
        {"annot-outfile", required_argument, 0, 'a'},
        {"newick-outfile", required_argument, 0, 'n'},
        {"cluster-outfile", required_argument, 0, 'c'},
        {"num-threads", required_argument, 0, 't'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while ((c = getopt_long(argc, argv, "a:c:fn:p:r:t:u:hx", long_options, &i)) != -1) {
        switch (c) {
            case 0:
                break;
            case 'r':
                args.nrate = atoi(optarg);
                break;
            case 'a':
                args.annot_outfile = fopen(optarg, "w");
                break;
            case 'f':
                args.tip_pdf = 1;
                break;
            case 'n':
                args.newick_outfile = fopen(optarg, "w");
                break;
            case 'p':
                if (strcmp(optarg, "one"))
                    args.tpb = ONE;
                else if (strcmp(optarg, "any"))
                    args.tpb = ANY;
            case 'c':
                args.cluster_outfile = fopen(optarg, "w");
                break;
            case 't':
                args.nthread = atoi(optarg);
                break;
            case 'u':
                if (strcmp(optarg, "yes") == 0)
                    args.ut = YES;
                else if (strcmp(optarg, "trans") == 0)
                    args.ut = TRANS_ONLY;
                else if (strcmp(optarg, "no") == 0)
                    args.ut = NO;
                break;
            case 'x':
                args.cheating = 1;
                break;
            case 'h':
                usage();
            case '?':
                break;
            default:
                usage();
        }
    }
    if (optind == argc)
        usage();

    return args;
}

double fit(pllNewickTree *tree, double *result, int nrates, threadpool thpool, int nthread, command_args *cmd_args)
{
    int i, j, cur, dimension = nrates * nrates;
    cmaes_t evo;
    cmaes_boundary_transformation_t boundaries;
	pcbr_workspace *w;
	double *const *pop, *xfinal, *arFunvals, *x_in_bounds;
    double loglik;
    thread_data *tdata = create_thread_data(tree, nrates, nthread, CMAES_POP_SIZE, cmd_args->ut, cmd_args->tpb, cmd_args->tip_pdf, cmd_args->cheating);
    double **args = malloc(CMAES_POP_SIZE * sizeof(double*));

    double *lowerBounds = malloc(dimension * sizeof(double));
    double *upperBounds = malloc(dimension * sizeof(double));

	w = pcbr_create(tree, nrates);
    for(i = 0; i < CMAES_POP_SIZE; ++i)
        args[i] = malloc(dimension * sizeof(double));
	guess_parameters(w, args[0]);

    // bounds on branching rates
    for (i = 0; i < nrates; ++i) {
        lowerBounds[i] = -100;
        upperBounds[i] = 10;
    }

    // bounds on transition rates
    for (i = nrates; i < dimension; ++i) {
        lowerBounds[i] = -100;
        upperBounds[i] = 10;
    }

    cmaes_boundary_transformation_init(&boundaries, lowerBounds, upperBounds, dimension);
    x_in_bounds = cmaes_NewDouble(dimension);

    fprintf(stderr, "fitting %d state model\n", nrates);

	fprintf(stderr, "Initial guess:\n");
    print_rates(args[0], nrates, w->branch_scale);

    double_log(args[0], args[0], dimension);

	// initialize CMA-ES 
	arFunvals = cmaes_init(&evo, dimension, args[0], init_sd, 0, CMAES_POP_SIZE, NULL);
    /*
	evo.sp.stopTolFun = 1e-6;
	evo.sp.stopTolFunHist = 1e-7;
	evo.sp.stopTolX = 1e-4;
    */

	// optimize likelihood function
	while (!cmaes_TestForTermination(&evo)) {

		pop = cmaes_SamplePopulation(&evo);
		for (i = 0; i < CMAES_POP_SIZE; ++i) {
            cmaes_boundary_transformation(&boundaries, pop[i], x_in_bounds, dimension);
			double_exp(x_in_bounds, args[i], dimension);
        }
        set_args(tdata, args, nthread, dimension);

		for (i = 0; i < nthread; ++i)
            thpool_add_work(thpool, do_likelihood, &tdata[i]);
        thpool_wait(thpool);

        cur = 0;
        for (i = 0; i < nthread; ++i) {
            for (j = 0; j < tdata[i].neval; ++j)
                arFunvals[cur++] = tdata[i].arFunvals[j];
        }

		cmaes_UpdateDistribution(&evo, arFunvals);
        fprintf(stderr, ".");
	}
    fprintf(stderr, "\n%s\n",  cmaes_TestForTermination(&evo));
    
	xfinal = cmaes_GetNew(&evo, "xbest");
    cmaes_boundary_transformation(&boundaries, xfinal, x_in_bounds, dimension);
	double_exp(x_in_bounds, result, dimension);
    free(xfinal);

	fprintf(stderr, "Best solution:\n");
    print_rates(result, nrates, w->branch_scale);

	cmaes_exit(&evo);
    loglik = likelihood(w, result, 0, cmd_args->ut, cmd_args->tpb, cmd_args->tip_pdf, cmd_args->cheating);

	pcbr_free(w);
    destroy_thread_data(tdata, nthread);
    for (i = 0; i < CMAES_POP_SIZE; ++i)
        free(args[i]);
    free(args);
    free(x_in_bounds);

    fprintf(stderr, "log10 likelihood of the %d state model is %f\n", nrates, loglik);
    return loglik;
}

int main(int argc, char *argv[])
{
    int i = 0;
	pllNewickTree *tree;
	pllStack *stack;
    double *cur_args, *prev_args;
    double loglik, prev_loglik;//, pvalue;
    pcbr_workspace *w;
    char *newick_out;
    threadpool thpool;
    double cur_bic, prev_bic, pvalue;

	// read arguments
    command_args args = parse_args(argc, argv);

    // start thread pool
    thpool = thpool_init(args.nthread);

    // read tree
	tree = pllNewickParseFile(argv[optind]);
    newick_out = malloc(tree->nodes * 100);

    // set initial standard deviations (could be better)
    for (i = 0; i < MAX_NRATES*MAX_NRATES+1; ++i)
        init_sd[i] = 1;

    cur_args = malloc((MAX_NRATES*MAX_NRATES+1)*sizeof(double));
    prev_args = malloc((MAX_NRATES*MAX_NRATES+1)*sizeof(double));

    if (args.nrate == 0) {
        prev_loglik = fit(tree, prev_args, 1, thpool, args.nthread, &args);
        if (prev_loglik != prev_loglik) {
            fprintf(stderr, "failed to fit %d state model", 1);
            exit(EXIT_FAILURE);
        }

        prev_bic = bic(prev_loglik, 1, tree->nodes-1);
        fprintf(stderr, "BIC for 2 state model is %f\n", prev_bic);
        for (args.nrate = 2; args.nrate <= MAX_NRATES; ++args.nrate) {
            loglik = fit(tree, cur_args, args.nrate, thpool, args.nthread, &args);
            if (loglik != loglik) {
                fprintf(stderr, "failed to fit %d state model", args.nrate);
                break;
            }

            cur_bic = bic(loglik, args.nrate * args.nrate, tree->nodes-1);
            fprintf(stderr, "BIC for %d state model is %f\n", args.nrate, cur_bic);
            pvalue = lrt(prev_loglik, loglik, (args.nrate-1)*(args.nrate-1), 
                         args.nrate*args.nrate);
            fprintf(stderr, "P-value for %d state model is %f\n", args.nrate, pvalue);
            //if (pvalue > 0.05) {
            if (cur_bic > prev_bic) {
                fprintf(stderr, "%d state model is not supported\n", args.nrate);
                break;
            }
            fprintf(stderr, "%d state model is supported\n", args.nrate);
            prev_loglik = loglik;
            prev_bic = cur_bic;
            memcpy(prev_args, cur_args, ((MAX_NRATES*MAX_NRATES+1)*sizeof(double)));
        }
        --args.nrate;
    } else {
        loglik = fit(tree, prev_args, args.nrate, thpool, args.nthread, &args);
    }

	// do ancestral reconstruction
	w = pcbr_create(tree, args.nrate);
	likelihood(w, prev_args, 1, args.ut, args.tpb, args.tip_pdf, args.cheating);
    get_clusters(w, prev_args);

    // write Newick output
    output_tree(w, prev_args, w->nnode, newick_out);
    fprintf(args.newick_outfile, "%s\n", newick_out);
    if (args.newick_outfile != stdout)
        fclose(args.newick_outfile);

    // write CSV output
    if (args.annot_outfile) {
    	stack = tree->tree;
    	for (i = 0; i < w->nnode; ++i) {
            fprintf(args.annot_outfile, "%s\t%f\n", 
                    ((pllNewickNodeInfo *) stack->item)->name,
                    prev_args[w->A[w->nnode - 1 - i]] / w->branch_scale);
    		stack = stack->next;
    	}
        fclose(args.annot_outfile);
    }

    if (args.cluster_outfile) {
    	stack = tree->tree;
    	for (i = 0; i < w->nnode; ++i) {
            fprintf(args.cluster_outfile, "%s\t%d\n", 
                    ((pllNewickNodeInfo *) stack->item)->name,
                    w->cluster[w->nnode - 1 - i]);
    		stack = stack->next;
    	}
        fclose(args.cluster_outfile);
    }

	// clean up
	pllNewickParseDestroy(&tree);
    pcbr_free(w);
    free(cur_args);
    free(prev_args);
    free(newick_out);
    thpool_destroy(thpool);
	return 0;
}
