#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#define NO_CHILD -1

#include <pll/pll.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv2.h>

typedef struct {
    double *branch_lengths;
    double *t;
    int *topology;
	int nnode;
    int nrates;
    double branch_scale;

    // dynamic programming matrices and scratch space
	double *L;
    double *Lz;
	int *C;
	int *A;
	int *scale;
    int *cluster;

    // ordering of branch lengths
    double *sorted_branch_lengths;
    int *brlen_order_inv;

    // parameters
    gsl_matrix *Q;
    double *Qbar;
    double *P;
    double *pi;

    // scratch objects for GSL
    gsl_matrix *Qt;
    gsl_matrix *expQt;
    gsl_vector_complex *eval;
    gsl_matrix_complex *evec;
    gsl_eigen_nonsymmv_workspace *w;
} pcbr_workspace;

typedef enum {ONE, ANY} trans_per_branch;
typedef enum {YES, TRANS_ONLY, NO} use_tips;

/** Set up to calculate the tree likelihood.
 *
 * \param tree the tree whose likelihood we will calculate
 * \param nrates number of branching rates
 * \return a pcbr_workspace structure for use in likelihood()
 * \sa pupo()
 */
pcbr_workspace *pcbr_create(const pllNewickTree * tree, int nrates);

/** Destroy a pcbr_workspace structure.
 * 
 * \param p the structure to destroy
 */
void pcbr_free(pcbr_workspace * p);

/** Create a flat representation of a tree.
 *
 * This creates two outputs. topology is an array of size 2*nnodes, where
 * topology[2*i] is the index of the left child of node i, and topology[2*i+1]
 * is the index of the right child of i. branch_lengths is an array of size
 * nnodes where branch_lengths[i] is the length of the branch ancestral to node
 * i.
 *
 * \param tree the tree to flatten, the "tree" attribute of a pllNewickTree
 * \param nnode number of nodes in the tree
 * \param topology array to store the topology in
 * \param branch_lengths array to store branch lengths in
 * \param next you must pass zero for this parameter (it's recursive)
 * \return the number of nodes processed
 */
int flatten(const pllStack *tree, int nnode, int *topology, double *branch_lengths, int next);

void scale_branches(pcbr_workspace *w);

/** Guess initial parameters for likelihood function optimization.
 *
 * This uses R and mixdist to fit a mixture of exponentials distribution to the
 * internal branch lengths. Note that this is a guess only and can't be taken
 * to represent the "true" mixture distribution, since mixtures of exponentials
 * are not in general identifiable by their parameters.
 *
 * \param tree a pcbr_workspace representing the tree to operate on
 * \param args place to store the estimated lambdas and q's
 * \todo R is reducing the numerical precision of the branch lengths for some
 * reason
 */
void guess_parameters(pcbr_workspace *tree, double *args);

/** Guess rate parameters for an exponential mixture model.
 *
 * This is used internally by guess_parameters() and isn't meant to be called
 * directly.
 *
 * \param x a vector of REALSXP
 * \param nrates the number of rates in the miuture model
 */
void guess_means(const double *x, int nrates, int nbranch, double *guess);

/** Fill in model parameters in a pcbr_workspace.
 *
 * This fills in Q and pi in the pcbr_workspace object.
 *
 * \param tree the pcbr_workspace to fill
 * \param q transition rates in row-major order, excluding diagonals
 */
void fill_parameters(pcbr_workspace *tree, double *q);

/** Calculate the likelihood of a tree.
 *
 * \param tree a pcbr_workspace produced by pcbr_create()
 * \param args the vector of model parameters, first lambdas, then q's
 * \param options options controlling the algorithm
 * \return the likelihood or negative log10 likelihood
 * \sa pcbr_create()
 */
double likelihood(pcbr_workspace *tree, double *args, int reconstruct,
                  use_tips ut, trans_per_branch tpb, int tip_pdf, int cheating);

/** Get the cluster assignment of each node.
 *
 * \param tree a pcbr_workspace which has been passed through likelihood() with
 * the reconstruct flag
 */
void get_clusters(pcbr_workspace *tree, double *rates);

/** Calculate a likelihood ratio test.
 *
 * \param dll log likelihood difference between alternative and null model
 * \param ddf difference in number of parameters
 * \return P-value of the test
 */
double lrt(double log10lik_null, double log10lik_alt, int nparam_null, int nparam_alt);

/** Calculate the Bayesian Information Criterion.
 *
 * \param log10lik log10 likelihood of the fitted model
 * \param nparam number of parameters of the model
 * \param ndata number of data points
 */
double bic(double log10lik, int nparam, int ndata);

#endif
