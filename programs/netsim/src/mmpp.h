/** \file mmpp.h
 *
 * Functions relating to fitting an MMPP to a tree.
 */

#ifndef MMPP_H
#define MMPP_H

#include <igraph/igraph.h>

typedef struct mmpp_workspace mmpp_workspace;

/** Fit an MMPP.
 *
 * \param[in] tree tree to fit MMPP to
 * \param[in] nrates number of states in the Markov chain
 * \param[out] theta fitted parameters will be stored here
 * \param[in] trace if 1, display parameter values as they are tried
 * \param[in] cmaes_settings file containing CMAES settings
 * \param[out] states if non-NULL, perform ancestral reconstruction and put the result here
 * \return 0 if the fit was successful, 1 otherwise
 */
int fit_mmpp(const igraph_t *tree, int nrates, double *theta, int trace,
        const char *cmaes_settings, int *states);

/** Guess initial parameters for an MMPP.
 *
 * To guess branching rates, partition the internal branch lengths into nrates
 * groups by length, take the inverse of the mean of each group. Guess
 * transition rates all equal, such that the expected number of transitions is
 * 1% of the total number of branches (eg. with 500 branches, we expect 5
 * transitions). The parameters are returned in the vector theta, which must
 * have enough space allocated for nrates*nrates doubles. The first nrates
 * elements will be the estimated branching rates, and the last
 * nrates*(nrates-1) will be the estimated transition rates, stored in
 * row-major order with the diagonals ommitted.
 *
 * \param[in] tree tree to estimate parameters for
 * \param[in] nrates number of rates in the MMPP
 * \param[out] theta estimated parameters will be stored here
 */
void guess_parameters(const igraph_t *tree, int nrates, double *theta);

/** Calculate the likelihood of a tree under an MMPP.
 *
 * \param[in] tree the tree to calcluate the likelihood for
 * \param[in] nrates number of rates of the MMPP
 * \param[in] theta MMPP parameters (see guess_parameters)
 * \param[in] w workspace for calculations, made by mmpp_workspace_create
 * \param[in] reconstruct if 1, perform ancestral reconstruction
 */
double likelihood(const igraph_t *tree, int nrates, const double *theta,
        mmpp_workspace *w, int reconstruct);

/** Perform ancestral reconstruction.
 *
 * \param[in] tree tree to reconstruct ancestral states for
 * \param[in] nrates number of states in MMPP
 * \param[in] theta fitted MMPP parameters
 * \param[in] w workspace created by mmpp_workspace_create
 * \param[out] A ancestral states will be stored here
 * \return the likelihood of the assignment of ancestral states to nodes
 */
double reconstruct(const igraph_t *tree, int nrates, const double *theta,
        mmpp_workspace *w, int *A);

/** Create a workspace for MMPP likelihood calculations.
 *
 * This pre-allocates space and initializes structures needed to calculate MMPP
 * likelihoods.
 *
 * \param[in] nrates number of MMPP rates
 * \return a workspace for use in likelihood()
 */
mmpp_workspace *mmpp_workspace_create(const igraph_t *tree, int nrates);

/** Free memory associated with an MMPP workspace.
 *
 * \param[in] w workspace to destory
 */
void mmpp_workspace_free(mmpp_workspace *w);

#endif
