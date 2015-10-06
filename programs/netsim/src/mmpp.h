/** \file mmpp.h
 *
 * Functions relating to fitting an MMPP to a tree.
 */

#ifndef MMPP_H
#define MMPP_H

#include <igraph/igraph.h>

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
 */
double likelihood(const igraph_t *tree, int nrates, const double *theta);

#endif
