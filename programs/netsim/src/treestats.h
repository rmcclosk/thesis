/** \file treestats.h
 * \brief Functions for calculating metrics and summary statistics of trees.
 */

#ifndef TREESTATS_H
#define TREESTATS_H

/** Calculate the tree kernel.
 *
 * Uses the fast algorithm from \cite moschitti2006making.
 * See \cite poon2015phylodynamic for the meanings of the parameters lambda and
 * sigma. The dot product of coalescent times is currently not implemented, so
 * the coal parameter is ignored.
 *
 * \param[in] t1,t2 trees to compare
 * \param[in] lambda decay factor in [0, 1] penalizing large matches
 * \param[in] sigma variance of Gaussian radial basis function of branch lengths
 * \param[in] coal if non-zero, multiply the kernel by a dot product of
 * coalescence times
 */
double kernel(const igraph_t *t1, const igraph_t *t2, double lambda, double sigma, int coal);

#endif
