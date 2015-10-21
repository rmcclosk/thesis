/** \file treestats.h
 * \brief Functions for calculating metrics and summary statistics of trees.
 */

#ifndef TREESTATS_H
#define TREESTATS_H

/** Calculate the tree kernel.
 *
 * Uses the fast algorithm from \cite moschitti2006making.
 * See \cite poon2015phylodynamic for the meanings of the parameters decay_factor and
 * rbf_variance.
 *
 * \param[in] t1,t2 trees to compare
 * \param[in] decay_factor decay factor in [0, 1] penalizing large matches
 * \param[in] rbf_variance variance of Gaussian radial basis function of branch lengths
 * \param[in] sst_control between 0 and 1, where 0 is a pure subtree kernel and
 * 1 is a pure subset tree kernel (see \cite moschitti2006making)
 * \param[in] coal if 1, multiply by the nLTT (Janzen et al. 2015)
 */
double kernel(const igraph_t *t1, const igraph_t *t2, double decay_factor,
        double rbf_variance, double sst_control, int coal);

#endif
