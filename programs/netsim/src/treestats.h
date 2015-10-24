/** \file treestats.h
 * \brief Functions for calculating metrics and summary statistics of trees.
 */

#ifndef TREESTATS_H
#define TREESTATS_H

/** Calculate the tree kernel.
 *
 * Uses the fast algorithm from \cite moschitti2006making.
 * See \cite poon2015phylodynamic for the meanings of the parameters
 * decay_factor and rbf_variance.
 *
 * \param[in] t1,t2 trees to compare
 * \param[in] decay_factor decay factor in [0, 1] penalizing large matches
 * \param[in] rbf_variance variance of Gaussian radial basis function of branch lengths
 * \param[in] sst_control between 0 and 1, where 0 is a pure subtree kernel and
 * 1 is a pure subset tree kernel (see \cite moschitti2006making)
 */
double kernel(const igraph_t *t1, const igraph_t *t2, double decay_factor,
        double rbf_variance, double sst_control);

/** Compute the normalized lineages-through-time statistic.
 *
 * Slightly modified from the version in \cite janzen2015approximate, using the
 * trapezoid rule instead of the rectangle method.
 *
 * \param[in] t1 first tree to compare
 * \param[in] t2 second tree to compare
 * \return the nLTT statistic
 */
double nLTT(const igraph_t *t1, const igraph_t *t2);

typedef enum {
    TREESHAPE_NORM_NONE,
    TREESHAPE_NORM_YULE,
    TREESHAPE_NORM_PDA
} treeshape_norm;

/** Compute Sackin's index.
 *
 * The norm parameter specifies a null model. The result will be divided by the
 * Sackin's index of the null model (see the sackin function in the apTreeshape
 * R package). It doesn't really make sense to use normalization other than
 * NONE with use_branch_lengths = 1, but it's allowed.
 *
 * \param[in] t tree to compute Sackin's index for
 * \param[in] use_branch_lengths if 0, treat all branches as if they had unit
 * length
 * \param[in] norm null model to use for normalizing the result
 * \return the average path length from tips to the root
 */
double sackin(const igraph_t *t, int use_branch_lengths, treeshape_norm norm);

/** Compute Colless' index.
 *
 * \param[in] t tree t compute Colless' index for
 * \param[in] norm null model to use for normalizing the result
 * \return Colless' index
 */
double colless(const igraph_t *t, treeshape_norm norm);

/** Compute the total cophenetic index.
 *
 * See Mir, Arnau, and Francesc Rosselló. "A new balance index for phylogenetic
 * trees." Mathematical biosciences 241.1 (2013): 125-136.
 *
 * \param[in] tree tree to compute index for
 * \param[in] norm null model to use for normalizing the result
 * \return the sum of most-recent common ancestor depths for each pair of tips
 */
double cophenetic(const igraph_t *tree, treeshape_norm norm);

/** Compute the maximum ladder length.
 *
 * See \cite colijn2014phylogenetic.
 *
 * \param[in] tree tree to compute ladder length for
 * \return the maximum number of branches from the root to a tip, divided by
 * the number of tips
 */
double ladder_length(const igraph_t *tree);

/** Compute the number of IL nodes.
 *
 * See \cite colijn2014phylogenetic.
 *
 * \param[in] tree the tree to compute the number of IL nodes for
 * \return the proportion of internal nodes which have exactly one leaf child
 */
double il_nodes(const igraph_t *tree);

/** Compute the BMI of a tree.
 *
 * It's the maximum width divided by the maximum height. Haha, I'm so clever.
 * See \cite colijn2014phylogenetic.
 *
 * \param[in] tree tree to compute BMI for
 * \return maximum width divided by maximum height
 */
double bmi(const igraph_t *tree);

#endif
