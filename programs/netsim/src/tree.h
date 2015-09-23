/** \file tree.h
 * \brief Functions for handling phylogenetic trees.
 */

#ifndef TREE_H
#define TREE_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <igraph/igraph.h>

typedef enum {
    MEAN,
    MEDIAN,
    NONE
} scaling;

/** Parse a Newick tree.
 *
 * \param[in] f open file handle to a file containing a Newick tree string
 * \return the tree represented by the Newick string
 */
igraph_t *parse_newick(FILE *f);

/** Output a tree in Newick format.
 *
 * \param[in] tree the tree to output
 * \param[in] f file handle open for writing
 */
void write_tree_newick(const igraph_t *tree, FILE *f);

/** Get the root of a tree.
 * 
 * The root is defined as the unique vertex with zero in-degree.
 *
 * \param[in] tree the tree to find the root of
 * \return the index of the tree's root
 */
int root(const igraph_t *tree);

/** Calculate the height of a tree.
 *
 * The tree's edges must have the "length" attribute defined on them. This
 * function is recursive internally, so if the passed graph is not a tree (ie.
 * if there are cycles), it will loop forever.
 *
 * \param[in] tree tree to calculate the height of
 * \return the height of the tree
 */
double height(const igraph_t *tree);

/** Get the depths of all nodes in a tree.
 *
 * If the tree has n nodes, the depths array must have enough space allocated
 * for n doubles.
 *
 * \param[in] tree tree to find depths for
 * \param[out] depths the result will be stored here
 */
void depths(const igraph_t *tree, double *depths);

/** Ladderize a tree.
 *
 * This renumbers the vertices of the tree such that:
 *
 * * children always have larger indices than their parents,
 * * the sibling with the most descendants has the larger index,
 * * if two siblings have an equal number of descendants, the one with the
 * largest branch length to its parent has the larger index.
 *
 * \param[in,out] tree the tree to ladderize
 */
void ladderize(igraph_t *tree);

/** Scale the branches of a tree.
 *
 * Scale down the branches in a tree according to the specified scaling
 * mode.
 *
 * \param[in,out] tree the tree to scale
 * \param[in] mode how to scale the branches
 */
void scale_branches(igraph_t *tree, scaling mode);

/** Cut a tree at a specified time.
 *
 * Always deletes nodes which were born after the cutoff. Optionally, tips
 * which were removed prior to the cutoff are also deleted.
 *
 * \param[in] tree the tree to shorten
 * \param[in] t the time to cut the tree (measured forward from the root)
 * \param[in] extant_only if true, only include tips which are extant at the
 * cutoff time
 */
void cut_at_time(igraph_t *tree, double t, int extant_only);

/** Collapse single nodes in a tree.
 *
 * This deletes nodes which have in-degree and out-degree both equal to 1, and
 * connects their parents to their children.
 *
 * \param[in] tree tree to collapse singles in
 */
void collapse_singles(igraph_t *tree);

/** Subsample tips from a tree.
 * 
 * Randomly deletes tips from a tree until there are only ntip tips remaining.
 * If there are ntip or fewer tips in the tree already, nothing is done.
 *
 * \param[in] tree the tree to subsample
 * \param[in] ntip number of tips to leave in the tree
 * \param[in] rng GSL random number generator object
 */
void subsample_tips(igraph_t *tree, int ntip, const gsl_rng *rng);

#endif
