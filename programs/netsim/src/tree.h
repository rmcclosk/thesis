/** \file tree.h
 * \brief Functions for handling phylogenetic trees.
 */

#ifndef TREE_H
#define TREE_H

#include <stdio.h>

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
void write_tree_newick(igraph_t *tree, FILE *f);

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
