/** \file simulate.h
 *
 * Functions for simulating phylogenies on contact networks.
 */

#ifndef SIMULATE_H
#define SIMULATE_H

#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

/** Simulate a phylogenetic tree from a contact network.
 *
 * This uses a simple algorithm. Let R be the sum of the rates of all
 * sero-discordant edges in the network, that is, edges where one endpoint is
 * infected and the other is susceptible. Whenever a new node is infected, we
 * update R accordingly.
 *
 * The epidemic is started by infecting a random node. The time to the next
 * infection is chosen from Exponential(R), a sero-discordant edge is chosen at
 * random, and the susceptible node along it is infected. The phylogeny is
 * constructed by keeping track of who transmitted to who.
 *
 * If the contact network is not connected, the phylogeny will only be
 * simulated on a (random) connected component. The stop_time and stop_tips
 * parameters are maximums, so the returned tree may have fewer nodes and/or be
 * shorter. This function produces a complete phylogeny (no tip sampling;
 * removed tips are included at their removal time).
 *
 * The contact network should have the following edge attributes:
 * * "transmit": the transmission rate along the edge
 *
 * The contact network should have the following vertex attributes:
 * * "remove": the removal rate of the node once it becomes infected (if the
 * node can recover or die, the removal rate is the sum of the recovery and
 * death rates)
 *
 * \param[in] net the contact network, with rates defined on each edge
 * \param[in] rng the GSL random generator object
 * \param[in] stop_time maximum amount of time to run the simulation for
 * \param[in] stop_nodes maximum number of nodes to infect
 * \return a phylogeny simulated over the contact network, or NULL
 */
igraph_t *simulate_phylogeny(igraph_t *net, gsl_rng *rng, 
        double stop_time, int stop_nodes);

/** Get the largest strongly connected component of a graph.
 *
 * \param[in] graph the graph to decompose
 * \return the graph's largest strongly connected component
 */
igraph_t *largest_component(igraph_t *graph);

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
