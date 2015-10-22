#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <Judy.h>
#include <igraph/igraph.h>
#include <gsl/gsl_randist.h>

#include "simulate.h"
#define NDEBUG

void print_node(const igraph_t *net, char *buf, int node, int numeric_ids);

void simulate_phylogeny(igraph_t *tree, igraph_t *net, gsl_rng *rng,
        double stop_time, int stop_nodes, int numeric_ids)
{
    long i;
    int inode, snode, e, v, n, head, tail, nnode_tree = 0;
    int Rc_int, ndiscordant, ninfected = 1;
    double t, trans_rate = 0., remove_rate = 0., time = 0.;
    char buf[128];
    igraph_vector_int_t *incident;
    igraph_vector_t edges, branch_lengths;
    igraph_inclist_t inclist_in, inclist_out;
    Word_t Bytes, Index, *PValue;

    // set of infected nodes
    Pvoid_t infected = (Pvoid_t) NULL;

    // set of removed nodes
    Pvoid_t removed = (Pvoid_t) NULL;

    // set of discordant edges
    Pvoid_t discordant = (Pvoid_t) NULL;

    // map from nodes in the network to extant tips in the tree
    Pvoid_t tip_map = (Pvoid_t) NULL;

    // map from nodes (internal + terminal) in the tree to nodes in the network
    Pvoid_t node_map = (Pvoid_t) NULL;

    igraph_inclist_init(net, &inclist_in, IGRAPH_IN);
    igraph_inclist_init(net, &inclist_out, IGRAPH_OUT);

    igraph_vector_init(&edges, 0);
    igraph_vector_init(&branch_lengths, 0);

    // start the epidemic
    inode = rand() % igraph_vcount(net);
    J1S(Rc_int, infected, inode);
    //fprintf(stderr, "start epidemic at node %d\n", inode);
    igraph_vector_push_back(&branch_lengths, 0.);
    JLI(PValue, tip_map, inode); *PValue = nnode_tree++;
    JLI(PValue, node_map, nnode_tree-1); *PValue = inode;

    // the initial node's incident edges are the discordant edges now
    incident = igraph_inclist_get(&inclist_out, inode);

    for (e = 0; e < igraph_vector_int_size(incident); ++e) {
        J1S(Rc_int, discordant, VECTOR(*incident)[e]);
    }

    ndiscordant = igraph_vector_int_size(incident);
    Index = 0; J1F(Rc_int, discordant, Index);
    for (e = 0; e < ndiscordant; ++e) {
        trans_rate += EAN(net, "transmit", Index);
        J1N(Rc_int, discordant, Index);
    }
    remove_rate = VAN(net, "remove", inode);

    // simulate until either we reach the time goal, or everybody is infected
    while (ndiscordant > 0 && time < stop_time && (nnode_tree + 1) / 2 < stop_nodes)
    {
        // choose the next event time
        t = gsl_ran_exponential(rng, 1. / (trans_rate + remove_rate));
        if (t + time < stop_time)
        {
            // increase all extant branch lengths by t
            Index = 0; J1F(Rc_int, infected, Index);
            for (v = 0; v < ninfected; ++v) {
                JLG(PValue, tip_map, Index);
                VECTOR(branch_lengths)[*PValue] += t;
                J1N(Rc_int, infected, Index);
            }

            // next event is a transmission
            if ((double) rand() / (double) RAND_MAX < trans_rate / (trans_rate + remove_rate))
            {
                // choose the next edge
                i = rand() % ndiscordant;
                J1BC(Rc_int, discordant, i+1, Index);
                igraph_edge(net, Index, &inode, &snode);
                //fprintf(stderr, "infect node %d->%d\n", inode, snode);

                // mark the new node as infected
                ++ninfected;
                J1S(Rc_int, infected, snode);
                JLI(PValue, tip_map, snode); *PValue = nnode_tree++;
                JLI(PValue, node_map, nnode_tree-1); *PValue = snode;

                // add edges in the tree for the new transmission
                JLG(PValue, tip_map, inode);
                assert(PValue != NULL);
                igraph_vector_push_back(&edges, *PValue);
                igraph_vector_push_back(&edges, nnode_tree-1);
                igraph_vector_push_back(&edges, *PValue);
                igraph_vector_push_back(&edges, nnode_tree);

                JLI(PValue, tip_map, inode); *PValue = nnode_tree++;
                JLI(PValue, node_map, nnode_tree-1); *PValue = inode;
                remove_rate += VAN(net, "remove", snode);
    
                // the new edges get branch length zero
                igraph_vector_push_back(&branch_lengths, 0.);
                igraph_vector_push_back(&branch_lengths, 0.);

                // outgoing edges to susceptible nodes are discordant now
                incident = igraph_inclist_get(&inclist_out, snode);
                n = igraph_vector_int_size(incident);
                for (e = 0; e < n; ++e)
                {
                    igraph_edge(net, VECTOR(*incident)[e], &tail, &head);
                    J1T(Rc_int, infected, head);
                    if (!Rc_int) 
                    {
                        J1T(Rc_int, removed, head);
                        if (!Rc_int)
                        {
                            //printf("add edge %d->%d\n", tail, head);
                            J1S(Rc_int, discordant, VECTOR(*incident)[e]);
                            trans_rate += EAN(net, "transmit", VECTOR(*incident)[e]);
                            ++ndiscordant;
                        }
                    }
                }

                // incoming edges which were discordant before aren't anymore
                incident = igraph_inclist_get(&inclist_in, snode);
                n = igraph_vector_int_size(incident);
                for (e = 0; e < n; ++e)
                {
                    igraph_edge(net, VECTOR(*incident)[e], &tail, &head);
                    J1T(Rc_int, discordant, VECTOR(*incident)[e]);
                    if (Rc_int)
                    {
                        //printf("remove edge %d->%d\n", tail, head);
                        J1U(Rc_int, discordant, VECTOR(*incident)[e]);
                        trans_rate -= EAN(net, "transmit", VECTOR(*incident)[e]);
                        --ndiscordant;
                    }
                }
            }

            // next event is a removal
            else
            {
                // choose and remove node
                i = rand() % ninfected;
                J1BC(Rc_int, infected, i+1, Index);
                J1U(Rc_int, infected, Index);
                J1S(Rc_int, removed, Index);
                remove_rate -= VAN(net, "remove", Index);
                --ninfected;
                //printf("remove node %d\n", Index);

                // outgoing discordant edges aren't discordant anymore
                incident = igraph_inclist_get(&inclist_out, Index);
                for (e = 0; e < igraph_vector_int_size(incident); ++e)
                {
                    J1T(Rc_int, discordant, VECTOR(*incident)[e]);
                    if (Rc_int)
                    {
                        igraph_edge(net, VECTOR(*incident)[e], &tail, &head);
                        //printf("remove edge %d->%d\n", tail, head);
                        J1U(Rc_int, discordant, VECTOR(*incident)[e]);
                        trans_rate -= EAN(net, "transmit", VECTOR(*incident)[e]);
                        --ndiscordant;
                    }
                }
            }
            time += t;
        }
        else
        {
            Index = 0; J1F(Rc_int, infected, Index);
            for (v = 0; v < ninfected; ++v) {
                JLG(PValue, tip_map, Index);
                VECTOR(branch_lengths)[*PValue] += (stop_time - time);
                J1N(Rc_int, infected, Index);
            }
            time = stop_time;
        }
    }

    // assemble the tree
    igraph_empty(tree, nnode_tree, IGRAPH_DIRECTED);
    igraph_add_edges(tree, &edges, 0);
    n = igraph_ecount(tree);
    for (e = 0; e < n; ++e)
    {
        igraph_edge(tree, e, &head, &tail);
        SETEAN(tree, "length", e, VECTOR(branch_lengths)[tail]);
    }
    for (v = 0; v < nnode_tree; ++v)
    {
        JLG(PValue, node_map, v);
        print_node(net, buf, *PValue, numeric_ids);
        SETVAS(tree, "id", v, buf);
    }

    // clean up
    JLFA(Bytes, tip_map);
    JLFA(Bytes, node_map);
    J1FA(Bytes, infected);
    J1FA(Bytes, discordant);
    J1FA(Bytes, removed);
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&branch_lengths);
    igraph_inclist_destroy(&inclist_in);
    igraph_inclist_destroy(&inclist_out);
}

/* Private */

void print_node(const igraph_t *net, char *buf, int node, int numeric_ids)
{
    if (igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_VERTEX, "id")) {
        if (numeric_ids)
            sprintf(buf, "%d", (int) VAN(net, "id", node));
        else
            sprintf(buf, "%s", VAS(net, "id", node));
    }
    else {
        sprintf(buf, "%d", node);
    }
}
