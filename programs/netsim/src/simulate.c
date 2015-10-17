#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <igraph/igraph.h>
#include <gsl/gsl_randist.h>

#include "simulate.h"

void add_vertex_id(const igraph_t *net, igraph_strvector_t *v, int node, int numeric_ids);

void simulate_phylogeny(igraph_t *tree, igraph_t *net, gsl_rng *rng,
        double stop_time, int stop_nodes, int numeric_ids)
{
    long i;
    int inode, snode, e, v, head, tail, nnode_tree = 0;
    double t, trans_rate = 0., remove_rate = 0., time = 0.;
    char buf[128];
    igraph_vector_t discordant, incident, edges, branch_lengths, tip_map;
    igraph_vector_int_t infected, removed; 
    igraph_strvector_t node_map;

    igraph_vector_int_init(&infected, 0);
    igraph_vector_int_init(&removed, 0);
    igraph_vector_init(&tip_map, 0);
    igraph_strvector_init(&node_map, 0);
    igraph_vector_init(&discordant, 0);
    igraph_vector_init(&incident, 0);
    igraph_vector_init(&edges, 0);
    igraph_vector_init(&branch_lengths, 0);

    // start the epidemic
    inode = rand() % igraph_vcount(net);
    igraph_vector_int_push_back(&infected, inode);
    igraph_vector_push_back(&branch_lengths, 0.);
    igraph_vector_push_back(&tip_map, nnode_tree++);
    add_vertex_id(net, &node_map, inode, numeric_ids);

    // the initial node's incident edges are the discordant edges now
    igraph_incident(net, &discordant, inode, IGRAPH_OUT);
    igraph_vector_sort(&discordant);

    for (e = 0; e < igraph_vector_size(&discordant); ++e)
        trans_rate += EAN(net, "transmit", (int) VECTOR(discordant)[e]);
    remove_rate = VAN(net, "remove", inode);

    // simulate until either we reach the time goal, or everybody is infected
    while (!igraph_vector_empty(&discordant) && time < stop_time &&
            (nnode_tree + 1) / 2 < stop_nodes)
    {
        // choose the next event time
        t = gsl_ran_exponential(rng, 1. / (trans_rate + remove_rate));
        if (t + time < stop_time)
        {
            // increase all extant branch lengths by t
            for (v = 0; v < igraph_vector_int_size(&infected); ++v)
                VECTOR(branch_lengths)[(int) VECTOR(tip_map)[v]] += t;

            // next event is a transmission
            if ((double) rand() / (double) RAND_MAX < trans_rate / (trans_rate + remove_rate))
            {
                // choose the next edge
                i = rand() % igraph_vector_size(&discordant);
                igraph_edge(net, VECTOR(discordant)[i], &inode, &snode);
                //fprintf(stderr, "infect node %s->%s\n", VAS(net, "id", inode), 
                                //VAS(net, "id", snode));

                // mark the new node as infected
                igraph_vector_int_binsearch(&infected, snode, &i);
                igraph_vector_int_insert(&infected, i, snode);
                igraph_vector_insert(&tip_map, i, nnode_tree++);
                add_vertex_id(net, &node_map, snode, numeric_ids);

                // add edges in the tree for the new transmission
                igraph_vector_int_binsearch(&infected, inode, &i);
                igraph_vector_push_back(&edges, VECTOR(tip_map)[i]);
                igraph_vector_push_back(&edges, nnode_tree-1);
                igraph_vector_push_back(&edges, VECTOR(tip_map)[i]);
                igraph_vector_set(&tip_map, i, nnode_tree++);
                igraph_vector_push_back(&edges, nnode_tree-1);
                remove_rate += VAN(net, "remove", snode);
                add_vertex_id(net, &node_map, inode, numeric_ids);
    
                // the new edges get branch length zero
                igraph_vector_push_back(&branch_lengths, 0.);
                igraph_vector_push_back(&branch_lengths, 0.);

                // outgoing edges to susceptible nodes are discordant now
                igraph_incident(net, &incident, snode, IGRAPH_OUT);
                for (e = 0; e < igraph_vector_size(&incident); ++e)
                {
                    igraph_edge(net, VECTOR(incident)[e], &tail, &head);
                    if (!igraph_vector_int_binsearch2(&removed, head) &&
                        !igraph_vector_int_binsearch2(&infected, head))
                    {
                        //printf("add edge %d->%d\n", tail, head);
                        igraph_vector_binsearch(&discordant, VECTOR(incident)[e], &i);
                        igraph_vector_insert(&discordant, i, VECTOR(incident)[e]);
                        trans_rate += EAN(net, "transmit", (int) VECTOR(incident)[e]);
                    }
                }

                // incoming edges which were discordant before aren't anymore
                igraph_incident(net, &incident, snode, IGRAPH_IN);
                for (e = 0; e < igraph_vector_size(&incident); ++e)
                {
                    igraph_edge(net, VECTOR(incident)[e], &tail, &head);
                    if (igraph_vector_binsearch(&discordant, VECTOR(incident)[e], &i))
                    {
                        //printf("remove edge %d->%d\n", tail, head);
                        igraph_vector_remove(&discordant, i);
                        trans_rate -= EAN(net, "transmit", (int) VECTOR(incident)[e]);
                    }
                }
            }

            // next event is a removal
            else
            {
                // choose and remove node
                i = rand() % igraph_vector_int_size(&infected);
                inode = VECTOR(infected)[i];
                igraph_vector_int_remove(&infected, i);
                igraph_vector_int_binsearch(&removed, inode, &i);
                igraph_vector_int_insert(&removed, i, inode);
                remove_rate -= VAN(net, "remove", inode);
                //printf("remove node %d\n", inode);

                // outgoing discordant edges aren't discordant anymore
                igraph_incident(net, &incident, inode, IGRAPH_OUT);
                for (e = 0; e < igraph_vector_size(&incident); ++e)
                {
                    if (igraph_vector_binsearch(&discordant, VECTOR(incident)[e], &i))
                    {
                        igraph_edge(net, VECTOR(incident)[e], &tail, &head);
                        //printf("remove edge %d->%d\n", tail, head);
                        igraph_vector_remove(&discordant, i);
                        trans_rate -= EAN(net, "transmit", (int) VECTOR(incident)[e]);
                    }
                }
            }
            time += t;
        }
        else
        {
            // increase all extant branch lengths up to stop_time
            for (v = 0; v < igraph_vector_int_size(&infected); ++v)
                VECTOR(branch_lengths)[(int) VECTOR(tip_map)[v]] += (stop_time - time);
            time = stop_time;
        }
    }

    // assemble the tree
    igraph_empty(tree, nnode_tree, IGRAPH_DIRECTED);
    igraph_add_edges(tree, &edges, 0);
    for (e = 0; e < igraph_ecount(tree); ++e)
    {
        igraph_edge(tree, e, &head, &tail);
        SETEAN(tree, "length", e, VECTOR(branch_lengths)[tail]);
    }
    for (v = 0; v < igraph_strvector_size(&node_map); ++v)
    {
        SETVAS(tree, "id", v, STR(node_map, v));
    }

    // clean up
    igraph_vector_int_destroy(&infected);
    igraph_vector_int_destroy(&removed);
    igraph_strvector_destroy(&node_map);
    igraph_vector_destroy(&tip_map);
    igraph_vector_destroy(&discordant);
    igraph_vector_destroy(&incident);
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&branch_lengths);
}

/* Private */

void add_vertex_id(const igraph_t *net, igraph_strvector_t *v, int node, int numeric_ids)
{
    char buf[128];
    if (igraph_cattribute_has_attr(net, IGRAPH_ATTRIBUTE_VERTEX, "id")) {
        if (numeric_ids)
            sprintf(buf, "%d", (int) VAN(net, "id", node));
        else
            sprintf(buf, "%s", VAS(net, "id", node));
    }
    else {
        sprintf(buf, "%d", node);
    }
    igraph_strvector_add(v, buf);
}
