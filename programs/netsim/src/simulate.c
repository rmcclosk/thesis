#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <igraph/igraph.h>
#include <gsl/gsl_randist.h>

#include "simulate.h"

void _cut_at_time(igraph_t *tree, double t, int root, double troot, 
        int extant_only, igraph_vector_t *work, igraph_vector_t *to_delete);
int _collapse_singles(igraph_t *tree, int root, igraph_vector_t *vdel,
        igraph_vector_t *eadd, igraph_vector_t *branch_length,
        igraph_vector_t *work, double *bl);

igraph_t *largest_component(igraph_t *graph)
{
    igraph_vector_ptr_t *comps = malloc(sizeof(igraph_vector_ptr_t));
    igraph_vector_ptr_init(comps, 1);
    igraph_t *subg = malloc(sizeof(igraph_t));
    int i, max, vmax = 0;

    igraph_decompose(graph, comps, IGRAPH_WEAK, -1, 0);
    for (i = 0; i < igraph_vector_ptr_size(comps); ++i)
    {
        if (igraph_vcount(VECTOR(*comps)[i]) > vmax)
        {
            vmax = igraph_vcount(VECTOR(*comps)[i]);
            max = i;
        }
    }
    igraph_copy(subg, VECTOR(*comps)[max]);
    igraph_decompose_destroy(comps);

    return subg;
}

void cut_at_time(igraph_t *tree, double t, int extant_only)
{
    int i;
    igraph_vector_t work, to_delete;
    igraph_vector_init(&work, 0);
    igraph_vector_init(&to_delete, 0);

    _cut_at_time(tree, t, root(tree), 0., extant_only, &work, &to_delete);
    igraph_vector_sort(&to_delete);
    igraph_delete_vertices(tree, igraph_vss_vector(&to_delete));
    collapse_singles(tree);

    igraph_vector_destroy(&work);
    igraph_vector_destroy(&to_delete);
}

void collapse_singles(igraph_t *tree)
{
    int i, edge;
    double bl = 0.;
    igraph_vector_t vdel, eadd, work, branch_length;

    igraph_vector_init(&vdel, 0);
    igraph_vector_init(&eadd, 0);
    igraph_vector_init(&work, 0);
    igraph_vector_init(&branch_length, 0);

    _collapse_singles(tree, root(tree), &vdel, &eadd, &branch_length, &work, &bl);
    igraph_add_edges(tree, &eadd, 0);
    for (i = 0; i < igraph_vector_size(&eadd)/2; ++i)
    {
        igraph_get_eid(tree, &edge, VECTOR(eadd)[2*i], VECTOR(eadd)[2*i+1], 1, 1);
        SETEAN(tree, "length", edge, VECTOR(branch_length)[i]);
    }
    igraph_delete_vertices(tree, igraph_vss_vector(&vdel));

    igraph_vector_destroy(&branch_length);
    igraph_vector_destroy(&vdel);
    igraph_vector_destroy(&eadd);
    igraph_vector_destroy(&work);
}

void subsample_tips(igraph_t *tree, int ntip, const gsl_rng *rng)
{
    int cur = 0, nnode = igraph_vcount(tree); 
    igraph_real_t i;
    int ndrop = (nnode + 1)/2 - ntip;
    igraph_real_t *tips;
    igraph_real_t *drop; 
    igraph_vector_t degree;
    igraph_vector_t view;
    
    if (ndrop <= 0)
        return;

    drop = malloc(ndrop * sizeof(igraph_real_t));
    tips = malloc((nnode + 1)/2 * sizeof(igraph_real_t));
    igraph_vector_init(&degree, 0);

    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);
    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        if ((int) VECTOR(degree)[(int) i] == 0)
            tips[cur++] = i;
    }
    gsl_ran_choose(rng, drop, ndrop, tips, (nnode + 1)/2, sizeof(igraph_real_t));
    igraph_vector_view(&view, drop, ndrop);
    igraph_delete_vertices(tree, igraph_vss_vector(&view));
    collapse_singles(tree);

    igraph_vector_destroy(&degree);
    free(tips);
    free(drop);
}

igraph_t *simulate_phylogeny(igraph_t *net, gsl_rng *rng, double stop_time, 
        int stop_nodes)
{
    long i;
    int inode, snode, e, v, head, tail, nnode_tree = 0;
    double t, trans_rate = 0., remove_rate = 0., time = 0.;
    int tmp;
    igraph_t *tree = malloc(sizeof(igraph_t));
    igraph_vector_t discordant, incident, edges, branch_lengths, node_map;
    igraph_vector_int_t infected, removed;

    igraph_vector_int_init(&infected, 0);
    igraph_vector_int_init(&removed, 0);
    igraph_vector_init(&node_map, 0);
    igraph_vector_init(&discordant, 0);
    igraph_vector_init(&incident, 0);
    igraph_vector_init(&edges, 0);
    igraph_vector_init(&branch_lengths, 0);

    // start the epidemic
    inode = rand() % igraph_vcount(net);
    igraph_vector_int_push_back(&infected, inode);
    igraph_vector_push_back(&branch_lengths, 0.);
    igraph_vector_push_back(&node_map, nnode_tree++);

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
                VECTOR(branch_lengths)[(int) VECTOR(node_map)[v]] += t;

            // next event is a transmission
            if ((double) rand() / (double) RAND_MAX < trans_rate / (trans_rate + remove_rate))
            {
                // choose the next edge
                i = rand() % igraph_vector_size(&discordant);
                igraph_edge(net, VECTOR(discordant)[i], &inode, &snode);

                // mark the new node as infected
                igraph_vector_int_binsearch(&infected, snode, &i);
                igraph_vector_int_insert(&infected, i, snode);
                igraph_vector_insert(&node_map, i, nnode_tree++);

                // add edges in the tree for the new transmission
                igraph_vector_int_binsearch(&infected, inode, &i);
                igraph_vector_push_back(&edges, VECTOR(node_map)[i]);
                igraph_vector_push_back(&edges, nnode_tree-1);
                igraph_vector_push_back(&edges, VECTOR(node_map)[i]);
                igraph_vector_set(&node_map, i, nnode_tree++);
                igraph_vector_push_back(&edges, nnode_tree-1);
                remove_rate += VAN(net, "remove", snode);
    
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

                // outgoing discordant edges aren't discordant anymore
                igraph_incident(net, &incident, inode, IGRAPH_OUT);
                for (e = 0; e < igraph_vector_size(&incident); ++e)
                {
                    if (igraph_vector_binsearch(&discordant, VECTOR(incident)[e], &i))
                    {
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
                VECTOR(branch_lengths)[(int) VECTOR(node_map)[v]] += (stop_time - time);
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

    // clean up
    igraph_vector_int_destroy(&infected);
    igraph_vector_int_destroy(&removed);
    igraph_vector_destroy(&node_map);
    igraph_vector_destroy(&discordant);
    igraph_vector_destroy(&incident);
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&branch_lengths);

    return tree;
}

/* Private */

void _cut_at_time(igraph_t *tree, double t, int root, double troot, 
        int extant_only, igraph_vector_t *work, igraph_vector_t *to_delete)
{
    int i, lc, rc;
    double tnode;

    igraph_incident(tree, work, root, IGRAPH_IN);
    if (igraph_vector_size(work) > 0)
        tnode = EAN(tree, "length", VECTOR(*work)[0]);
    else
        tnode = 0.;

    if (tnode + troot > t)
    {
        // adjust my branch length
        SETEAN(tree, "length", VECTOR(*work)[0], t - troot);

        // delete all descendents
        igraph_subcomponent(tree, work, root, IGRAPH_OUT);
        for (i = 0; i < igraph_vector_size(work); ++i) 
        {
            if ((int) VECTOR(*work)[i] != root)
                igraph_vector_push_back(to_delete, VECTOR(*work)[i]);
        }
    }
    else
    {
        igraph_neighbors(tree, work, root, IGRAPH_OUT);
        if (igraph_vector_size(work) > 0)
        {
            lc = (int) VECTOR(*work)[0];
            rc = (int) VECTOR(*work)[1];

            _cut_at_time(tree, t, lc, troot + tnode, extant_only, work, to_delete);
            _cut_at_time(tree, t, rc, troot + tnode, extant_only, work, to_delete);
        }

        else if (extant_only && troot + tnode < t)
        {
            igraph_vector_push_back(to_delete, (igraph_real_t) root);
        }
    }
}

int _collapse_singles(igraph_t *tree, int root, igraph_vector_t *vdel,
        igraph_vector_t *eadd, igraph_vector_t *branch_length, 
        igraph_vector_t *work, double *bl)
{
    int lc, rc, new_lc, new_rc;
    double lbl, rbl;

    igraph_incident(tree, work, root, IGRAPH_IN);
    if (igraph_vector_size(work) > 0)
        bl[0] = EAN(tree, "length", (int) VECTOR(*work)[0]);
    else
        bl[0] = 0;

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    if (igraph_vector_size(work) == 0)
    {
        return root;
    }
    else if (igraph_vector_size(work) == 1)
    {
        igraph_vector_push_back(vdel, root);
        lc = (int) VECTOR(*work)[0];
        new_lc = _collapse_singles(tree, lc, vdel, eadd, branch_length, work, &lbl);
        bl[0] += lbl;
        return new_lc;
    }
    else
    {
        lc = (int) VECTOR(*work)[0];
        rc = (int) VECTOR(*work)[1];
        new_lc = _collapse_singles(tree, lc, vdel, eadd, branch_length, work, &lbl);
        new_rc = _collapse_singles(tree, rc, vdel, eadd, branch_length, work, &rbl);

        if (new_lc != lc)
        {
            igraph_vector_push_back(eadd, root);
            igraph_vector_push_back(eadd, new_lc);
            igraph_vector_push_back(branch_length, lbl);
        }
        if (new_rc != rc)
        {
            igraph_vector_push_back(eadd, root);
            igraph_vector_push_back(eadd, new_rc);
            igraph_vector_push_back(branch_length, rbl);
        }

        return root;
    }
}
