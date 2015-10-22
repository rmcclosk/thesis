#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <igraph/igraph.h>
#include <Judy.h>

#include "tree.h"
#include "util.h"
#include "treestats.h"

#define NDEBUG

/* helper functions for kernel */
int *production(const igraph_t *tree);
int count_node_pairs(const int *production1, const int *production2, int nnode1, 
        int nnode2);
int *get_node_pairs(const igraph_t *t1, const igraph_t *t2,
        int *production1, int *production2, int npairs);
int *children(const igraph_t *tree);
double *branch_lengths(const igraph_t *tree);
double Lp_norm(const double *x1, const double *x2, const double *y1, 
        const double *y2, int n1, int n2, double p);

double kernel(const igraph_t *t1, const igraph_t *t2, double decay_factor, 
        double rbf_variance, double sst_control)
{
    int i, c1, c2, n1, n2, coord, npairs;
    int *production1, *production2, *children1, *children2;
    double val, tmp, K = 0;
    double *bl1, *bl2;
    Pvoid_t delta = (Pvoid_t) NULL;
    int *pairs, *map, cur = 0;
    PWord_t Pvalue;
    Word_t bytes = 0;

    // preconditions
    assert(decay_factor > 0.0 && decay_factor <= 1.0);
    assert(rbf_variance > 0.0);
    assert(igraph_vcount(t1) < 65535);

    production1 = production(t1);
    production2 = production(t2);
    children1 = children(t1);
    children2 = children(t2);
    bl1 = branch_lengths(t1);
    bl2 = branch_lengths(t2);

    npairs = count_node_pairs(production1, production2, igraph_vcount(t1), igraph_vcount(t2));
    pairs = get_node_pairs(t1, t2, production1, production2, npairs);
    
    for (cur = 0; cur < npairs; ++cur)
    {
        val = decay_factor;
        n1 = pairs[cur] >> 16;
        n2 = pairs[cur] & 65535;

        // branch lengths
        tmp = pow(bl1[2*n1] - bl2[2*n2], 2) + pow(bl1[2*n1+1] - bl2[2*n2+1], 2);
        val *= exp(-tmp/rbf_variance);

        for (i = 0; i < 2; ++i)  // assume tree is binary
        {
            c1 = children1[2*n1+i];
            c2 = children2[2*n2+i];

            if (production1[c1] == production2[c2])
            {
                // children are leaves
                if (production1[c1] == 0)
                {
                    val *= (sst_control + decay_factor);
                }

                // children are not leaves
                else
                {
                    JLG(Pvalue, delta, (c1 << 16) | c2);
                    /* don't visit parents before children */
                    assert(Pvalue != NULL);
                    memcpy(&tmp, Pvalue, sizeof(double));
                    val *= (sst_control + tmp);
                }
            }
        }

        JLI(Pvalue, delta, pairs[cur]);
        if (Pvalue == PJERR) exit(EXIT_FAILURE); 
        memcpy(Pvalue, &val, sizeof(double));

        K += val;
    }

    free(production1);
    free(production2);
    free(children1);
    free(children2);
    free(bl1);
    free(bl2);
    free(pairs);
    JLFA(bytes, delta);
    return K;
}

double nLTT(const igraph_t *t1, const igraph_t *t2)
{
    int itree, i, cur;
    double h, k, prev;
    const igraph_t *trees[2] = {t1, t2};

    int n[2] = {(igraph_vcount(t1) + 1) / 2, (igraph_vcount(t2) + 1) / 2};
    double *x[2] = {calloc(n[0], sizeof(double)), calloc(n[1], sizeof(double))};
    double *y[2] = {calloc(n[0], sizeof(double)), calloc(n[1], sizeof(double))};
    double *buf = malloc(fmax(igraph_vcount(t1), igraph_vcount(t2)) * sizeof(double));
    int *node_order = malloc(fmax(igraph_vcount(t1), igraph_vcount(t2)) * sizeof(int));

    igraph_vector_t vec;
    igraph_vector_init(&vec, igraph_vcount(t1));

    for (itree = 0; itree < 2; ++itree)
    {
        depths(trees[itree], buf);
        order(buf, node_order, sizeof(double), igraph_vcount(trees[itree]),
                compare_doubles);
        igraph_degree(trees[itree], &vec, igraph_vss_all(), IGRAPH_OUT, 0);

        prev = 0; cur = 0; h = 0;
        y[itree][0] = 1.0 / n[itree];
        for (i = 0; i < igraph_vcount(trees[itree]); ++i)
        {
            if (VECTOR(vec)[node_order[i]] > 0) {
                if (buf[node_order[i]] == prev) {
                    fprintf(stderr, "%f same as prev\n", buf[node_order[i]]);
                    y[itree][cur] += 1.0 / n[itree];
                }
                else {
                    fprintf(stderr, "%f different from prev\n", buf[node_order[i]]);
                    x[itree][++cur] = buf[node_order[i]];
                    y[itree][cur] = y[itree][cur-1] + 1.0 / n[itree];
                    h = fmax(h, buf[node_order[i]]);
                    prev = buf[node_order[i]];
                }
            }
        }

        for (i = 0; i < n[itree]; ++i) {
            x[itree][i] /= h;
        }
    }

    k = Lp_norm(x[0], x[1], y[0], y[1], n[0], n[1], 1.0);

    free(x[0]);
    free(x[1]);
    free(y[0]);
    free(y[1]);
    free(buf);
    free(node_order);
    igraph_vector_destroy(&vec);
    return k;
}

/* Private. */

double Lp_norm(const double *x1, const double *x2, const double *y1, 
        const double *y2, int n1, int n2, double p)
{
    double a, b, fa1, fa2, fb1, fb2, m1, m2, intersect, fintersect, area, norm = 0;
    int i1 = 0, i2 = 0;

    if (x1[0] < x2[0])
    {
        b = x1[0];
        fb1 = y1[0];
        fb2 = ((y2[1] - y2[0]) / (x2[1] - x2[0])) * (b - x2[0]);
    }
    else if (x1[0] == x2[0])
    {
        b = x1[0];
        fb1 = y1[0];
        fb2 = y2[0];
    }
    else
    {
        b = x2[0];
        fb2 = y2[0];
        fb1 = ((y1[1] - y1[0]) / (x1[1] - x1[0])) * (b - x1[0]);
    }

    while (i1 < n1 && i2 < n2)
    {
        a = b;
        b = fmin(x1[i1], x2[i2]);
        fa1 = fb1;
        fa2 = fb2;

        if (x1[i1] < x2[i2])
        {
            fb1 = y1[i1++];
            fb2 = fa2 + (y2[i2] - fa2) / (x2[i2] - a) * (b-a);
        }
        else if (x1[i1] == x2[i2])
        {
            fb1 = y1[i1++];
            fb2 = y2[i2++];
        }
        else
        {
            fb2 = y2[i2++];
            fb1 = fa1 + (y1[i1] - fa1) / (x1[i1] - a) * (b-a);
        }

        m1 = (fb1 - fa1) / (b - a);
        m2 = (fb2 - fa2) / (b - a);

        if (m1 != m2)
        {
            intersect = a + (fa1 - fa2) / (m2 - m1);
        }
        else
        {
            intersect = a;
        }

        if (intersect > a && intersect < b)
        {
            fintersect = fa1 + m1 * intersect;
            area  = fabs(((fa1 + fintersect) / 2 - (fa2 + fintersect) / 2) * (intersect-a));
            area += fabs(((fintersect + fb1) / 2 - (fintersect + fb2) / 2) * (b-intersect));
        }
        else
        {
            area = fabs(((fa1 + fb1) / 2 - (fa2 + fb2) / 2) * (b-a));
        }
        norm += area;
    }
    return pow(norm, p);
}

/* get production rules for each node */
int *production(const igraph_t *tree)
{
    int i, nnode = igraph_vcount(tree);
    int *p = malloc(nnode * sizeof(int));
    igraph_vector_t vec, nbr;
    igraph_vector_init(&nbr, 2);
    igraph_vector_init(&vec, 2);

    igraph_degree(tree, &vec, igraph_vss_all(), IGRAPH_OUT, 0);
    for (i = 0; i < nnode; ++i)
    {
        if ((int) VECTOR(vec)[i] == 0)
        {
            p[i] = 0;
        }
        else
        {
            igraph_neighbors(tree, &nbr, i, IGRAPH_OUT);
            p[i] = ((int) VECTOR(vec)[(int) VECTOR(nbr)[0]] == 0) +
                   ((int) VECTOR(vec)[(int) VECTOR(nbr)[1]] == 0) + 1;
        }
    }
    igraph_vector_destroy(&nbr);
    igraph_vector_destroy(&vec);
    return p;
}

/* find how many pairs of nodes we need to evaluate */
int count_node_pairs(const int *production1, const int *production2, int nnode1, 
        int nnode2)
{
    int i, npairs = 0, table1[4] = {0}, table2[4] = {0};

    for (i = 0; i < nnode1; ++i) {
        table1[production1[i]] += 1;
    }
    for (i = 0; i < nnode2; ++i) {
        table2[production2[i]] += 1;
    }

    for (i = 1; i < 4; ++i) {
        npairs += table1[i] * table2[i];
    }
    return npairs;
}

/* get all pairs of nodes with the same production */
int *get_node_pairs(const igraph_t *t1, const igraph_t *t2,
        int *production1, int *production2, int npairs)
{
    int coord, n1, n2, i1 = 1, i2 = 1, start = 0;
    int nnode1 = igraph_vcount(t1), nnode2 = igraph_vcount(t2);
    int *order1, *order2;
    int *pairs, i = 0;

    order1 = malloc(nnode1 * sizeof(int));
    order2 = malloc(nnode2 * sizeof(int));
    order(production1, order1, sizeof(int), nnode1, compare_ints);
    order(production2, order2, sizeof(int), nnode2, compare_ints);

    pairs = malloc(npairs * sizeof(int));

    n1 = order1[0]; n2 = order2[0];

    // find pairs of nodes with equal productions
    while (i1 < nnode1 || i2 < nnode2) 
    {
        while (production1[n1] != production2[n2])
        {
            if (production1[n1] > production2[n2])
            {
                if (i2 < nnode2)
                    n2 = order2[i2++];
                else
                    break;
            }
            else if (production1[n1] < production2[n2])
            {
                if (i1 < nnode1)
                    n1 = order1[i1++];
                else
                    break;
            }
        } 
        start = i2 - 1;

        if (production1[n1] != production2[n2])
            break;

        while (production1[n1] == production2[n2])
        {
            while (production1[n1] == production2[n2])
            {
                coord = (n1 << 16) | n2;

                // internal nodes
                if (production1[n1] > 0) 
                    pairs[i++] = coord;

                if (i2 < nnode2) 
                    n2 = order2[i2++];
                else 
                    break;
            }
            if (i1 < nnode1)
            {
                n1 = order1[i1++];
                i2 = start;
                n2 = order2[i2++];
            }
            else
            {
                break;
            }
        }
    }
    free(order1);
    free(order2);

    qsort(pairs, npairs, sizeof(int), compare_ints);
    return pairs;
}

/* get indices of children from each node */
int *children(const igraph_t *tree)
{
    igraph_adjlist_t al;
    igraph_vector_int_t *nbr;
    int i;
    int *children = malloc(igraph_vcount(tree) * 2 * sizeof(int));

    igraph_adjlist_init(tree, &al, IGRAPH_OUT);
    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        nbr = igraph_adjlist_get(&al, i);
        if (igraph_vector_int_size(nbr) > 0)
        {
            children[2*i] = VECTOR(*nbr)[0];
            children[2*i+1] = VECTOR(*nbr)[1];
        }
    }

    igraph_adjlist_destroy(&al);
    return children;
}

/* get branch lengths leading out of each node */
double *branch_lengths(const igraph_t *tree)
{
    igraph_inclist_t il;
    int i;
    igraph_vector_int_t *edge;
    double *branch_lengths = malloc(2 * igraph_vcount(tree) * sizeof(double));

    igraph_inclist_init(tree, &il, IGRAPH_OUT);
    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        edge = igraph_inclist_get(&il, i);
        if (igraph_vector_int_size(edge) > 0)
        {
            branch_lengths[2*i] = EAN(tree, "length", VECTOR(*edge)[0]);
            branch_lengths[2*i+1] = EAN(tree, "length", VECTOR(*edge)[1]);
        }
    }

    igraph_inclist_destroy(&il);
    return branch_lengths;
}

