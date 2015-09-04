#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <gsl/gsl_randist.h>
#include <igraph/igraph.h>
#include <Judy.h>

#include "newick_parser.h"
#include "tree.h"
#include "util.h"

#define NDEBUG

void yyrestart(FILE *f);
int _ladderize(igraph_t *tree, igraph_vector_t *work, int root, int *perm);
double _height(const igraph_t *tree, igraph_vector_t *work, int root);
int _write_tree_newick(igraph_t *tree, char *out, int root, igraph_vector_t *vec);

igraph_t *parse_newick(FILE *f)
{
    igraph_t *tree = malloc(sizeof(igraph_t));
    char buf[BUFSIZ];
    char *cur;
    int *size;
    double *branch_length;
    int i, error, nnode = 0;
    igraph_vector_t *vec = malloc(sizeof(igraph_vector_t));

    // count the nodes in the tree
    while (fgets(buf, BUFSIZ, f) != NULL) 
    {
        cur = buf;
        while ((cur = strpbrk(cur, ",);")) != NULL)
        {
            ++nnode;
            ++cur;
        }
    }

    error = igraph_empty(tree, nnode, 1);
    igraph_vector_init(vec, 2*(nnode-1));
    branch_length = calloc(nnode, sizeof(double));
    size = calloc(nnode, sizeof(int));

    fseek(f, 0, SEEK_SET);
    yyrestart(f);
    yyparse(vec, size, branch_length);

    igraph_add_edges(tree, vec, 0);
    for (i = 0; i < nnode-1; ++i)
    {
        error = igraph_incident(tree, vec, i, IGRAPH_IN);
        SETEAN(tree, "length", (int) VECTOR(*vec)[0], branch_length[i]);
        SETVAN(tree, "extant", i, (double) size[i] == 1);
        SETVAN(tree, "size", i, (double) size[i]);
    }

    free(size);
    free(branch_length);
    igraph_vector_destroy(vec);
    return tree;
}

int root(const igraph_t *tree)
{
    long r = -1;
    igraph_vector_t work;

    igraph_vector_init(&work, 1);
    igraph_degree(tree, &work, igraph_vss_all(), IGRAPH_IN, 0);
    igraph_vector_search(&work, 0, 0, &r);
    igraph_vector_destroy(&work);
    return (int) r;
}

void ladderize(igraph_t *tree)
{
    igraph_vector_t work, permvec;
    igraph_t new_tree;
    int i, *perm = malloc(igraph_vcount(tree) * sizeof(int));

    igraph_vector_init(&work, 2);
    igraph_vector_init(&permvec, igraph_vcount(tree));

    _ladderize(tree, &work, root(tree), perm);
    for (i = 0; i < igraph_vcount(tree); ++i)
        igraph_vector_set(&permvec, perm[i], (double) i);
    igraph_permute_vertices(tree, &new_tree, &permvec);

    // this maybe is questionable
    memcpy(tree, &new_tree, sizeof(*tree));

    free(perm);
    igraph_vector_destroy(&work);
    igraph_vector_destroy(&permvec);
}

double height(const igraph_t *tree)
{
    igraph_vector_t work;
    igraph_vector_init(&work, 2);
    double ht = igraph_vcount(tree) > 1 ? _height(tree, &work, root(tree)) : 0.;
    igraph_vector_destroy(&work);
    return ht;
}

void write_tree_newick(igraph_t *tree, FILE *f)
{
    // TODO: should get a more accurate estimate of the space we need
    char *s = calloc(igraph_vcount(tree) * 100, sizeof(char));
    igraph_vector_t work;

    igraph_vector_init(&work, 0);
    _write_tree_newick(tree, s, root(tree), &work);
    fprintf(f, "%s", s);

    igraph_vector_destroy(&work);
    free(s);
}

int *production(const igraph_t *tree, igraph_vector_t *vec)
{
    int i, nnode = igraph_vcount(tree);
    int *p = malloc(nnode * sizeof(int));
    igraph_vector_t nbr;
    igraph_vector_init(&nbr, 2);

    igraph_degree(tree, vec, igraph_vss_all(), IGRAPH_OUT, 0);
    for (i = 0; i < nnode; ++i)
    {
        if ((int) VECTOR(*vec)[i] == 0)
        {
            p[i] = 0;
        }
        else
        {
            igraph_neighbors(tree, &nbr, i, IGRAPH_OUT);
            p[i] = ((int) VECTOR(*vec)[(int) VECTOR(nbr)[0]] == 0) +
                   ((int) VECTOR(*vec)[(int) VECTOR(nbr)[1]] == 0) + 1;
        }
    }
    return p;
}

double kernel(const igraph_t *t1, const igraph_t *t2, double lambda, double sigma, int coal)
{
    int i1 = 1, i2 = 1, start = 0, i, c1, c2, n1, n2;
    double val, tmp, K = 0;
    Pvoid_t delta = (Pvoid_t) NULL, pairs = (Pvoid_t) NULL;
    PWord_t Pvalue, cdelta;
    Word_t bytes = 0;
    int coord, nnode1 = igraph_vcount(t1), nnode2 = igraph_vcount(t2);
    int *production1, *production2;
    int *order1 = malloc(nnode1 * sizeof(int)), *order2 = malloc(nnode2 * sizeof(int));
    igraph_vector_t vec, v1, v2;

    // preconditions
    assert(lambda > 0.0 && lambda <= 1.0);
    assert(sigma > 0.0);

    igraph_vector_init(&vec, 2);
    igraph_vector_init(&v1, 2);
    igraph_vector_init(&v2, 2);

    production1 = production(t1, &vec); production2 = production(t2, &vec);
    order(production1, order1, nnode1); order(production2, order2, nnode2);

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
                {
                    JLI(Pvalue, pairs, coord);
                    if (Pvalue == PJERR) exit(EXIT_FAILURE); 
                    *Pvalue = coord;
                }

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
    
    JLF(Pvalue, pairs, bytes);
    while (Pvalue != NULL)
    {
        val = lambda;
        n1 = *Pvalue >> 16;
        n2 = *Pvalue & 65535;

        // branch lengths
        igraph_incident(t1, &v1, n1, IGRAPH_OUT);
        igraph_incident(t2, &v2, n2, IGRAPH_OUT);

        tmp = 0;
        for (i = 0; i < 2; ++i)
            tmp += pow(VECTOR(v1)[i] - VECTOR(v2)[i], 2);
        val *= gsl_ran_gaussian_pdf(sqrt(tmp), sigma);

        igraph_neighbors(t1, &v1, n1, IGRAPH_OUT);
        igraph_neighbors(t2, &v2, n2, IGRAPH_OUT);

        for (i = 0; i < 2; ++i)  // assume tree is binary
        {

            c1 = (int) VECTOR(v1)[i];
            c2 = (int) VECTOR(v2)[i];

            if (production1[c1] == production2[c2])
            {
                // children are leaves
                if (production1[c1] == 0)
                {
                    val *= (1 + lambda);
                }

                // children are not leaves
                else
                {
                    JLG(cdelta, delta, (c1 << 16) | c2);
                    /* don't visit children before parents */
                    assert(cdelta != NULL);
                    memcpy(&tmp, cdelta, sizeof(double));
                    val *= (1 + tmp);
                }
            }
            fflush(stdout);
        }

        JLI(Pvalue, delta, *Pvalue);
        if (Pvalue == PJERR) exit(EXIT_FAILURE); 
        memcpy(Pvalue, &val, sizeof(double));

        K += val;
        JLN(Pvalue, pairs, bytes);
    }

    free(production1);
    free(production2);
    free(order1);
    free(order2);
    JLFA(bytes, delta);
    JLFA(bytes, pairs);
    return K;
}

int _ladderize(igraph_t *tree, igraph_vector_t *work, int root, int *perm)
{
    igraph_es_t es;
    int do_swap, lc, rc, lsize, rsize;

    // if a leaf, just add the next number
    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_OUT, 0);
    if ((int) VECTOR(*work)[0] == 0)
    {
        perm[0] = root;
        return 1;
    }

    // get the children
    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    lc = (int) VECTOR(*work)[0];
    rc = (int) VECTOR(*work)[1];

    // order first by size
    lsize = _ladderize(tree, work, lc, perm);
    rsize = _ladderize(tree, work, rc, &perm[lsize]);
    do_swap = lsize > rsize;

    // order next by branch length
    if (!do_swap)
    {
        igraph_es_incident(&es, root, IGRAPH_OUT);
        igraph_cattribute_EANV(tree, "length", es, work);
        do_swap = (VECTOR(*work)[0] > VECTOR(*work)[1]);
    }

    // swap order of children if necessary
    if (do_swap)
        rotl(perm, (rsize + lsize) * sizeof(int), lsize * sizeof(int));

    // now add the root
    perm[lsize + rsize] = root;
    return lsize + rsize + 1;
}

double _height(const igraph_t *tree, igraph_vector_t *work, int root)
{
    int lc, rc;
    double height;
    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_OUT, 0);
    if ((int) VECTOR(*work)[0] == 0)
    {
        igraph_incident(tree, work, root, IGRAPH_IN);
        return EAN(tree, "length", (int) VECTOR(*work)[0]);
    }

    igraph_neighbors(tree, work, root, IGRAPH_OUT);
    lc = (int) VECTOR(*work)[0]; rc = (int) VECTOR(*work)[1];
    height = fmax(_height(tree, work, lc), _height(tree, work, rc));

    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_IN, 0);
    if ((int) VECTOR(*work)[0] == 0)
        return height;
    igraph_incident(tree, work, root, IGRAPH_IN);
    return EAN(tree, "length", (int) VECTOR(*work)[0]) + height;
}

int _write_tree_newick(igraph_t *tree, char *out, int root, 
        igraph_vector_t *work)
{
    double length;
    int nchar, left_child, right_child, is_root = 0;

    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_IN, 0);
    if ((int) VECTOR(*work)[0] == 0)
    {
        length = 0.;
        is_root = 1;
    }
    else
    {
        igraph_incident(tree, work, root, IGRAPH_IN);
        length = EAN(tree, "length", VECTOR(*work)[0]);
    }

    igraph_degree(tree, work, igraph_vss_1(root), IGRAPH_OUT, 0);
    if ((int) VECTOR(*work)[0] == 0)
        return sprintf(out, "%d:%f", root, length);

    igraph_neighbors(tree, work, root, IGRAPH_OUT);

    left_child = VECTOR(*work)[0];
    right_child = VECTOR(*work)[1];

    nchar = sprintf(out, "(");
    nchar += _write_tree_newick(tree, &out[nchar], left_child, work);
    nchar += sprintf(&out[nchar], ",");
    nchar += _write_tree_newick(tree, &out[nchar], right_child, work);
    nchar += sprintf(&out[nchar], ")");
    nchar += sprintf(&out[nchar], "%d:%f", root, length);
    if (is_root)
        nchar += sprintf(&out[nchar], ";");
    return nchar;
}
