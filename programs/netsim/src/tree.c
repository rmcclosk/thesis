#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <igraph/igraph.h>

#include "newick_parser.h"
#include "tree.h"
#include "util.h"

#define NDEBUG

void yyrestart(FILE *f);

/* recursive functions */
int _ladderize(igraph_t *tree, igraph_vector_t *work, int root, int *perm);
double _height(const igraph_t *tree, igraph_vector_t *work, int root);
int _write_tree_newick(const igraph_t *tree, char *out, int root, igraph_vector_t *vec);
void _cut_at_time(igraph_t *tree, double t, int root, double troot, 
        int extant_only, igraph_vector_t *work, igraph_vector_t *to_delete);
int _collapse_singles(igraph_t *tree, int root, igraph_vector_t *vdel,
        igraph_vector_t *eadd, igraph_vector_t *branch_length,
        igraph_vector_t *work, double *bl);

igraph_t *parse_newick(FILE *f)
{
    igraph_t *tree;
    int i;
    extern int yynode;
    igraph_vector_t edge, branch_length, size;

    igraph_vector_init(&edge, 0);
    igraph_vector_init(&size, 0);
    igraph_vector_init(&branch_length, 0);

    yynode = 0;
    yyrestart(f);
    yyparse(&edge, &size, &branch_length);

    tree = malloc(sizeof(igraph_t));
    igraph_empty(tree, igraph_vector_size(&size), 1);
    igraph_add_edges(tree, &edge, 0);

    for (i = 0; i < igraph_vector_size(&size); ++i)
    {
        igraph_incident(tree, &edge, i, IGRAPH_IN);
        if (igraph_vector_size(&edge) > 0) {
            SETEAN(tree, "length", (int) VECTOR(edge)[0], VECTOR(branch_length)[i]);
        }
    }

    igraph_vector_destroy(&edge);
    igraph_vector_destroy(&size);
    igraph_vector_destroy(&branch_length);
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

void write_tree_newick(const igraph_t *tree, FILE *f)
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

void scale_branches(igraph_t *tree, scaling mode)
{
    int i;
    double scale;
    igraph_vector_t vec;
    igraph_vector_init(&vec, igraph_vcount(tree));

    EANV(tree, "length", &vec);
    switch (mode)
    {
        case MEAN:
            scale = gsl_stats_mean(VECTOR(vec), 1, igraph_ecount(tree));
            break;
        case MEDIAN:
            igraph_vector_sort(&vec);
            scale = gsl_stats_median_from_sorted_data(VECTOR(vec), 1, igraph_ecount(tree));
            break;
        default:
            scale = 1.;
            break;
    }

    for (i = 0; i < igraph_ecount(tree); ++i) {
        SETEAN(tree, "length", i, EAN(tree, "length", i) / scale);
    }
    igraph_vector_destroy(&vec);
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
    int i, j, orig_ntip = (igraph_vcount(tree) + 1)/2;
    igraph_vector_t tips, keep_tips, keep_all, drop, degree;
    igraph_vector_ptr_t nbhd;
    igraph_vector_t *elem;

    if (orig_ntip <= ntip)
        return;

    igraph_vector_init(&tips, 0);
    igraph_vector_init(&keep_tips, ntip);
    igraph_vector_init(&keep_all, 0);
    igraph_vector_init(&degree, 0);
    igraph_vector_init(&drop, 0);
    igraph_vector_ptr_init(&nbhd, 0);

    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        if ((int) VECTOR(degree)[i] == 0) {
            igraph_vector_push_back(&tips, i);
        }
    }
    gsl_ran_choose(rng, VECTOR(keep_tips), ntip, VECTOR(tips), orig_ntip, sizeof(igraph_real_t));

    igraph_neighborhood(tree, &nbhd, igraph_vss_vector(&keep_tips), INT_MAX, IGRAPH_IN, 0);

    for (i = 0; i < igraph_vector_ptr_size(&nbhd); ++i) {
        elem = (igraph_vector_t *) igraph_vector_ptr_e(&nbhd, i);
        for (j = 0; j < igraph_vector_size(elem); ++j) {
            igraph_vector_push_back(&keep_all, VECTOR(*elem)[j]);
        }
    }

    igraph_vector_sort(&keep_all);
    for (i = 0; i < igraph_vcount(tree); ++i) {
        if (!igraph_vector_binsearch2(&keep_all, i))
            igraph_vector_push_back(&drop, i);
    }

    igraph_delete_vertices(tree, igraph_vss_vector(&drop));
    collapse_singles(tree);

    igraph_vector_destroy(&tips);
    igraph_vector_destroy(&keep_tips);
    igraph_vector_destroy(&keep_all);
    igraph_vector_destroy(&degree);
    igraph_vector_destroy(&drop);
    igraph_vector_ptr_destroy_all(&nbhd);
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

int _write_tree_newick(const igraph_t *tree, char *out, int root,
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
