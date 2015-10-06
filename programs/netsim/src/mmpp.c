#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <igraph/igraph.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_cblas.h>
#include "mmpp.h"
#include "tree.h"
#include "util.h"

struct ck_params {
    int nrates;
    double *Q_minus_Lambda;
};
int ckdiff(double t, const double y[], double f[], void *params);
void calculate_P(const igraph_t *tree, int nrates, const double *theta, double *P);
void calculate_pi(const igraph_t *tree, int nrates, const double *theta, double *pi);

void guess_parameters(const igraph_t *tree, int nrates, double *theta)
{
    int i, nbranch = (igraph_vcount(tree) - 3)/2, cur = 0;
    int part_size = (int) ceil((float) nbranch / (float) nrates);
    igraph_inclist_t il;
    igraph_vector_t degree;
    igraph_vector_int_t *edge;
    double *int_edges = malloc(nbranch * sizeof(double));
    double intsum = 0;

    // find out-degree of all nodes
    igraph_vector_init(&degree, igraph_vcount(tree));
    igraph_degree(tree, &degree, igraph_vss_all(), IGRAPH_OUT, 0);

    // collect the internal branch lengths
    igraph_inclist_init(tree, &il, IGRAPH_IN);
    for (i = 0; i < igraph_vcount(tree); ++i) {
        edge = igraph_inclist_get(&il, i);
        if (igraph_vector_int_size(edge) > 0 && (int) VECTOR(degree)[i] > 0) {
            int_edges[cur++] = EAN(tree, "length", VECTOR(*edge)[0]);
            intsum += int_edges[cur-1];
        }
    }

    // for branching rates, partition internal branch lengths, and take means of partitions
    qsort(int_edges, nbranch, sizeof(double), compare_doubles);
    memset(theta, 0, nrates * sizeof(double));
    for (i = 0; i < nbranch; ++i)
        theta[nrates - (i / part_size) - 1] += int_edges[i] / part_size;

    // the guessed rates are the inverse of the means
    for (i = 0; i < nrates; ++i)
        theta[i] = 1.0 / fmax(theta[i], 1e-4);

    // choose transition rates all equal such that the expected number of
    // transitions is 1% of the number of branches
    for (i = 0; i < nrates * (nrates - 1); ++i)
        theta[nrates+i] = nbranch / nrates / intsum / 100.0;

    // clean up
    igraph_inclist_destroy(&il);
    free(int_edges);
}

double likelihood(const igraph_t *tree, int nrates, const double *theta)
{ 
    // TODO: some of this stuff (malloc'ing, the adjacency list, ...) can be
    // done once beforehand
    int i, rt = root(tree);
    double *P = malloc(nrates * nrates * igraph_vcount(tree) * sizeof(double));
    double *L = malloc(nrates * igraph_vcount(tree) * sizeof(double));
    double *pi = malloc(nrates * sizeof(double));
    double lik = 0;
    int *scale = malloc(igraph_vcount(tree) * sizeof(int));
    int lchild, rchild, pstate, cstate, new_scale;
    igraph_adjlist_t al;
    igraph_vector_int_t *children;

    igraph_adjlist_init(tree, &al, IGRAPH_OUT);
    calculate_P(tree, nrates, theta, P);
    calculate_pi(tree, nrates, theta, pi);

    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        children = igraph_adjlist_get(&al, i);
        if (igraph_vector_int_size(children) == 0)
        {
            lchild = -1; rchild = -1;
            scale[i] = 0;
        }
        else
        {
            lchild = VECTOR(*children)[0];
            rchild = VECTOR(*children)[1];
            scale[i] = scale[lchild] + scale[rchild];
        }

        for (pstate = 0; pstate < nrates; ++pstate)
        {
            if (i == rt) {
                L[i * nrates + pstate] = pi[pstate] *
                                         L[rchild * nrates + pstate] *
                                         L[lchild * nrates + pstate];
            } 
            else {
                L[i * nrates + pstate] = 0;
                for (cstate = 0; cstate < nrates; ++cstate) {
                    if (lchild == -1) {
                        L[i * nrates + pstate] += P[pstate * nrates + cstate];
                    }
                    else {
                        L[i * nrates + pstate] += P[pstate * nrates + cstate] *
                                                  L[lchild * nrates + cstate] *
                                                  L[rchild * nrates + cstate] *
                                                  theta[cstate];
                    }
                }
            }
        }
        new_scale = get_scale(&L[i * nrates], nrates);
        for (pstate = 0; pstate < nrates; ++pstate)
        {
            L[i * nrates + pstate] /= pow(10, new_scale);
        }
        scale[i] += new_scale;
    }

    for (pstate = 0; pstate < nrates; ++pstate)
    {
        lik += L[rt * nrates + pstate];
    }

    igraph_adjlist_destroy(&al);
    free(L);
    free(P);
    free(scale);
    free(pi);
    return log10(lik) + scale[rt];
}

/* Private. */

void calculate_P(const igraph_t *tree, int nrates, const double *theta, double *P)
{
    int i, j, from, to, cur = nrates; 
    int *bl_order = malloc((igraph_vcount(tree)-1) * sizeof(int));
    double rowsum, t = 0.0;
    double y[nrates * nrates];
    struct ck_params ckpar;
    gsl_odeiv2_system sys;
    gsl_odeiv2_driver *d;
    igraph_vector_t branch_lengths;

    // set up ODE
    ckpar.nrates = nrates;
    ckpar.Q_minus_Lambda = malloc(nrates * nrates * sizeof(double));
    for (i = 0; i < nrates; ++i)
    {
        rowsum = 0;
        for (j = 0; j < nrates; ++j)
        {
            if (j != i)
            {
                y[i * nrates + j] = 0;
                ckpar.Q_minus_Lambda[i * nrates + j] = theta[cur];
                rowsum += theta[cur++];
            }
        }
        ckpar.Q_minus_Lambda[i * nrates + i] = -rowsum - theta[i];
        y[i * nrates + i] = 1;
    }
    sys.function = ckdiff;
    sys.jacobian = NULL;
    sys.dimension = nrates * nrates;
    sys.params = &ckpar;
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    // collect and order branch lengths
    igraph_vector_init(&branch_lengths, igraph_vcount(tree)-1);
    EANV(tree, "length", &branch_lengths);

    order(VECTOR(branch_lengths), bl_order, sizeof(double), igraph_ecount(tree), compare_doubles);
    for (i = 0; i < igraph_ecount(tree); ++i)
        fprintf(stderr, "%d\t%f\n", bl_order[i], VECTOR(branch_lengths)[bl_order[i]]);
    fprintf(stderr, "\n");

    // set the values at the root of the tree to the identity
    memcpy(&P[root(tree) * nrates * nrates], y, nrates * nrates * sizeof(double));

    // calculate values for each branch
    for (i = 0; i < igraph_ecount(tree); ++i)
    {
        fprintf(stderr, "%d\t%f\n", bl_order[i], VECTOR(branch_lengths)[bl_order[i]]);
        if (VECTOR(branch_lengths)[bl_order[i]] != t)
            gsl_odeiv2_driver_apply(d, &t, VECTOR(branch_lengths)[bl_order[i]], y);
        igraph_edge(tree, bl_order[i], &from, &to);
        memcpy(&P[to * nrates * nrates], y, nrates * nrates * sizeof(double));
    }

    // clean up
    free(bl_order);
    free(ckpar.Q_minus_Lambda);
    igraph_vector_destroy(&branch_lengths);
    gsl_odeiv2_driver_free(d);
}

void calculate_pi(const igraph_t *tree, int nrates, const double *theta, double *pi)
{
    // TODO: malloc beforehand
    gsl_matrix *Q = gsl_matrix_alloc(nrates, nrates);
    int i, j, cur = nrates;
    double sum;
    gsl_vector_complex *eval = gsl_vector_complex_alloc(nrates);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc(nrates, nrates);
    gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(nrates);

    // calculate Q
    for (i = 0; i < nrates; ++i) {
        sum = 0;
        for (j = 0; j < nrates; ++j) {
            if (i != j) {
                gsl_matrix_set(Q, i, j, theta[cur]);
                sum += theta[cur++];
            }
        }
        gsl_matrix_set(Q, i, i, -sum);
    }

    // calculate pi
    gsl_matrix_transpose(Q);
    gsl_eigen_nonsymmv(Q, eval, evec, w);
    gsl_eigen_nonsymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

    sum = 0;
    for (i = 0; i < nrates; ++i)
    {
        pi[i] = GSL_REAL(gsl_matrix_complex_get(evec, i, 0));
        sum += pi[i];
    }
    for (i = 0; i < nrates; ++i)
        pi[i] /= sum;

    gsl_matrix_free(Q);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
    gsl_eigen_nonsymmv_free(w);
}

int ckdiff(double t, const double y[], double f[], void *params)
{
    struct ck_params *ckpar = (struct ck_params *) params;
    int n = ckpar->nrates;
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, y, n, 
                ckpar->Q_minus_Lambda, n, 0.0, f, n);
    return GSL_SUCCESS;
}
