#define NDEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <pll/pll.h>
#include <math.h>
#include <assert.h>
//#include <Rinternals.h>
//#include <Rembedded.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "likelihood.h"
//#include "rwrapper.h"
#include "util.h"

// helpers for dynamic programming
#define L(i, j) (w->L[ (i) * nrates + (j) ])
#define C(i, j) (w->C[ (i) * nrates + (j) ])
#define A(i) (w->A[ (i) ])
#define pi(i) (w->pi[ (i) ])
#define scale(i) (w->scale[ (i) ])
#define left_child(i) (w->topology[ (i) * 2 ])
#define right_child(i) (w->topology[ (i) * 2 + 1])
#define t(i) (w->t[ (i) ])

#define ROOT w->nnode - 1

#define USE_TIPS TRANS_ONLY
#define TRANS_PER_BRANCH ONE
#define CHEATING 0
#define TIP_PDF 0

void calculate_P(pcbr_workspace *w, double *lambda, trans_per_branch tpb, 
                 use_tips ut);

pcbr_workspace *pcbr_create(const pllNewickTree *tree, int nrates)
{
    pcbr_workspace *w = malloc(sizeof(pcbr_workspace));

    w->nnode = tree->nodes;
    w->nrates = nrates;

    w->topology = malloc(2 * tree->nodes * sizeof(int));
    w->branch_lengths = malloc(tree->nodes * sizeof(double));
    w->sorted_branch_lengths = malloc(tree->nodes * sizeof(double));
    flatten(tree->tree, tree->nodes, w->topology, w->branch_lengths, 0);

    w->t = malloc(tree->nodes * sizeof(double));
    w->brlen_order_inv = malloc(tree->nodes*sizeof(int));
    scale_branches(w);

    // dynamic programming matrices and scratch space
    w->L = malloc(nrates * sizeof(double) * tree->nodes);
    w->Lz = malloc(nrates * sizeof(double));
    w->C = malloc(nrates * sizeof(int) * tree->nodes);
    w->A = malloc(sizeof(int) * tree->nodes);
    w->scale = malloc(sizeof(int) * tree->nodes);
    w->cluster = malloc(tree->nodes * sizeof(int));

    // space for parameters
    w->pi = malloc(nrates*sizeof(double));
    w->Q = gsl_matrix_alloc(nrates, nrates);
    w->P = malloc(nrates*nrates*tree->nodes*sizeof(double));
    w->Qbar = malloc(nrates * sizeof(double));

    // GSL stuff
    w->Qt = gsl_matrix_alloc(nrates, nrates);
    w->expQt = gsl_matrix_alloc(nrates, nrates);
    w->eval = gsl_vector_complex_alloc(nrates);
    w->evec = gsl_matrix_complex_alloc(nrates, nrates);
    w->w = gsl_eigen_nonsymmv_alloc(nrates);

    return w;
}

void pcbr_free(pcbr_workspace *p)
{
    free(p->L);
    free(p->Lz);
    free(p->C);
    free(p->A);
    free(p->scale);
    free(p->pi);
    free(p->P);
    free(p->topology);
    free(p->branch_lengths);
    free(p->sorted_branch_lengths);
    free(p->brlen_order_inv);
    free(p->t);
    free(p->cluster);
    free(p->Qbar);
    gsl_matrix_free(p->Q);
    gsl_matrix_free(p->Qt);
    gsl_matrix_free(p->expQt);
    gsl_vector_complex_free(p->eval);
    gsl_matrix_complex_free(p->evec);
    gsl_eigen_nonsymmv_free(p->w);
    free(p);
}

void scale_branches(pcbr_workspace *w)
{
    int i;
    int *brlen_order = malloc(w->nnode*sizeof(int)); 

    // scale the branch lengths to be in [0, 1]
    w->branch_scale = 0;
    for (i = 0; i < w->nnode; ++i)
        w->branch_scale = fmax(w->branch_scale, w->branch_lengths[i]);
    for (i = 0; i < w->nnode; ++i)
        w->t[i] = w->branch_lengths[i] / w->branch_scale;

    // store the order of branch lengths
    order(w->t, w->sorted_branch_lengths, brlen_order, w->nnode);
    inverse_perm(brlen_order, w->brlen_order_inv, w->nnode);
    free(brlen_order);
}


int flatten(const pllStack *tree, int nnode, int *order, double *bl, int next)
{
    pllNewickNodeInfo *info = (pllNewickNodeInfo *) tree->item;
    int size;

    if (next > 0)
        return flatten(tree->next, nnode-1, order, bl, next-1) + 1;

    bl[nnode-1] = atof(info->branch);
    if (info->leaf == PLL_TRUE) {
        order[(nnode-1)*2] = order[(nnode-1)*2+1] = NO_CHILD;
        return 1;
    }

    size = flatten(tree->next, nnode-1, order, bl, 0);
    order[(nnode-1)*2] = nnode - 2 - size;
    order[(nnode-1)*2+1] = nnode - 2;

    return flatten(tree->next, nnode-1, order, bl, size) + 1;
}

void guess_parameters(pcbr_workspace *w, double *args)
{
    int nbranch = (w->nnode - 3)/2, nrates = w->nrates;
    double intsum = 0;
    double *int_edges = malloc(nbranch * sizeof(double));
    double *guess = malloc(nrates * sizeof(double));
    int i, intcur = 0;

    //R_library("mixdist");

    // collect the branch lengths
    int_edges = malloc(nbranch * sizeof(double));
    for (i = 0; i < ROOT; ++i) {
        if (left_child(i) != NO_CHILD) {
            int_edges[intcur++] = w->t[i];
            intsum += w->t[i];
        }
    }

    // guess initial rates for mixed distribution
    guess_means(int_edges, w->nrates, intcur, guess);
    for (i = 0; i < w->nrates; ++i)
        guess[i] = fmax(guess[i], 1e-4);

    // branching rates are the inverses of the fitted means
    for (i = 0; i < nrates; ++i)
        args[i] = 1.0/guess[nrates - i - 1];

    // choose transition rates all equal such that the expected number of
    // transitions is 1% of the number of branches
    // TODO: make this informed by the equilibrium frequencies
    for (i = 0; i < nrates * (nrates - 1); ++i)
        args[nrates+i] = nbranch / w->nrates / intsum / 100;

    free(int_edges);
    free(guess);
}

void guess_means(const double *x, int nrates, int nbranch, double *guess)
{
    int i, part_size = (int) ceil((float) nbranch/ (float) nrates);
    double *copy = malloc(nbranch*sizeof(double));

    // make sure branch lengths have been scaled properly
    for (i = 0; i < nbranch; ++i)
        assert(x[i] >= 0 && x[i] <= 1);

    memcpy(copy, x, nbranch*sizeof(double));
    qsort(copy, nbranch, sizeof(double), compare_doubles);

    for (i = 0; i < nrates; ++i)
        guess[i] = 0;
    for (i = 0; i < nbranch; ++i) 
        guess[i / part_size] += copy[i] / part_size;

    free(copy);
}

void fill_parameters(pcbr_workspace *w, double *q)
{
    int i, j, cur = -1;
    double sum;

    // calculate Q
    for (i = 0; i < w->nrates; ++i) {
        sum = 0;
        for (j = 0; j < w->nrates; ++j) {
            if (i != j) {
                gsl_matrix_set(w->Q, i, j, q[++cur]);
                sum += q[cur];
            }
        }
        gsl_matrix_set(w->Q, i, i, -sum);
        w->Qbar[i] = sum;
    }

    // calculate pi
    gsl_matrix_memcpy(w->Qt, w->Q);
    gsl_matrix_transpose(w->Qt);
    gsl_eigen_nonsymmv(w->Qt, w->eval, w->evec, w->w);
    gsl_eigen_nonsymmv_sort(w->eval, w->evec, GSL_EIGEN_SORT_ABS_ASC);
    // TODO: should exit with error if the model isn't stationary

    for (i = 0; i < w->nrates; ++i)
        w->pi[i] = GSL_REAL(gsl_matrix_complex_get(w->evec, i, 0));
    sum = double_sum(w->pi, w->nrates);
    for (i = 0; i < w->nrates; ++i)
        w->pi[i] /= sum;
}

void calculate_P(pcbr_workspace *w, double *lambda, 
        trans_per_branch tpb, use_tips ut)
{
    int i, j, node, nrates = w->nrates, terminal;
    double t = 0.0;
    double *P;
    gsl_matrix *Q = w->Q;
    double *Qbar = w->Qbar;

    for (node = 0; node < w->nnode; ++node) {
        t = w->t[node];
        P = &w->P[node * nrates * nrates];
        terminal = (left_child(node) == NO_CHILD);

        if (terminal && ut == NO) {
            memset(P, 0, nrates * nrates * sizeof(double));
            for (i = 0; i < nrates; ++i)
                P[i * nrates + i] = 1;
        }

        else if (tpb == ONE) {
            for (i = 0; i < nrates; ++i) {
                for (j = 0; j < nrates; ++j) {
                    if (i != j) {
                        if ((!terminal) || (ut == YES)) {
                            P[i * nrates + j] = gsl_matrix_get(Q, i, j) / 
                                (Qbar[i] - Qbar[j] + lambda[i] - lambda[j]) *
                                (exp(-(Qbar[j] + lambda[j]) * t) - 
                                 exp(-(Qbar[i] + lambda[i]) * t));
                        } else if (ut == TRANS_ONLY) {
                            P[i * nrates + j] = gsl_matrix_get(Q, i, j) /
                                (Qbar[i] - Qbar[j]) *
                                (exp(-Qbar[j] * t) - exp(-Qbar[i] * t));
                        } 
                    } else { /* i == j */
                        if ((!terminal) || (ut == YES)) {
                            P[i * nrates + j] = exp(-(Qbar[i] + lambda[i]) * t);
                        } else if (ut == TRANS_ONLY) {
                            P[i * nrates + j] = exp(-Qbar[i] * t);
                        } 
                    }
                    if (P[i * nrates + j] != P[i * nrates + j])
                        P[i * nrates + j] = 0;
                }
            }
        } 

        else if (tpb == ANY) {
            gsl_matrix_memcpy(w->Qt, Q);
            if ((!terminal) || (ut == YES)) {
                for (i = 0; i < nrates; ++i)
                    gsl_matrix_set(w->Qt, i, i, -Qbar[i] - lambda[i]);
            }
            gsl_matrix_scale(w->Qt, t);
            gsl_linalg_exponential_ss(w->Qt, w->expQt, GSL_PREC_DOUBLE);
            memcpy(P, w->expQt->data, nrates * nrates * sizeof(double));
        }
    }
}

double likelihood(pcbr_workspace *w, double *args, int reconstruct, use_tips ut,
                  trans_per_branch tpb, int tip_pdf, int cheating)
{
    double *Lz = w->Lz;
    int nrates = w->nrates;
    int i, left_child, right_child, parent_state, child_state;
    double lik;
    int new_scale;
    double *P;

    fill_parameters(w, &args[nrates]);
    calculate_P(w, args, tpb, ut);
    
    for (i = 0; i <= ROOT; ++i) {
        left_child = left_child(i);
        right_child = right_child(i);
        P = &w->P[i * nrates * nrates];

        if (left_child == NO_CHILD)
            scale(i) = 0;
        else
            scale(i) = scale(left_child) + scale(right_child);

        for (parent_state = 0; parent_state < nrates; ++parent_state) {
            if (i == ROOT) {
                L(i, parent_state) = pi(parent_state) * 
                                     L(right_child, parent_state) *
                                     L(left_child,  parent_state);
                C(i, parent_state) = parent_state;
            } else {
                for (child_state = 0; child_state < nrates; ++child_state) {
                    Lz[child_state] = P[parent_state * nrates + child_state];
                    if (left_child != NO_CHILD)   // not a tip
                        Lz[child_state] *= L(left_child, child_state) * 
                                           L(right_child, child_state) *
                                           args[child_state];
                    else if (tip_pdf)
                        Lz[child_state] *= args[child_state];
                }

                if (reconstruct) {
                    C(i, parent_state) = which_max(Lz, nrates);
                    L(i, parent_state) = Lz[C(i, parent_state)];

                    if (cheating && left_child == NO_CHILD)
                        C(i, parent_state) = parent_state;
                } else {
                    L(i, parent_state) = double_sum(Lz, nrates);
                }
            }
        }

        new_scale = get_scale(&L(i, 0), nrates);
        for (parent_state = 0; parent_state < nrates; ++parent_state)
            L(i, parent_state) /= pow(10, new_scale);
        scale(i) += new_scale;
    }

    for (parent_state = 0; parent_state < nrates; ++parent_state)
        Lz[parent_state] = L(ROOT, parent_state);

    if (reconstruct) {
        A(ROOT) = which_max(Lz, nrates);
        for (i = ROOT; i >= 0; --i) {
            if (left_child(i) != NO_CHILD) {
                A(left_child(i)) = C(left_child(i), A(i));
                A(right_child(i)) = C(right_child(i), A(i));
            }
        }
    }

    lik = (reconstruct ? double_max : double_sum)(Lz, nrates);
    /*
    for (i = 0; i < nrates * nrates; ++i)
        fprintf(stderr, "%f\t", args[i]);
    fprintf(stderr, "%f\n", log10(lik) + scale(ROOT));
    */
    return log10(lik) + scale(ROOT);
}

void get_clusters (pcbr_workspace *w, double *rates)
{
    int i, lc, rc;
    int ccount = 0;
    double minrate = double_min(rates, w->nrates);

    memset(w->cluster, 0, w->nnode * sizeof(int));

    if (rates[A(ROOT)] != minrate)
        w->cluster[ROOT] = ++ccount;

    for (i = w->nnode-1; i >= 0; --i) {
        lc = left_child(i);
        rc = right_child(i);
        if (lc != NO_CHILD) {
            if (rates[A(lc)] == minrate)
                w->cluster[lc] = 0;
            else if (rates[A(lc)] != rates[A(i)])
                w->cluster[lc] = ++ccount;
            else
                w->cluster[lc] = w->cluster[i];

            if (rates[A(rc)] == minrate)
                w->cluster[rc] = 0;
            else if (rates[A(lc)] != rates[A(i)])
                w->cluster[rc] = ++ccount;
            else
                w->cluster[rc] = w->cluster[i];
        }
    }
}

double lrt(double log10lik_null, double log10lik_alt, int nparam_null, int nparam_alt)
{
    double teststat = - 2 * log10lik_null / log10(exp(1)) + 2 * log10lik_alt / log10(exp(1));
    int df = nparam_alt - nparam_null;
    return 1 - gsl_cdf_chisq_P(teststat, df);
}

double bic(double log10lik, int nparam, int ndata)
{
    return -2 * log10lik / log10(exp(1)) + nparam * log(ndata);
}
