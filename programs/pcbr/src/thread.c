#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pll/pll.h>

#include "thread.h"
#include "likelihood.h"

/*
typedef struct {
    pcbr_workspace *w;
    double **args;
    double *arFunvals;
    int neval;
} thread_data;
*/

thread_data *create_thread_data(const pllNewickTree *tree, int nrate, int
        nthread, int neval, use_tips ut, trans_per_branch tpb, int tip_pdf, 
        int cheating)
{
    int i, j;
    int narg = nrate * nrate;
    thread_data *tdata = malloc(nthread * sizeof(thread_data));

    for (i = 0; i < nthread; ++i) {
        tdata[i].w = pcbr_create(tree, nrate);
        tdata[i].neval = neval / nthread + (i == 0) * (neval % nthread);

        tdata[i].args = malloc(tdata[i].neval * sizeof(double*));
        for (j = 0; j < tdata[i].neval; ++j)
            tdata[i].args[j] = malloc(narg * sizeof(double));

        tdata[i].arFunvals = malloc(tdata[i].neval * sizeof(double));

        tdata[i].ut = ut;
        tdata[i].tpb = tpb;
        tdata[i].tip_pdf = tip_pdf;
        tdata[i].cheating = cheating;
    }

    return tdata;
}

void destroy_thread_data(thread_data *tdata, int nthread)
{
    int i, j;
    for (i = 0; i < nthread; ++i) {
        pcbr_free(tdata[i].w);
        for (j = 0; j < tdata[i].neval; ++j)
            free(tdata[i].args[j]);
        free(tdata[i].args);
        free(tdata[i].arFunvals);
    }
    free(tdata);
}

void set_args(thread_data *tdata, double **args, int nthread, int narg)
{
    int i, j, cur = 0;

    for (i = 0; i < nthread; ++i) {
        for (j = 0; j < tdata[i].neval; ++j)
            memcpy(tdata[i].args[j], args[cur++], narg * sizeof(double));
    }
}

void *do_likelihood(void *tdata)
{
    int i;
    thread_data *t = (thread_data *) tdata;
    for (i = 0; i < t->neval; ++i)
        t->arFunvals[i] = -likelihood(t->w, t->args[i], 0, t->ut, t->tpb,
                t->tip_pdf, t->cheating);
    return 0;
}
