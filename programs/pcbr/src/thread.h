#ifndef THREAD_H
#define THREAD_H

#include <pll/pll.h>
#include "config.h"
#include "likelihood.h"

typedef struct {
    pcbr_workspace *w;
    double **args;
    double *arFunvals;
    int neval;
    use_tips ut;
    trans_per_branch tpb;
    int tip_pdf;
    int cheating;
} thread_data;

thread_data *create_thread_data(const pllNewickTree *tree, int nrate, 
        int nthread, int neval, use_tips ut, trans_per_branch tpb, int tip_pdf,
        int cheating);

void set_args(thread_data *tdata, double **args, int nthread, int narg);

void *do_likelihood(void *tdata);

void destroy_thread_data(thread_data *tdata, int nthread);

#endif
