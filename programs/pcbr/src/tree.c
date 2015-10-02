#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tree.h"

int output_tree(const pcbr_workspace *w, const double *rates, int nnode, char *buf)
{
    int left_child = w->topology[(nnode-1)*2];
    int right_child = w->topology[(nnode-1)*2+1];
    int nchar;
    double rate = rates[w->A[nnode-1]] / w->branch_scale;
    double ll = log10(w->L[(nnode-1) * w->nrates + w->A[nnode-1]]) + \
                w->scale[nnode-1];
    double t = w->branch_lengths[nnode-1];
    int cluster = w->cluster[nnode-1];
    const char comment[] = "%d[&rate=%f,cluster=%d,log10lik=%f]:%f";

    if (left_child == NO_CHILD)
        return sprintf(buf, comment, nnode, rate, cluster, ll, t);

    nchar = sprintf(buf, "(");
    nchar += output_tree(w, rates, left_child+1, &buf[nchar]);
    nchar += sprintf(&buf[nchar], ",");
    nchar += output_tree(w, rates, right_child+1, &buf[nchar]);
    nchar += sprintf(&buf[nchar], ")");
    nchar += sprintf(&buf[nchar], comment, nnode, rate, cluster, ll, t);
    if (nnode == w->nnode)
        nchar += sprintf(&buf[nchar], ";");
    return nchar;
}
