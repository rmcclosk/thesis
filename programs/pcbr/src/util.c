#include <math.h>
#include <stdlib.h>
#include "util.h"

void swap(double *x, int i, int j)
{
    double tmp = x[i];

    x[i] = x[j];
    x[j] = tmp;
}

int which_max(double *x, int n)
{
    double max = x[0];
    int which_max = 0;
    int i;

    for (i = 1; i < n; ++i) {
        if (max < x[i]) {
            max = x[i];
            which_max = i;
        }
    }
    return which_max;
}

double double_sum(double *x, int n)
{
    int i;
    double sum = 0;

    for (i = 0; i < n; ++i)
        sum += x[i];
    return sum;
}

double double_max(double *x, int n)
{
    int i;
    double max = x[0];

    for (i = 1; i < n; ++i)
        max = fmax(max, x[i]);
    return max;
}

double double_min(double *x, int n)
{
    int i;
    double min = x[0];

    for (i = 1; i < n; ++i)
        min = fmin(min, x[i]);
    return min;
}

void double_exp(double *x, double *result, int n)
{
    int i;

    for (i = 0; i < n; ++i)
        result[i] = exp(x[i]);
}

void double_log(double *x, double *result, int n)
{
    int i;

    for (i = 0; i < n; ++i)
        result[i] = log(x[i]);
}

int get_scale(double *x, int n)
{
    int i, nkeep = 0;
    double ilog, logsum = 0, logmax = LOG_ZERO;

    for (i = 0; i < n; ++i)
        logmax = fmax(x[i] == 0 ? LOG_ZERO : log10(x[i]), logmax);

    for (i = 0; i < n; ++i) {
        if (x[i] > 0) {
            ilog = log10(x[i]);
            if (ilog - logmax > LOG_ZERO) {
                logsum += ilog;
                nkeep += 1;
            }
        }
    }
    if (nkeep == 0)
        return 0;
    return (int) ceil(logsum / nkeep);
}

int compare_doubles(const void *a, const void *b)
{
    return ( *(double*) a > *(double*) b ? 1 : -1);
}

void order(const double *x, double *sorted, int *order, int n)
{
    // TODO: malloc beforehand only once and pass a buffer
    // (for now it's fine since it's only called once when the tree is initialized)
    double *buf = malloc(n*2*sizeof(double));
    int i;

    for (i = 0; i < n; ++i) {
        buf[2*i] = x[i];
        buf[2*i+1] = (double) i;
    }
    qsort(buf, n, 2*sizeof(double), compare_doubles);
    for (i = 0; i < n; ++i) {
        sorted[i] = buf[2*i];
        order[(int) buf[2*i+1]] = i;
    }

    free(buf);
}

void inverse_perm(const int *p, int *pinv, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        pinv[p[i]] = i;
}
