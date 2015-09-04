#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <igraph/igraph.h>
#include "util.h"

int compare_pairs(const void *a, const void *b);

gsl_rng *set_seed(int seed)
{
    gsl_rng *rng;

    // C
    srand(seed);

    // GSL
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, seed);

    // igraph
    igraph_rng_seed(igraph_rng_default(), seed);

    // GSL requires an object, so return it
    return rng;
}

void order(const int *x, int *order, int n)
{
    int *buf = malloc(n*2*sizeof(int));
    int i;

    for (i = 0; i < n; ++i) {
        buf[2*i] = x[i];
        buf[2*i+1] = i;
    }
    qsort(buf, n, 2*sizeof(int), compare_pairs);
    for (i = 0; i < n; ++i)
        order[i] = buf[2*i+1];

    free(buf);
}

void rotl(void *x, size_t nx, size_t n)
{
    void *tmp = malloc(n);
    memcpy(tmp, x, n); // copy first n bytes of x
    memmove(x, x + n, nx - n); // shift x left by nx - n bytes
    memcpy(x + nx - n, tmp, n); // copy first n bytes into last n bytes
    free(tmp);
}

/* Private */

int compare_pairs(const void *a, const void *b)
{
    int *x = (int*) a;
    int *y = (int*) b;

    if (x[0] < y[0])
        return -1;
    else if (x[0] > y[0])
        return 1;
    else
        return x[1] < y[1] ? -1 : 1;
}
