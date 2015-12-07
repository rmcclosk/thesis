#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <igraph/igraph.h>
#include "util.h"

#define LOG_ZERO DBL_MIN_10_EXP / 2

gsl_rng *set_seed(int seed)
{
    gsl_rng *rng;
    if (seed > 0)
        seed = time(NULL);

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

int compare_ints (const void * a, const void * b)
{
    return ( *(int*)a - *(int*)b );
}

int compare_doubles (const void * a, const void * b)
{
    if ( *(double*)a < *(double*)b ) return -1;
    if ( *(double*)a > *(double*)b ) return 1;
    return 0;
}

void order(const void *base, int *order, size_t size, int nitems,
        int (*compar) (const void *, const void *))
{
    int i; 
    size_t pair_size = size + sizeof(int);
    char *buf = malloc(nitems * pair_size);

    for (i = 0; i < nitems; ++i) {
        memcpy(&buf[i * pair_size], &((char*) base)[i * size], size);
        memcpy(&buf[i * pair_size + size], &i, sizeof(int));
    }

    qsort(buf, nitems, pair_size, compar);
    for (i = 0; i < nitems; ++i)
        memcpy(&order[i], &buf[i * pair_size + size], sizeof(int));

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

double sum_doubles(double *x, int n)
{
    int i;
    double sum = 0;

    for (i = 0; i < n; ++i)
        sum += x[i];
    return sum;
}

double max_doubles(double *x, int n)
{
    int i;
    double max = x[0];

    for (i = 1; i < n; ++i)
        max = fmax(max, x[i]);
    return max;
}

void *safe_realloc(void *ptr, size_t size)
{
    void *tmp = realloc(ptr, size);
    if (tmp == NULL)
    {
        fprintf(stderr, "Aborting: out of memory\n");
        abort();
    }
    return tmp;
}

void permute(void *v, size_t size, int nitems, const int *perm, 
             void (get) (const void *, int, void *),
             void (set) (void *, int, const void *))
{
    int cyc_start = 0, cur, next;
    int *done = calloc(nitems, sizeof(int));
    char *cur_item = malloc(size), *next_item = malloc(size);

    while (cyc_start < nitems) {
        cur = cyc_start;
        next = perm[cyc_start];
        get(v, cur, cur_item);
        do {
            get(v, next, next_item);
            set(v, next, cur_item);
            memcpy(cur_item, next_item, size);
    
            done[next] = 1;
            cur = next;
            next = perm[next];
        } while (cur != cyc_start);

        while (done[cyc_start] == 1 && cyc_start < nitems) {
            ++cyc_start;
        }
    }

    free(done);
    free(cur_item);
    free(next_item);
}

void set_igraph_vector_t(void *v, int n, const void *value)
{
    igraph_vector_set((igraph_vector_t *) v, n, * (igraph_real_t *) value);
}

void get_igraph_vector_t(const void *v, int n, void *value)
{
    * (igraph_real_t *) value = VECTOR(*(igraph_vector_t *) v)[n];
}

void set_igraph_vector_bool_t(void *v, int n, const void *value)
{
    igraph_vector_bool_set((igraph_vector_bool_t *) v, n, 
                           * (igraph_bool_t *) value);
}

void get_igraph_vector_bool_t(const void *v, int n, void *value)
{
    * (igraph_bool_t *) value = VECTOR(*(igraph_vector_bool_t *) v)[n];
}

void set_igraph_strvector_t(void *v, int n, const void *value)
{
    igraph_strvector_set((igraph_strvector_t *) v, n, (const char *) value);
}

void get_igraph_strvector_t(const void *v, int n, void *value)
{
    strcpy(value, STR(* (igraph_strvector_t *) v, n));
}
