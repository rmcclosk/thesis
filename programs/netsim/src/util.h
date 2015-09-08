/** \file util.h
 * \brief Utility functions. 
 */

#ifndef UTIL_H
#define UTIL_H

#include <gsl/gsl_rng.h>

/** Set seeds for all random number generators.
 *
 * This sets the seeds of the random number generators used by libc, GSL, and
 * igraph all to the same value. Since GSL's random generation requires that
 * you pass the generator as a parameter, the generator is returned.
 *
 * \param[in] seed the new seed
 * \return the GSL random generator
 */
gsl_rng *set_seed(int seed);

/** Compare two doubles.
 *
 * This is intended for use as a comparator in qsort and order.
 * 
 * \param[in] a,b doubles to compare
 * \return an integer with the same sign as a - b, or zero if a == b
 * \sa order()
 * \sa http://man7.org/linux/man-pages/man3/qsort.3.html
 */
int compare_doubles (const void * a, const void * b);

/** Compare two integers.
 *
 * This is intended for use as a comparator in qsort and order.
 * 
 * \param[in] a,b integers to compare
 * \return an integer with the same sign as a - b, or zero if a == b
 * \sa order()
 * \sa http://man7.org/linux/man-pages/man3/qsort.3.html
 */
int compare_ints (const void * a, const void * b);

/** Get the sorting order of an array.
 *
 * This function returns a permutation of x which would put x in sorted order.
 * It's like the order function in R, only zero-indexed. The compar function is
 * like one which is used for qsort.
 *
 * \param[in] base an array of items to order
 * \param[out] order the result will be stored here
 * \param[in] nitems number of items to order
 * \param[in] compar function which compares two elements of base
 *
 * \sa https://stat.ethz.ch/R-manual/R-devel/library/base/html/order.html
 * \sa http://man7.org/linux/man-pages/man3/qsort.3.html
 */
void order(const void *base, int *order, size_t size, int nitems,
        int (*compar) (const void *, const void *));

/** Circular left shift a block of memory.
 *
 * \param[in] x the block of memory to shift
 * \param[in] nx the size of x, in bytes
 * \param[in] n the number of bytes to shift left
 */
void rotl(void *x, size_t nx, size_t n);

#endif
