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

/** Get the sorting order of a list of integers.
 *
 * This function returns a permutation of x which would put x in sorted order.
 * It's like the order function in R, only zero-indexed.
 *
 * \param[in] x an array of integers to order
 * \param[out] order the result will be stored here
 * \param[in] n number of integers to order
 *
 * \sa https://stat.ethz.ch/R-manual/R-devel/library/base/html/order.html
 */
void order(const int *x, int *order, int n);

/** Circular left shift a block of memory.
 *
 * \param[in] x the block of memory to shift
 * \param[in] nx the size of x, in bytes
 * \param[in] n the number of bytes to shift left
 */
void rotl(void *x, size_t nx, size_t n);

#endif
