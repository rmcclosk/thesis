#include <float.h>

#define LOG_ZERO DBL_MIN_10_EXP / 2

/** Find the index of the maximum element in a double array.
 *
 * \param x an array of doubles
 * \param n the length of x
 * \return the index of the largest element in x
 */
int which_max(double *x, int n);

/** Find the maximum element in a double array.
 *
 * \param x an array of doubles
 * \param n the length of x
 * \return the largest element in x
 */
double double_max(double *x, int n);
double double_min(double *x, int n);

/** Sum the elements of a double array.
 *
 * \param x an array of doubles
 * \param n the number of elements to sum, up to the length of x
 * \return the sum of the first n elements of x
 */
double double_sum(double *x, int n);

/** Get a log10 scale factor for some numbers.
 *
 * The numbers can be divided by 10 to the power of the returned value, to keep
 * them in a reasonable range.
 *
 * \param x the numbers to scale
 * \param n the length of x
 * \return the base 10 logarithm of a scaling factor to divide both numbers by
 */
int get_scale(double *x, int n);

/** Compare two doubles, for use in qsort.
 *
 * \param a a double to compare
 * \param b a second double to compare
 * \return -1 if a < b, otherwise 1
 */
int compare_doubles(const void *a, const void *b);

/** Compute the exponential of an array of doubles.
 *
 * \param x doubles to compute exp of
 * \param result place to write each exp(x[i])
 * \param n length of x
 */
void double_exp(double *x, double *result, int n);

/** Compute the natural logarithm of an array of doubles.
 *
 * \param x doubles to compute log of
 * \param result place to write each log(x[i])
 * \param n length of x
 */
void double_log(double *x, double *result, int n);

/** Sort an array of doubles, and keep track of the indexes.
 *
 * The sorted version of x will be stored in sorted. The sort order (like you
 * would get from R's order function) is stored in order.
 *
 * \param x array of doubles to sort
 * \param sorted place to store sorted x
 * \param order place to store sort order
 * \param n length of x
 */
void order(const double *x, double *sorted, int *order, int n);

/** Compute the inverse of a permutation.
 *
 * \param p a permutation of 1, ..., n
 * \param pinv plate to store inverse permutation
 * \param n length of p
 */
void inverse_perm(const int *p, int *pinv, int n);
