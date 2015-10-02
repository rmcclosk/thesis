#ifndef RWRAPPER_H
#define RWRAPPER_H

#include <Rinternals.h>

/** Initialize an embedded R session. */
void start_R(void);

/** Terminate an embedded R session. */
void stop_R(void);

/** Call an R function.
 *
 * \param func the name of the function to call
 * \param args arguments to supply to the function
 * \param tags names for the arguments, if any
 * \param narg number of arguments
 */
SEXP call_R(const char *func, const SEXP *args, const char **tags, int narg);

/** Import an R library into the embedded session. */
void R_library(const char *libname);

#endif /* RWRAPPER_H */
