#include <string.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include "rwrapper.h"

SEXP construct_call(const SEXP *args, const char **tags, int narg);

void start_R(void)
{
	char *argv[] = { "R", "--no-save", "--silent" };
    Rf_initEmbeddedR(3, argv);
}

void stop_R(void)
{
    Rf_endEmbeddedR(0);
}

SEXP call_R(const char *func, const SEXP *args, const char **tags, int narg)
{
    int Rerror;
    SEXP expr, result;
    expr = PROTECT(LCONS(install(func), construct_call(args, tags, narg)));
    result = R_tryEval(expr, R_GlobalEnv, &Rerror);
    UNPROTECT(1);
    return Rerror ? R_NilValue : result;
}

void R_library(const char *libname)
{
    SEXP expr;
    int Rerror;
    expr = PROTECT(lang2(install("library"), mkString(libname)));
    R_tryEval(expr, R_GlobalEnv, &Rerror);
    UNPROTECT(1);
}

// inspiration from Advanced R
SEXP construct_call(const SEXP *args, const char **tags, int narg)
{
    SEXP call;
    if (narg == 1) 
        call = PROTECT(LCONS(args[0], R_NilValue));
    else
        call = PROTECT(LCONS(args[0], construct_call(&args[1], &tags[1], narg-1)));
    if (tags[0])
        SET_TAG(call, install(tags[0]));
    UNPROTECT(1);
    return call;
}
