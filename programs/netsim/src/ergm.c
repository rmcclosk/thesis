#define CSTACK_DEFNS

#include <stdint.h>
#include <R.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <Rinterface.h>
#include <igraph/igraph.h>
#include "ergm.h"

struct ergm {
    SEXP model;
};

void R_library(const char *libname);

void start_R(void)
{
    char *argv[] = {"R", "--vanilla", "--silent"};
    int argc = sizeof(argv)/sizeof(argv[0]);
    SEXP call, opt, devnull;
    Rf_initEmbeddedR(argc, argv);
    R_SignalHandlers = 0;

    // turn off warnings
    PROTECT(opt = ScalarInteger(-1));
    PROTECT(call = lang2(install("options"), opt));
    SET_TAG(CDR(call), install("warn"));
    eval(call, R_GlobalEnv);
    UNPROTECT(2);

    // turn off all output (TODO: unportable)
    PROTECT(devnull = mkString("/dev/null"));
    PROTECT(call = lang2(install("sink"), devnull));
    eval(call, R_GlobalEnv);
    unprotect(2);

    // disable stack checking because R sucks at multithreading
    R_CStackLimit = (uintptr_t) - 1;

    R_library("network");
    R_library("ergm");
}

void stop_R(void)
{
    Rf_endEmbeddedR(0);
}

ergm *ergm_create(int nnode)
{
    SEXP call, n, directed;
    SEXP net;
    SEXP formula, model;
    ergm *m = malloc(sizeof(ergm));

    // create an empty network
    PROTECT(n = ScalarInteger(nnode));
    PROTECT(directed = ScalarLogical(0));
    PROTECT(call = lang3(install("network.initialize"), n, directed));
    PROTECT(net = eval(call, R_GlobalEnv));
    defineVar(install("net"), net, R_GlobalEnv);

    // create the model structure
    PROTECT(formula = mkString("net ~ edges + triangles"));
    PROTECT(call = lang2(install("as.formula"), formula));
    PROTECT(formula = eval(call, R_GlobalEnv));
    PROTECT(call = lang2(install("ergm"), formula));
    PROTECT(model = eval(call, R_GlobalEnv));

    m->model = model;

    // will model be GC'ed? I don't know how to handle this.
    UNPROTECT(9);
    return m;
}

void simulate_from_ergm(igraph_t *g, ergm *e, int nnode, int ncoef, double *coef)
{
    int i;
    SEXP call, net, el, nedge;
    SEXP f = PROTECT(ScalarLogical(0));
    igraph_vector_t edge;

    // copy the coefficients into an R vector
    SEXP vcoef = PROTECT(allocVector(REALSXP, ncoef));
    for (i = 0; i < ncoef; i++) {
        REAL(vcoef)[i] = coef[i];
    }

    // simulate a network from the model
    PROTECT(call = lang3(install("simulate.ergm"), e->model, vcoef));
    SET_TAG(CDDR(call), install("coef"));
    PROTECT(net = eval(call, R_GlobalEnv));

    // get the edgelist of the simulated network
    PROTECT(call = lang2(install("as.edgelist"), net));
    PROTECT(el = eval(call, R_GlobalEnv));

    // how many edges?
    PROTECT(call = lang2(install("nrow"), el));
    PROTECT(nedge = eval(call, R_GlobalEnv));

    // put the edges into an igraph vector (edge list is in column major order)
    igraph_vector_init(&edge, 0);
    for (i = 0; i < asInteger(nedge); ++i)
    {
        igraph_vector_push_back(&edge, INTEGER(el)[i]-1);
        igraph_vector_push_back(&edge, INTEGER(el)[i+asInteger(nedge)]-1);
    }

    UNPROTECT(8);

    igraph_empty(g, nnode, 0);
    igraph_add_edges(g, &edge, 0);
    igraph_vector_destroy(&edge);
}

void ergm_destroy(ergm *model)
{
    free(model);
}

void R_library(const char *libname)
{
    SEXP arg = PROTECT(mkString(libname));
    SEXP call = PROTECT(lang2(install("library"), arg));
    SEXP call2 = PROTECT(lang2(install("suppressPackageStartupMessages"), call));
    eval(call2, R_GlobalEnv);
    UNPROTECT(3);
}

