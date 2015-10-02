#include <stdlib.h>
#include <check.h>
#include <Rinternals.h>
#include "../src/rwrapper.h"
#include "config.h"

Suite *rwrapper_suite(void);

START_TEST (test_basic)
{
    int Rerror;
    SEXP expr, result, one;

    one = PROTECT(ScalarInteger(1));
    expr = PROTECT(lang3(install("+"), one, one));
    result = R_tryEval(expr, R_GlobalEnv, &Rerror);

    ck_assert(INTEGER(result)[0] == 2);

    UNPROTECT(2);
}
END_TEST

START_TEST (test_onearg_notag)
{
    SEXP ans;
    const char *tags[1] = { NULL };
    const SEXP args[] = { PROTECT(ScalarInteger(1)) };

    ans = call_R("exp", args, tags, 1);
    ck_assert(fabs(REAL(ans)[0] - 2.718282) < 1e-6);
    UNPROTECT(1);
}
END_TEST

START_TEST (test_twoarg_tag)
{
    SEXP args[2], ans;
    const char *tags[2] = { NULL, "each"};

    args[0] = PROTECT(allocVector(INTSXP, 2));
    INTEGER(args[0])[0] = 1;
    INTEGER(args[0])[1] = 2;
    args[1] = PROTECT(ScalarInteger(2));

    ans = call_R("rep", (SEXP *) &args, tags, 2);
    ck_assert(INTEGER(ans)[0] == 1);
    ck_assert(INTEGER(ans)[1] == 1);
    ck_assert(INTEGER(ans)[2] == 2);
    ck_assert(INTEGER(ans)[3] == 2);
    UNPROTECT(2);
}
END_TEST

Suite *rwrapper_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("rwrapper");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_basic);
    tcase_add_test(tc_core, test_onearg_notag);
    tcase_add_test(tc_core, test_twoarg_tag);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    setenv("R_HOME", R_HOME, 0);

    s = rwrapper_suite();
    sr = srunner_create(s);

    start_R();
    srunner_run_all(sr, CK_NORMAL);
    stop_R();

    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
