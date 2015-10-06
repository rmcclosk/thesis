#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/util.h"

Suite *util_suite(void);

START_TEST (test_get_scale_nonzero)
{
    double x[2] = { 1e5, 1e-10 };
    int scale = get_scale(x, 2);
    ck_assert_int_eq(scale, -2);
}
END_TEST

START_TEST (test_get_scale_zero)
{
    double x[2] = { 1e5, 0 };
    int scale = get_scale((double *) &x, 2);
    ck_assert_int_eq(scale, 5);
    x[0] = 0; x[1] = 1e5;
    scale = get_scale((double *) &x, 2);
    ck_assert_int_eq(scale, 5);
}
END_TEST

START_TEST (test_get_scale_tiny)
{
    double x[2] = { 1e5, 1e-200 };
    int scale = get_scale(x, 2);
    ck_assert_int_eq(scale, 5);
    x[0] = 1e-200; x[1] = 1e5;
    scale = get_scale(x, 2);
    ck_assert_int_eq(scale, 5);
}
END_TEST

START_TEST (test_set_seed)
{
    gsl_rng *rng = set_seed(0);
    ck_assert(rand() == 1804289383);
    ck_assert(gsl_rng_get(rng) == 4293858116);
    ck_assert(igraph_rng_get_integer(igraph_rng_default(), 0, INT_MAX) == 2146929058);
}
END_TEST

START_TEST (test_order_ints)
{
    int x[5] = {4, 2, 1, 5, 3};
    int ord[5];
    order(x, ord, sizeof(int), 5, compare_ints);

    ck_assert_int_eq(ord[0], 2);
    ck_assert_int_eq(ord[1], 1);
    ck_assert_int_eq(ord[2], 4);
    ck_assert_int_eq(ord[3], 0);
    ck_assert_int_eq(ord[4], 3);
}
END_TEST

START_TEST (test_order_doubles)
{
    double x[5] = {4, 2, 1, 5, 3};
    int ord[5];
    order(x, ord, sizeof(double), 5, compare_doubles);

    ck_assert_int_eq(ord[0], 2);
    ck_assert_int_eq(ord[1], 1);
    ck_assert_int_eq(ord[2], 4);
    ck_assert_int_eq(ord[3], 0);
    ck_assert_int_eq(ord[4], 3);
}
END_TEST

START_TEST (test_rotl_int)
{
    int x[5] = {0, 1, 2, 3, 4};
    rotl(x, sizeof(x), 2 * sizeof(int));
    ck_assert_int_eq(x[0], 2);
    ck_assert_int_eq(x[1], 3);
    ck_assert_int_eq(x[2], 4);
    ck_assert_int_eq(x[3], 0);
    ck_assert_int_eq(x[4], 1);
}
END_TEST

START_TEST (test_rotl_double)
{
    double x[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
    rotl(x, sizeof(x), 2 * sizeof(double));
    ck_assert(x[0] == 2.0);
    ck_assert(x[1] == 3.0);
    ck_assert(x[2] == 4.0);
    ck_assert(x[3] == 0.0);
    ck_assert(x[4] == 1.0);
}
END_TEST

Suite *util_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("util");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_get_scale_zero);
    tcase_add_test(tc_core, test_get_scale_nonzero);
    tcase_add_test(tc_core, test_get_scale_tiny);
    tcase_add_test(tc_core, test_set_seed);
    tcase_add_test(tc_core, test_order_ints);
    tcase_add_test(tc_core, test_order_doubles);
    tcase_add_test(tc_core, test_rotl_int);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = util_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
