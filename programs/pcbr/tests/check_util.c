#include <check.h>
#include "../src/util.h"

Suite *util_suite(void);

START_TEST (test_swap)
{
    double test[4] = {1, 2, 3, 4};
    swap(test, 0, 2);
    ck_assert(test[0] == 3);
    ck_assert(test[2] == 1);
}
END_TEST

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
    double x[2] = { 1e5, 1e-100 };
    int scale = get_scale(x, 2);
    ck_assert_int_eq(scale, 5);
    x[0] = 1e-100; x[1] = 1e5;
    scale = get_scale(x, 2);
    ck_assert_int_eq(scale, 5);
}
END_TEST

START_TEST (test_order)
{
    double x[5] = {-1, 8, -4, 7, 1};
    int result[5];
    order(x, x, result, 5);
    ck_assert_int_eq(result[0], 1);
    ck_assert_int_eq(result[1], 4);
    ck_assert_int_eq(result[2], 0);
    ck_assert_int_eq(result[3], 3);
    ck_assert_int_eq(result[4], 2);
}
END_TEST

START_TEST (test_inverse_perm)
{
    int p[5] = {1, 4, 0, 3, 2};
    int pinv[5];
    inverse_perm(p, pinv, 5);
    ck_assert_int_eq(pinv[0], 2);
    ck_assert_int_eq(pinv[1], 0);
    ck_assert_int_eq(pinv[2], 4);
    ck_assert_int_eq(pinv[3], 3);
    ck_assert_int_eq(pinv[4], 1);
}
END_TEST

Suite *util_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("util");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_swap);
    tcase_add_test(tc_core, test_get_scale_zero);
    tcase_add_test(tc_core, test_get_scale_nonzero);
    tcase_add_test(tc_core, test_get_scale_tiny);
    tcase_add_test(tc_core, test_order);
    tcase_add_test(tc_core, test_inverse_perm);
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
