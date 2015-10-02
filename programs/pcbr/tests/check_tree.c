#include <check.h>
#include <stdlib.h>
#include <pll/pll.h>
#include <Rinternals.h>
#include "../src/tree.h"
#include "../src/likelihood.h"
#include "config.h"

#define NO_CHILD -1

Suite *tree_suite(void);

const char *tree_3 = "(v1:0.5,(v3:0.5,v4:0.25)v2:0.25)v0;";

START_TEST (test_flatten)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    double *bl = malloc(5*sizeof(double));
    int *order = malloc(10*sizeof(double));

    flatten(tree->tree, tree->nodes, order, bl, 0);
    ck_assert(order[0] == NO_CHILD);
    ck_assert(order[1] == NO_CHILD);
    ck_assert(order[2] == NO_CHILD);
    ck_assert(order[3] == NO_CHILD);
    ck_assert(order[4] == NO_CHILD);
    ck_assert(order[5] == NO_CHILD);
    ck_assert(order[6] == 1);
    ck_assert(order[7] == 2);
    ck_assert(order[8] == 0);
    ck_assert(order[9] == 3);
    ck_assert(bl[0] == 0.5);
    ck_assert(bl[1] == 0.5);
    ck_assert(bl[2] == 0.25);
    ck_assert(bl[3] == 0.25);
    ck_assert(bl[4] == 0);
        
    free(order);
    free(bl);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_output)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *p = pcbr_create(tree, 2);
    char buf[1000];
    double rates[2] = {1, 2};
    int i;

    for (i = 0; i < 5; ++i) {
        p->A[i] = i % 2;
        p->L[2*i] = p->L[2*i+1] = 1;
        p->scale[i] = 0;
    }

    get_clusters(p, rates);
    output_tree(p, rates, p->nnode, buf);
    fprintf(stderr, "%s\n", buf);
    ck_assert(strcmp(buf, "(1[&rate=2.000000,cluster=0,log10lik=0.000000]:0.500000,(2[&rate=4.000000,cluster=1,log10lik=0.000000]:0.500000,3[&rate=2.000000,cluster=0,log10lik=0.000000]:0.250000)4[&rate=4.000000,cluster=1,log10lik=0.000000]:0.250000)5[&rate=2.000000,cluster=0,log10lik=0.000000]:0.000000;") == 0);

    pcbr_free(p);
    pllNewickParseDestroy(&tree);
}
END_TEST

Suite *tree_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("tree");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_flatten);
    tcase_add_test(tc_core, test_output);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    setenv("R_HOME", R_HOME, 0);

    s = tree_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
