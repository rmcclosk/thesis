#include <stdio.h>
#include <string.h>
#include <check.h>
#include <math.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_odeiv2.h>

#include "../src/tree.h"
#include "../src/mmpp.h"

Suite *mmpp_suite(void); 
// TODO: this is copypasta from check_tree
FILE *newick_file(const char *newick)
{
    FILE *f = tmpfile();
    fprintf(f, "%s", newick);
    fseek(f, 0, SEEK_SET);
    return f;
}

igraph_vector_t new_vector(void)
{
    igraph_vector_t v;
    igraph_vector_init(&v, 0);
    return v;
}

igraph_t *tree_from_newick(const char *newick)
{
    FILE *f = newick_file(newick);
    igraph_t *tree = parse_newick(f);
    fclose(f);
    return tree;
}

START_TEST(test_guess_parameters)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    double theta[4];
    guess_parameters(tree, 2, theta);
    ck_assert(fabs(theta[0] - 1.0/0.88) < 1e-5);
    ck_assert(fabs(theta[1] - 1.0/0.46) < 1e-5);
    ck_assert(fabs(theta[2] - 0.007463) < 1e-5);
    ck_assert(fabs(theta[3] - 0.007463) < 1e-5);
}
END_TEST

START_TEST(test_likelihood_toy)
{
    igraph_t *tree = tree_from_newick("(t1:1,t2:1);");
    double theta[4] = {1, 2, 1, 3};
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);
    ck_assert(fabs(pow(10, likelihood(tree, 2, theta, w, 0)) - 0.0883478) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_likelihood_toy3)
{
    igraph_t *tree = tree_from_newick("((t1:0.5,t2:1):0.75,t3:1);");
    double theta[4] = {1, 2, 1, 3};
    mmpp_workspace *w = mmpp_workspace_create(tree, 2);
    ck_assert(fabs(pow(10, likelihood(tree, 2, theta, w, 0)) - 0.02248663) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

Suite *mmpp_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("mmpp");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_guess_parameters);
    tcase_add_test(tc_core, test_likelihood_toy);
    tcase_add_test(tc_core, test_likelihood_toy3);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    s = mmpp_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
