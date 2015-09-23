#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/treestats.h"
#include "../src/util.h"

Suite *tree_suite(void);

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

START_TEST(test_kernel_no_coal)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    igraph_t *t2 = tree_from_newick("(((1:0.25,2:0.25)5:0.5,3:0.25)6:0.5,4:0.25)7;");
    ck_assert(kernel(t1, t2, 0.5, 1, 1, 0) == 1.125 * (1 + exp(-0.0625)));
    igraph_destroy(t1);
    igraph_destroy(t2);
}
END_TEST

START_TEST(test_kernel_coal_self)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    ck_assert(kernel(t1, t1, 0.5, 1, 1, 1) == 0.0);
    igraph_destroy(t1);
}
END_TEST

START_TEST(test_kernel_coal_L2)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    igraph_t *t2 = tree_from_newick("(((1:0.25,2:0.25)5:0.5,3:0.25)6:0.5,4:0.25)7;");
    ck_assert(fabs(kernel(t1, t2, 0.5, 1, 1, 2) - 1.125 / 144.0 * (1 + exp(-0.0625)) < 1e-5));
    igraph_destroy(t1);
    igraph_destroy(t2);
}
END_TEST

Suite *tree_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("tree");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_kernel_no_coal);
    tcase_add_test(tc_core, test_kernel_coal_self);
    tcase_add_test(tc_core, test_kernel_coal_L2);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    s = tree_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
