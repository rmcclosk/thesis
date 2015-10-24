#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/tree.h"
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

START_TEST(test_kernel)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    igraph_t *t2 = tree_from_newick("(((1:0.25,2:0.25)5:0.5,3:0.25)6:0.5,4:0.25)7;");

    ck_assert(fabs(kernel(t1, t2, 0.5, 1, 1) - 1.125 * (1 + exp(-0.0625))) < 1e-5);
    igraph_destroy(t1);
    igraph_destroy(t2);
}
END_TEST

START_TEST(test_nLTT)
{
    igraph_t *t1 = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    igraph_t *t2 = tree_from_newick("((((1:0.4,2:0.4)6:0.3,3:0.4)7:0.2,4:0.4)8:0.1,5:0.4)9;");
    ck_assert(fabs(nLTT(t1, t2) - 0.03095142) < 1e-5);
    igraph_destroy(t1);
    igraph_destroy(t2);
}
END_TEST

START_TEST(test_sackin)
{
    igraph_t *tree = tree_from_newick("(((t2:0.53,t3:0.33):0.09,(t1:0.9,t4:0.31):0.08):0.77,t5:0.32);");
    ck_assert(sackin(tree, 0, TREESHAPE_NORM_NONE) == 13);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_sackin_yule)
{
    igraph_t *tree = tree_from_newick("(((t2:0.53,t3:0.33):0.09,(t1:0.9,t4:0.31):0.08):0.77,t5:0.32);");
    ck_assert(fabs(sackin(tree, 0, TREESHAPE_NORM_YULE) - 0.0333333) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_sackin_pda)
{
    igraph_t *tree = tree_from_newick("(((t2:0.53,t3:0.33):0.09,(t1:0.9,t4:0.31):0.08):0.77,t5:0.32);");
    ck_assert(fabs(sackin(tree, 0, TREESHAPE_NORM_PDA) - 1.162755) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_sackin_branch_lengths)
{
    igraph_t *tree = tree_from_newick("(((t2:0.53,t3:0.33):0.09,(t1:0.9,t4:0.31):0.08):0.77,t5:0.32);");
    ck_assert(sackin(tree, 1, TREESHAPE_NORM_NONE) == 5.81);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_colless)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    ck_assert(colless(t, TREESHAPE_NORM_NONE) == 2);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_colless_yule)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    ck_assert(fabs(colless(t, TREESHAPE_NORM_YULE) - -0.09350639) < 1e-5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_colless_pda)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    ck_assert(fabs(colless(t, TREESHAPE_NORM_PDA) - 0.1788854) < 1e-5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_cophenetic)
{
    igraph_t *t = tree_from_newick("((t2:0.1,t5:0.39):0.93,((t4:0.75,t1:0.48):0.52,t3:0.15):0.53);");
    ck_assert(cophenetic(t, TREESHAPE_NORM_NONE) == 5);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_ladder_length)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert(ladder_length(t) == 5.0 / 8.0);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_il_nodes)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert(il_nodes(t) == 2.0 / 7.0);
    igraph_destroy(t);
}
END_TEST

START_TEST(test_bmi)
{
    igraph_t *t = tree_from_newick("(((t5:1,t1:1):1,(t7:1,t8:1):1):1,(t6:1,(t2:1,(t3:1,t4:1):1):1):1);");
    ck_assert(bmi(t) == 6.0 / 5.0);
    igraph_destroy(t);
}
END_TEST

Suite *tree_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("tree");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_kernel);
    tcase_add_test(tc_core, test_nLTT);
    tcase_add_test(tc_core, test_sackin);
    tcase_add_test(tc_core, test_sackin_yule);
    tcase_add_test(tc_core, test_sackin_pda);
    tcase_add_test(tc_core, test_colless);
    tcase_add_test(tc_core, test_colless_yule);
    tcase_add_test(tc_core, test_colless_pda);
    tcase_add_test(tc_core, test_cophenetic);
    tcase_add_test(tc_core, test_ladder_length);
    tcase_add_test(tc_core, test_il_nodes);
    tcase_add_test(tc_core, test_bmi);
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
