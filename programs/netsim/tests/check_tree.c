#include <check.h>
#include <limits.h>
#include <igraph/igraph.h>
#include <gsl/gsl_rng.h>

#include "../src/tree.h"
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

START_TEST (test_parse_newick_topology)
{
    FILE *f = newick_file("(t3,((t2,t1),t4));");
    igraph_vector_t vec = new_vector();
    int res;

    igraph_t *tree = parse_newick(f);

    ck_assert_int_eq(igraph_vcount(tree), 7);
    ck_assert_int_eq(igraph_ecount(tree), 6);

    igraph_is_connected(tree, &res, IGRAPH_WEAK); ck_assert_int_eq(res, 1);
    igraph_is_connected(tree, &res, IGRAPH_STRONG); ck_assert_int_eq(res, 0);
    igraph_is_dag(tree, &res); ck_assert_int_eq(res, 1);

    igraph_degree(tree, &vec, igraph_vss_all(), IGRAPH_OUT, 1);
    ck_assert(igraph_vector_sum(&vec) == 6.0);
    ck_assert(igraph_vector_min(&vec) == 0.0);
    ck_assert(igraph_vector_max(&vec) == 2.0);

    igraph_degree(tree, &vec, igraph_vss_all(), IGRAPH_IN, 1);
    ck_assert(igraph_vector_sum(&vec) == 6.0);
    ck_assert(igraph_vector_min(&vec) == 0.0);
    ck_assert(igraph_vector_max(&vec) == 1.0);

    igraph_vector_destroy(&vec);
    igraph_destroy(tree);
    fclose(f);
}
END_TEST

START_TEST (test_parse_newick_branch_lengths)
{
    int i;
    FILE *f = newick_file("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    igraph_t *tree = parse_newick(f);
    igraph_vector_t bl = new_vector();
    igraph_vector_t size = new_vector();
    igraph_vector_t edge = new_vector();

    EANV(tree, "length", &bl);
    ck_assert(igraph_vector_sum(&bl) == 2.84);

    igraph_neighborhood_size(tree, &size, igraph_vss_all(), INT_MAX, IGRAPH_OUT, 0);
    for (i = 0; i < igraph_vcount(tree); ++i)
    {
        igraph_incident(tree, &edge, i, IGRAPH_IN);
        switch( (int) VECTOR(size)[i] )
        {
            case 5:
                ck_assert(VECTOR(bl)[(int) VECTOR(edge)[0]] == 0.46);
                break;
            case 3:
                ck_assert(VECTOR(bl)[(int) VECTOR(edge)[0]] == 0.88);
                break;
            default:
                break;
        }
    }

    fclose(f);
    igraph_destroy(tree);
}
END_TEST

START_TEST (test_write_newick)
{
    int res;
    FILE *f = newick_file("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    igraph_t *tree = parse_newick(f), *copy;

    fseek(f, 0, SEEK_SET);
    write_tree_newick(tree, f);

    fseek(f, 0, SEEK_SET);
    copy = parse_newick(f);

    igraph_isomorphic(tree, copy, &res);
    ck_assert_int_eq(res, 1);

    fclose(f);
    igraph_destroy(tree);
}
END_TEST

START_TEST (test_root)
{
    igraph_t *tree = tree_from_newick("(t1,((t3,t2),t4));");
    int rt = root(tree);
    igraph_vector_t v = new_vector();

    igraph_degree(tree, &v, igraph_vss_1(rt), IGRAPH_OUT, 1);
    ck_assert_int_eq((int) VECTOR(v)[0], 2);
    igraph_degree(tree, &v, igraph_vss_1(rt), IGRAPH_IN, 1);
    ck_assert_int_eq((int) VECTOR(v)[0], 0);

    igraph_neighborhood_size(tree, &v, igraph_vss_1(rt), INT_MAX, IGRAPH_OUT, 0);
    ck_assert_int_eq((int) VECTOR(v)[0], igraph_vcount(tree));

    igraph_vector_destroy(&v);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_height)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    ck_assert(height(tree) == 2.14);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_depths)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    double calc_depths[7]; 
    int i;

    depths(tree, calc_depths);
    qsort(calc_depths, 7, sizeof(double), compare_doubles);
    for (i = 0; i < 7; ++i)
        fprintf(stderr, "%f\n", calc_depths[i]);

    ck_assert(fabs(calc_depths[0] - 0.0) < 1e-5);
    ck_assert(fabs(calc_depths[1] - 0.4) < 1e-5);
    ck_assert(fabs(calc_depths[2] - 0.46) < 1e-5);
    ck_assert(fabs(calc_depths[3] - 0.71) < 1e-5);
    ck_assert(fabs(calc_depths[4] - 1.34) < 1e-5);
    ck_assert(fabs(calc_depths[5] - 1.39) < 1e-5);
    ck_assert(fabs(calc_depths[6] - 2.14) < 1e-5);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_ladderize)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    igraph_vs_t vs;
    igraph_vector_t v = new_vector();
    int i;

    ladderize(tree);
    for (i = 0; i < 7; ++i)
    {
        igraph_degree(tree, &v, igraph_vss_1(i), IGRAPH_OUT, 1);
        if (VECTOR(v)[0]) {
            igraph_vs_adj(&vs, i, IGRAPH_OUT);
            igraph_neighborhood_size(tree, &v, vs, INT_MAX, IGRAPH_OUT, 0);
            if (VECTOR(v)[0] != VECTOR(v)[1]) {
                ck_assert(VECTOR(v)[0] < VECTOR(v)[1]);
            }
            else {
                igraph_incident(tree, &v, i, IGRAPH_OUT);
                ck_assert(EAN(tree, "length", (int) VECTOR(v)[0]) 
                          <= EAN(tree, "length", (int) VECTOR(v)[1])); 
            }
            memset(&vs, 0, sizeof(igraph_vs_t));
        }
    }

    igraph_vs_destroy(&vs);
    igraph_vector_destroy(&v);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_scale_branches_none)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    igraph_vector_t v1 = new_vector();
    igraph_vector_t v2 = new_vector();

    EANV(tree, "length", &v1);
    scale_branches(tree, NONE);
    EANV(tree, "length", &v2);

    ck_assert(igraph_vector_maxdifference(&v1, &v2) == 0);

    igraph_vector_destroy(&v1);
    igraph_vector_destroy(&v2);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_scale_branches_mean)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.8,t4:0.25):0.4);");
    igraph_vector_t v1 = new_vector();
    igraph_vector_t v2 = new_vector();

    EANV(tree, "length", &v1);
    igraph_vector_scale(&v1, 1./0.45);
    scale_branches(tree, MEAN);
    EANV(tree, "length", &v2);

    ck_assert(fabs(igraph_vector_maxdifference(&v1, &v2) < 1e-5));

    igraph_vector_destroy(&v1);
    igraph_vector_destroy(&v2);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_scale_branches_median)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.8,t4:0.25):0.4);");
    igraph_vector_t v1 = new_vector();
    igraph_vector_t v2 = new_vector();

    EANV(tree, "length", &v1);
    igraph_vector_scale(&v1, 1./0.4);
    scale_branches(tree, MEDIAN);
    EANV(tree, "length", &v2);

    ck_assert(fabs(igraph_vector_maxdifference(&v1, &v2) < 1e-5));

    igraph_vector_destroy(&v1);
    igraph_vector_destroy(&v2);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_cut_at_time_extinct)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    igraph_vector_t v = new_vector();
    cut_at_time(tree, 0.6, 0);

    ck_assert_int_eq(igraph_vcount(tree), 5);
    ck_assert_int_eq(igraph_ecount(tree), 4);

    EANV(tree, "length", &v);
    igraph_vector_sort(&v);

    ck_assert(fabs(VECTOR(v)[0] - 0.14) < 1e-5);
    ck_assert(fabs(VECTOR(v)[1] - 0.14) < 1e-5);
    ck_assert(fabs(VECTOR(v)[2] - 0.4) < 1e-5);
    ck_assert(fabs(VECTOR(v)[3] - 0.46) < 1e-5);

    igraph_vector_destroy(&v);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_cut_at_time_extant)
{
    igraph_t *tree = tree_from_newick("(t1:0.4,((t3:0.05,t2:0.8):0.88,t4:0.25):0.46);");
    igraph_vector_t v = new_vector();
    cut_at_time(tree, 0.6, 1);

    ck_assert_int_eq(igraph_vcount(tree), 3);
    ck_assert_int_eq(igraph_ecount(tree), 2);

    EANV(tree, "length", &v);
    igraph_vector_sort(&v);

    ck_assert(fabs(VECTOR(v)[0] - 0.14) < 1e-5);
    ck_assert(fabs(VECTOR(v)[1] - 0.14) < 1e-5);

    igraph_vector_destroy(&v);
    igraph_destroy(tree);
}
END_TEST

START_TEST(test_subsample_tips)
{
    igraph_t *tree = tree_from_newick("((1:0.5,2:0.25)5:0.5,(3:0.25,4:0.25)6:0.5)7;");
    gsl_rng *rng = set_seed(2);

    subsample_tips(tree, 2, rng);
    ck_assert_int_eq(igraph_vcount(tree), 3);
    ck_assert_int_eq(igraph_ecount(tree), 2);
    igraph_destroy(tree);
}
END_TEST

Suite *tree_suite(void)
{
    Suite *s;
    TCase *tc_io, *tc_tree;

    s = suite_create("tree");

    tc_io = tcase_create("IO");
    tcase_add_test(tc_io, test_parse_newick_topology);
    tcase_add_test(tc_io, test_parse_newick_branch_lengths);
    tcase_add_test(tc_io, test_write_newick);
    suite_add_tcase(s, tc_io);

    tc_tree = tcase_create("Tree");
    tcase_add_test(tc_tree, test_root);
    tcase_add_test(tc_tree, test_height);
    tcase_add_test(tc_tree, test_depths);
    tcase_add_test(tc_tree, test_ladderize);
    tcase_add_test(tc_tree, test_scale_branches_none);
    tcase_add_test(tc_tree, test_scale_branches_mean);
    tcase_add_test(tc_tree, test_scale_branches_median);
    tcase_add_test(tc_tree, test_cut_at_time_extinct);
    tcase_add_test(tc_tree, test_cut_at_time_extant);
    tcase_add_test(tc_tree, test_subsample_tips);
    suite_add_tcase(s, tc_tree);

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