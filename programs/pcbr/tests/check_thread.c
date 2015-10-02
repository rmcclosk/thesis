#include <check.h>
#include "../src/thread.h"
#include "../src/likelihood.h"
#include "../C-Thread-Pool/thpool.h"

Suite *thread_suite(void);

const char *tree_6 = "(((t5:0.4119189046,t6:0.7903786616):0.8664339224,t4:0.03210153081):0.5706496239,((t2:0.805275284,t1:0.7685161901):0.4054109394,t3:0.4271961174):0.3987961309);";

START_TEST (test_thread_data_create)
{
    int nthread = 3;
    int nrate = 2;
    int popsize = 100;

    pllNewickTree *tree = pllNewickParseString(tree_6);
    thread_data *tdata = create_thread_data(tree, nrate, nthread, popsize);

    ck_assert_int_eq(tdata[0].neval, 34);
    ck_assert_int_eq(tdata[1].neval, 33);
    ck_assert_int_eq(tdata[2].neval, 33);

    ck_assert_int_eq(memcmp(tdata[0].w->branch_lengths, tdata[1].w->branch_lengths, 11*sizeof(double)), 0);
    ck_assert_int_eq(memcmp(tdata[1].w->branch_lengths, tdata[2].w->branch_lengths, 11*sizeof(double)), 0);

    destroy_thread_data(tdata, nthread);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_set_args)
{
    int nthread = 3;
    int nrate = 2;
    int popsize = 10;
    int narg = nrate * nrate;
    int i, j;
    double **args = malloc(popsize * sizeof(double*));

    pllNewickTree *tree = pllNewickParseString(tree_6);
    thread_data *tdata = create_thread_data(tree, nrate, nthread, popsize);

    srand(0);
    for (i = 0; i < popsize; ++i) {
        args[i] = malloc(narg * sizeof(double));
        for (j = 0; j < narg; ++j)
            args[i][j] = (double) rand() / (double) RAND_MAX;
    }
    set_args(tdata, args, nthread, narg);

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < narg; ++j) {
            ck_assert(tdata[0].args[i][j] == args[i][j]);
            ck_assert(tdata[1].args[i][j] == args[i+4][j]);
            ck_assert(tdata[2].args[i][j] == args[i+7][j]);
        }
    }

    destroy_thread_data(tdata, nthread);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_likelihood_serial)
{
    int nthread = 3;
    int nrate = 2;
    int popsize = 10;
    int narg = nrate * nrate;
    int i, j;
    double **args = malloc(popsize * sizeof(double*));
    double *ans = malloc(popsize * sizeof(double));

    pllNewickTree *tree = pllNewickParseString(tree_6);
    pcbr_workspace *w = pcbr_create(tree, nrate);
    thread_data *tdata = create_thread_data(tree, nrate, nthread, popsize);

    srand(0);
    for (i = 0; i < popsize; ++i) {
        args[i] = malloc(narg * sizeof(double));
        for (j = 0; j < narg; ++j)
            args[i][j] = (double) rand() / (double) RAND_MAX;
    }
    set_args(tdata, args, nthread, narg);

    for (i = 0; i < popsize; ++i)
        ans[i] = likelihood(w, args[i], 0);

    for (i = 0; i < nthread; ++i)
        do_likelihood(&tdata[i]);

    ck_assert(tdata[0].arFunvals[0] = ans[0]);
    ck_assert(tdata[0].arFunvals[0] = ans[1]);
    ck_assert(tdata[0].arFunvals[0] = ans[2]);
    ck_assert(tdata[0].arFunvals[0] = ans[3]);
    ck_assert(tdata[1].arFunvals[0] = ans[4]);
    ck_assert(tdata[1].arFunvals[0] = ans[5]);
    ck_assert(tdata[1].arFunvals[0] = ans[6]);
    ck_assert(tdata[2].arFunvals[0] = ans[7]);
    ck_assert(tdata[2].arFunvals[0] = ans[8]);
    ck_assert(tdata[2].arFunvals[0] = ans[9]);

    free(ans);
    destroy_thread_data(tdata, nthread);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_likelihood_parallel)
{
    int nthread = 3;
    int nrate = 2;
    int popsize = 10;
    int narg = nrate * nrate;
    int i, j, iter;
    double **args = malloc(popsize * sizeof(double*));
    double *ans = malloc(popsize * sizeof(double));
    threadpool thpool = thpool_init(nthread);

    pllNewickTree *tree = pllNewickParseString(tree_6);
    pcbr_workspace *w = pcbr_create(tree, nrate);
    thread_data *tdata = create_thread_data(tree, nrate, nthread, popsize);

    srand(0);
    for (i = 0; i < popsize; ++i) {
        args[i] = malloc(narg * sizeof(double));
        for (j = 0; j < narg; ++j)
            args[i][j] = (double) rand() / (double) RAND_MAX;
    }

    for (iter = 0; iter < 2; ++iter) {
        fprintf(stderr, "starting iteration %d\n", iter);

        for (i = 0; i < popsize; ++i) {
            for (j = 0; j < narg; ++j)
                args[i][j] = (double) rand() / (double) RAND_MAX;
        }
        set_args(tdata, args, nthread, narg);
    
        for (i = 0; i < popsize; ++i)
            ans[i] = likelihood(w, args[i], 0);
    
        for (i = 0; i < nthread; ++i)
            thpool_add_work(thpool, do_likelihood, &tdata[i]);
        thpool_wait(thpool);
    
        ck_assert(tdata[0].arFunvals[0] = ans[0]);
        ck_assert(tdata[0].arFunvals[0] = ans[1]);
        ck_assert(tdata[0].arFunvals[0] = ans[2]);
        ck_assert(tdata[0].arFunvals[0] = ans[3]);
        ck_assert(tdata[1].arFunvals[0] = ans[4]);
        ck_assert(tdata[1].arFunvals[0] = ans[5]);
        ck_assert(tdata[1].arFunvals[0] = ans[6]);
        ck_assert(tdata[2].arFunvals[0] = ans[7]);
        ck_assert(tdata[2].arFunvals[0] = ans[8]);
        ck_assert(tdata[2].arFunvals[0] = ans[9]);
    }

    thpool_destroy(thpool);
    free(ans);
    destroy_thread_data(tdata, nthread);
    pllNewickParseDestroy(&tree);
}
END_TEST

Suite *thread_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("thread");

    tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_thread_data_create);
    tcase_add_test(tc_core, test_set_args);
    tcase_add_test(tc_core, test_likelihood_serial);
    tcase_add_test(tc_core, test_likelihood_parallel);
    tcase_set_timeout(tc_core, 10);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = thread_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
