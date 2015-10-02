#include <check.h>
#include <math.h>
#include <float.h>
#include <pll/pll.h>
#include <Rinternals.h>
#include <gsl/gsl_linalg.h>
#include "../src/likelihood.h"
#include "../src/rwrapper.h"
#include "config.h"

Suite *likelihood_suite(void);

const double x[100] = {0.0728086347350113, 0.602110493947618,
    0.0539702775039816, 0.123037593874584, 0.173192410891559,
    0.150060501350615, 0.0287769297555776, 0.255294705780355,
    0.081917778758616, 0.0435785336980695, 0.140959135802653,
    0.0107712528037133, 0.024278037816916, 0.0924111232999797,
    0.00219484240094387, 0.0121767540461289, 0.114765753806295,
    0.01170882371098, 0.0214211474907487, 0.259291104498715, 0.26757185192475,
    0.0157678332618751, 0.048926274984119, 0.0434253490942797,
    0.0804473073095657, 0.019793183477126, 0.0640305909809228,
    0.369301584555492, 0.160194032953655, 0.00608590489983153,
    0.0799107236015473, 0.564968985513339, 0.089209360045027,
    0.0122622598338913, 0.00697178835952705, 0.570939519847807,
    0.0132819668817668, 0.273229690249203, 0.228839048124922,
    0.0670438044965008, 0.166398110916178, 0.00289262732177717,
    0.00122990820321939, 0.415458348641762, 0.168542802986151,
    0.434607616874106, 1, 0.0421101252572851, 0.0925862260471381,
    0.0624358836292899, 0.211621250667168, 0.27843909739698, 0.665821166185745,
    0.478493943026058, 0.169366005273815, 0.215942261710302, 0.212763655045616,
    0.0111670546426562, 0.0274315213357667, 0.00304265855655133,
    0.267393279213092, 0.0254419551992687, 0.223447750658408,
    0.0873699162848675, 0.393113377436661, 0.325273692470244,
    0.0135479413703481, 0.323877921189978, 0.0388187018915403,
    0.00902308135654281, 0.202241604924605, 0.059094154785108,
    0.168301212751618, 0.011974651076497, 0.0218204847359255,
    0.0296987576209323, 0.208627254975382, 0.0185552561306158,
    0.229302626703307, 0.00771156019359826, 0.00301494666554995,
    0.020625937893032, 0.0256083489278167, 0.0721915712064602,
    0.00380806356435934, 0.0545724152683169, 0.0486316071567948,
    0.115425402211494, 0.106392348141476, 0.0915395330054313, 0.45082908924029,
    0.235231836051702, 0.155508340877168, 0.0623411276785619,
    0.205792328520891, 0.214181099483814, 0.0511636321110566,
    0.472572308008228, 0.327696553469504, 0.0599023527351818};

const char *tree_100 = "(((t101:0.6037924015,(t29:0.1621090819,(t85:0.001932834741,(t9:0.2609297372,t87:0.9384137455):0.8370064907):0.5946176518):0.260828245):2.909887269,(((t46:0.690256723,(t11:0.9551467383,t32:0.7125401224):0.3958932852):1.233791508,t79:0.3971479333):0.1390735129,(((t26:0.4359764017,(t55:0.6919276721,(t77:0.1554012226,(((t7:0.5022508153,(t63:0.001836858224,t62:0.8775780622):0.05658655246):0.5546413977,t41:0.1341113378):0.05884797215,t38:0.02274122415):0.01060726233):0.4466056545):0.1173312106):0.05205544783,((t48:0.1643265744,t78:0.3991025558):1.253105354,(t20:0.4340308486,t19:0.5170098264):1.293124656):0.1035243946):0.681229108,(((t35:0.6882752376,t27:0.6580575537):0.2098665806,((t69:0.9695281661,(t14:0.8495521038,t86:0.7566449074):0.3094478562):0.09565674937,t33:0.5326012194):0.3887867721):0.2364515253,((t43:0.008128455374,((t49:0.187426371,t72:0.6460672768):0.3861935635,(t17:0.3353207598,t18:0.6379087232):2.730389314):0.02941203876):0.7741877641,t56:0.8292010638):1.784765405):0.07620298555):0.2106068931):0.7252143034):0.3518704979,(((t58:0.3880784889,t8:0.9281775483):0.03369334764,(((t65:0.9939138787,(t10:0.02937716548,(t66:0.2776580867,t60:0.1171975501):0.3240101528):1.105936268):1.320467929,(t95:0.3703097857,t67:0.3368783088):0.8041709112):0.06418925882,(t22:0.6217732795,t40:0.3978436275):0.01397952619):2.759243788):0.05926120561,(((t82:0.1971467135,t21:0.1153423248):0.8145358063,((t99:0.5619880927,t37:0.7327180177):4.832812745,t96:0.8708055529):2.10037723):2.007832402,((((((t52:0.0496535839,t100:0.8211623181):3.217789018,t50:0.8293243044):1.345644019,(((t70:0.731371579,t59:0.9072914161):1.043608515,t39:0.6961969994):0.8185141888,t1:0.2415792223):2.312471626):1.022725877,(((((t45:0.2520097974,t97:0.09439026774):1.292261648,t73:0.8277177697):0.01470459905,(((t76:0.8425899034,t81:0.7373054679):0.4222424449,(t30:0.9489381774,(t83:0.03527776687,t98:0.5964484557):1.571986847):1.899843341):1.079881137,t13:0.4153180015):0.1229562053):0.1325714059,((t80:0.9623333085,t31:0.7087400511):1.565241345,t16:0.5534757036):0.06547466372):0.053968284,(t54:0.778042529,(t23:0.8302457039,t5:0.6485509474):0.04360686258):0.1876035172):1.028246904):0.3017409341,(((t12:0.4504854372,t74:0.8142517505):0.813368246,(t3:0.1474810452,(t28:0.9756573497,t24:0.9747924635):0.1054543167):0.05787124634):0.2855909844,((((t102:0.9347601165,t25:0.3461620999):1.108176657,t2:0.5330606059):0.08967407831,t88:0.5387942994):1.008256457,(t51:0.4057905001,(t42:0.3402327613,t47:0.6266548475):0.01457067267):0.03726852639):0.1435285343):0.9773958058):0.4474518932,((((t15:0.1364875862,(t53:0.6210762938,t57:0.2559822486):0.263738264):0.01840365813,t36:0.6348757988):0.3488883454,((t75:0.8575015371,(t44:0.31420183,t68:0.8285343598):0.5141742961):0.5578293549,((t84:0.09780853754,t4:0.06490053656):2.178772568,(t64:0.6680506039,t89:0.9045466478):1.136831415):0.4423934218):0.2350274509):0.1237603551,(((t71:0.792378756,(t94:0.0307565711,t34:0.8620340412):1.035097147):0.9945557881,((t6:0.6758537604,t61:0.8431201484):2.283853473,t93:0.361894185):0.2472642533):0.3012829964,(t90:0.5676873976,(t92:0.1937840716,t91:0.588066393):0.2894968537):1.58369608):0.7515426917):0.09968129552):0.20351035):0.00594391604):0.4311321322);";

const char *tree_3 = "(v1:0.5,(v3:0.5,v4:0.25)v2:0.25)v0;";

const char *tree_6 = "(((t5:0.4119189046,t6:0.7903786616):0.8664339224,t4:0.03210153081):0.5706496239,((t2:0.805275284,t1:0.7685161901):0.4054109394,t3:0.4271961174):0.3987961309);";

double _tree3_likelihood(double *lambda, double *q, int terminals, int *A, int n);
double _Pmut(int i, int j, double t, double *q, int n);
double likassign(pcbr_workspace *tree, int *v, double *lambda, double *q, double *pi, int terminals, int n);
double liklong(pcbr_workspace *tree, double *lambda, double *q, int terminals, int nrates);
void reclong(pcbr_workspace *tree, double *lambda, double *q, int *A, int terminals, int n);

START_TEST (test_guess_means)
{
    double guess[2];
    guess_means(x, 2, 100, guess);
    ck_assert(fabs(guess[0] - 0.03218813 < 1e-5));
    ck_assert(fabs(guess[1] - 0.27499070 < 1e-5));
}
END_TEST

START_TEST (test_guess_parameters)
{
    pllNewickTree *tree = pllNewickParseString(tree_100);
    pcbr_workspace *w = pcbr_create(tree, 2);
    double args[4];

    guess_parameters(w, args);
    ck_assert(fabs(args[0] - 5.097274) < 1e-5);
    ck_assert(fabs(args[1] - 37.580371) < 1e-5);
    ck_assert(fabs(args[2] - 0.032554) < 1e-5);
    ck_assert(fabs(args[3] - 0.032554) < 1e-5);

    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_scale_branches)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *p = pcbr_create(tree, 2);
    scale_branches(p);

    ck_assert(fabs(p->t[0] - 1) < 1e-6);
    ck_assert(fabs(p->t[1] - 1) < 1e-6);
    ck_assert(fabs(p->t[2] - 0.5) < 1e-6);
    ck_assert(fabs(p->t[3] - 0.5) < 1e-6);
    ck_assert(fabs(p->t[4]) < 1e-6);
}
END_TEST

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

START_TEST (test_pi)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *w = pcbr_create(tree, 3);
    double q[6] = { 0.560869, 0.436442, 0.008089, 0.532517, 0.266445, 0.568684 };
    fill_parameters(w, q);

    ck_assert(fabs(w->pi[0] - 0.106184) < 1e-5);
    ck_assert(fabs(w->pi[1] - 0.511907) < 1e-5);
    ck_assert(fabs(w->pi[2] - 0.381908) < 1e-5);

    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_pmut)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *w = pcbr_create(tree, 2);
    double p, t = 0.5, q[2] = {1, 2};
    int n = 2;

    ck_assert(w->t[n] == t);
    fill_parameters(w, q);
    kb(w);
    p = w->P[n*4];

    ck_assert(fabs(p - 0.7410434) < 1e-5);
    p = w->P[n*4+1];
    ck_assert(fabs((1-p) - 0.7410434) < 1e-5);

    p = _Pmut(0, 0, t, q, 2);
    ck_assert(fabs(p - 0.7410434) < 1e-5);
    p = _Pmut(0, 1, t, q, 2);
    ck_assert(fabs((1-p) - 0.7410434) < 1e-5);

    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_pmut_3)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *w = pcbr_create(tree, 3);
    double t = 0.5, q[6] = {1, 2, 3, 1, 2, 3};
    int n = 2;

    ck_assert(w->t[n] == t);
    fill_parameters(w, q);
    kb(w);

    ck_assert(fabs(w->P[n*9] - 0.4888829) < 1e-5);
    ck_assert(fabs(w->P[n*9+1] - 0.2655994) < 1e-5);
    ck_assert(fabs(w->P[n*9+2] - 0.2455177) < 1e-5);
    ck_assert(fabs(w->P[n*9+3] - 0.4451907) < 1e-5);
    ck_assert(fabs(w->P[n*9+4] - 0.3331609) < 1e-5);
    ck_assert(fabs(w->P[n*9+5] - 0.2216485) < 1e-5);
    ck_assert(fabs(w->P[n*9+6] - 0.4213215) < 1e-5);
    ck_assert(fabs(w->P[n*9+7] - 0.3133378) < 1e-5);
    ck_assert(fabs(w->P[n*9+8] - 0.2653407) < 1e-5);

    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_pt_nonterminal)
{
    double t = 0.5, lambda[2] = {1, 2};
    double p = Pt(0, t, lambda, 0);
    ck_assert(fabs(p - 0.6065307) < 1e-5);
    p = Pt(1, t, lambda, 0);
    ck_assert(fabs(p - 0.7357589) < 1e-5);
}
END_TEST

START_TEST (test_pt_nonterminal_3)
{
    double t = 0.5, lambda[3] = {1, 2, 3};
    double p = Pt(0, t, lambda, 0);
    ck_assert(fabs(p - 0.6065307) < 1e-5);
    p = Pt(1, t, lambda, 0);
    ck_assert(fabs(p - 0.7357589) < 1e-5);
}
END_TEST

START_TEST (test_pt_terminal)
{
    double t = 0.5, lambda[2] = {1, 2};
    double p = Pt(0, t, lambda, 1);
    ck_assert(fabs(p - 0.6065307) < 1e-5);
    p = Pt(1, t, lambda, 1);
    ck_assert(fabs(p - 0.3678794) < 1e-5);
}
END_TEST

START_TEST (test_pt_terminal_3)
{
    double t = 0.5, lambda[3] = {1, 2, 3};
    double p = Pt(0, t, lambda, 1);
    ck_assert(fabs(p - 0.6065307) < 1e-5);
    p = Pt(1, t, lambda, 1);
    ck_assert(fabs(p - 0.3678794) < 1e-5);
}
END_TEST

double _Pmut(int i, int j, double t, double *q, int n)
{
    double res;
    int k, l, cur = 0;
    double rowsum = 0;
    gsl_matrix *Q = gsl_matrix_alloc(n, n);
    gsl_matrix *P = gsl_matrix_alloc(n, n);

    for (k = 0; k < n; ++k) {
        rowsum = 0;
        for (l = 0; l < n; ++l) {
            if (k != l) {
                gsl_matrix_set(Q, k, l, q[cur] * t);
                rowsum += q[cur++];
            }
        }
        gsl_matrix_set(Q, k, k, -rowsum * t);
    }

    gsl_linalg_exponential_ss(Q, P, 0);
    res = gsl_matrix_get(P, i, j);
    gsl_matrix_free(P);
    return res;
}

double _tree3_likelihood(double *lambda, double *q, int terminals, int *A, int n)
{
    int v0, v1, v2, v3, v4;
    double t1 = 1, t2 = 0.5, t3 = 1, t4 = 0.5;
    double lik, likmax = 0, lik_total = 0;
    int i;
    double *pi = malloc(n*sizeof(double));

    for (i = 0; i < n; ++i)
        pi[i] = _Pmut(0, i, 100, q, n);

    for (v0 = 0; v0 < n; ++v0) {
    for (v1 = 0; v1 < n; ++v1) {
    for (v2 = 0; v2 < n; ++v2) {
    for (v3 = 0; v3 < n; ++v3) {
    for (v4 = 0; v4 < n; ++v4) {
        lik = pi[v0];
        lik *= _Pmut(v0, v1, t1, q, n);
        lik *= _Pmut(v0, v2, t2, q, n);
        lik *= _Pmut(v2, v3, t3, q, n);
        lik *= _Pmut(v2, v4, t4, q, n);
        lik *= Pt(v0, t2, lambda, 0);
        if (terminals) {
            lik *= Pt(v0, t1, lambda, 1);
            lik *= Pt(v2, t3, lambda, 1);
            lik *= Pt(v2, t4, lambda, 1);
        }
        lik_total += lik;
        if (likmax < lik) {
            likmax = lik;
            A[0] = v0; A[1] = v1; A[2] = v2; A[3] = v3; A[4] = v4;
        }
    }}}}}

    free(pi);
    return lik_total;
}

double likassign(pcbr_workspace *tree, int *v, double *lambda, double *q, double *pi, int terminals, int n)
{
    double lik = pi[v[tree->nnode-1]];
    int i, lc, rc, lterm, rterm;

    for (i = 0; i < tree->nnode; ++i) {
        lc = tree->topology[i*2];
        rc = tree->topology[i*2+1];
        if (lc != -1) {
            lterm = tree->topology[lc*2] == NO_CHILD;
            rterm = tree->topology[rc*2] == NO_CHILD;

            lik *= _Pmut(v[i], v[lc], tree->t[lc], q, n);
            lik *= _Pmut(v[i], v[rc], tree->t[rc], q, n);
            if (!lterm || terminals)
                lik *= Pt(v[i], tree->t[lc], lambda, lterm);
            if (!rterm || terminals)
                lik *= Pt(v[i], tree->t[rc], lambda, rterm);
        }
    }
    return lik;
}

double liklong(pcbr_workspace *tree, double *lambda, double *q, int terminals, int nrates)
{
    int i, j, icopy;
    int *v = malloc(tree->nnode*sizeof(int));
    double lik = 0;
    fill_parameters(tree, q);

    for (i = 0; i < (int) pow(nrates, tree->nnode); ++i) {
        icopy = i;
        j = 0;
        memset(v, 0, tree->nnode*sizeof(int));
        while (icopy > 0) {
            v[j] = icopy % nrates;
            icopy /= nrates;
            ++j;
        }
        lik += likassign(tree, v, lambda, q, tree->pi, terminals, nrates);
    }
    free(v);
    return lik;
}

void reclong(pcbr_workspace *tree, double *lambda, double *q, int *A, int terminals, int n)
{
    int i, j, icopy;
    int *v = malloc(tree->nnode*sizeof(int));
    double lik = 0, likmax = 0;
    fill_parameters(tree, q);

    for (i = 0; i < (int) pow(n, tree->nnode); ++i) {
        icopy = i;
        j = 0;
        memset(v, 0, tree->nnode*sizeof(int));
        while (icopy > 0) {
            v[j++] = icopy % n;
            icopy /= n;
        }
        lik = likassign(tree, v, lambda, q, tree->pi, terminals, n);
        if (likmax < lik) {
            likmax = lik;
            memcpy(A, v, tree->nnode*sizeof(int));
        }
    }
    free(v);
}

START_TEST(test_likelihood_terminals)
{
    double args[4] = {1, 2, 1, 2};
    pllNewickTree *tree = pllNewickParseString(tree_6);
    pcbr_workspace *w = pcbr_create(tree, 2);

    double lik = liklong(w, args, &args[2], 1, 2);
    double p = pow(10, likelihood(w, args, 0));

    ck_assert(fabs(p - lik) < 1e-5);
    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST(test_likelihood_3)
{
    double args[9] = {1, 2, 3, 1, 2, 3, 1, 2, 3};

    pllNewickTree *tree = pllNewickParseString(tree_6);
    pcbr_workspace *w = pcbr_create(tree, 3);

    double p = pow(10, likelihood(w, args, 0));
    double lik = liklong(w, args, &args[3], 1, 3);

    ck_assert(fabs(p - lik) < 1e-5);
    pcbr_free(w);
    pllNewickParseDestroy(&tree);

}
END_TEST

START_TEST(test_reconstruct_terminals)
{
    int i;
    double args[4] = {1, 2, 1, 2};
    pllNewickTree *tree = pllNewickParseString(tree_6);
    pcbr_workspace *w = pcbr_create(tree, 2);
    int *A = malloc(w->nnode*sizeof(int));

    likelihood(w, args, 1);
    reclong(w, args, &args[2], A, 1, 2);

    for (i = 0; i < w->nnode; ++i)
        ck_assert(w->A[w->nnode-1-i] == A[i]);

    free(A);
    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST(test_reconstruct_3)
{
    int i;
    double args[10] = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    pllNewickTree *tree = pllNewickParseString(tree_6);
    pcbr_workspace *w = pcbr_create(tree, 3);
    int *A = malloc(w->nnode*sizeof(int));

    likelihood(w, args, 1);
    reclong(w, args, &args[3], A, 1, 3);

    for (i = 0; i < w->nnode; ++i)
        ck_assert(w->A[w->nnode-1-i] == A[i]);

    free(A);
    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_kb)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *w = pcbr_create(tree, 2);
    double t, q[2] = {0.5, 0.8};
    int i, j, node;

    fill_parameters(w, q);
    kb(w);

    for (node = 0; node < 4; ++node) {
        t = w->t[node];

        for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j)
                ck_assert(fabs(w->P[(node*4)+i*2+j] - _Pmut(i, j, t, q, 2) < 1e-6));
        }
    }
    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_kb_large)
{
    pllNewickTree *tree = pllNewickParseString(tree_100);
    pcbr_workspace *w = pcbr_create(tree, 2);
    double t, q[2] = {0.5, 0.8};
    int i, j, node;

    fill_parameters(w, q);
    kb(w);

    for (node = 0; node < tree->nodes-1; ++node) {
        t = w->t[node];

        for (i = 0; i < 2; ++i) {
            for (j = 0; j < 2; ++j) {
                ck_assert(fabs(w->P[(node*4)+i*2+j] - _Pmut(i, j, t, q, 2) < 1e-6));
            }
        }
    }
    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_liklong)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *w = pcbr_create(tree, 2);
    double lambda[2] = {0.5, 0.8};
    double q[2] = {0.5, 0.8};
    int A[5];
    double lik1 = liklong(w, lambda, q, 1, 2);
    double lik2 = _tree3_likelihood(lambda, q, 1, A, 2);

    ck_assert(fabs(lik1 - lik2) < 1e-6);

    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_liklong_3)
{
    pllNewickTree *tree = pllNewickParseString(tree_3);
    pcbr_workspace *w = pcbr_create(tree, 3);
    double lambda[3] = {1, 2, 3};
    double q[6] = {1, 2, 3, 1, 2, 3};
    double A[5];
    double lik1 = liklong(w, lambda, q, 1, 3);
    double lik2 = _tree3_likelihood(lambda, q, 1, (int *) &A, 3);
    ck_assert(fabs(lik1 - lik2) < 1e-6);

    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_likelihood_random)
{
    int ntip = 5;
    SEXP rargs[] = { PROTECT(ScalarInteger(ntip)) }, rtree, treestr;
    const char *tags[] = { NULL };
    pllNewickTree *tree;
    pcbr_workspace *w;
    double args[4] = {1.5, 1.8, 1.5, 1.8};
    double lik1, lik2;

    R_library("ape");

    rtree = call_R("rtree", rargs, tags, 1);
    rargs[0] = rtree;
    treestr = call_R("write.tree", rargs, tags, 1);
    tree = pllNewickParseString(CHAR(asChar(treestr)));
    w = pcbr_create(tree, 2);
    fill_parameters(w, &args[2]);
    
    lik1 = liklong(w, args, &args[2], 1, 2);
    lik2 = pow(10, likelihood(w, args, 0));
    ck_assert(fabs(lik1 - lik2) < 1e-6);

    UNPROTECT(1);
    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

START_TEST (test_likelihood_random_3)
{
    int ntip = 6;
    SEXP rargs[] = { PROTECT(ScalarInteger(ntip)) }, rtree, treestr;
    const char *tags[] = { NULL };
    pllNewickTree *tree;
    pcbr_workspace *w;
    double args[9];
    int i;
    double lik1, lik2;

    for (i = 0; i < 9; ++i)
        args[i] = (double) rand() / (double) RAND_MAX;

    R_library("ape");

    rtree = call_R("rtree", rargs, tags, 1);
    rargs[0] = rtree;
    treestr = call_R("write.tree", rargs, tags, 1);
    tree = pllNewickParseString(CHAR(asChar(treestr)));
    w = pcbr_create(tree, 3);
    fill_parameters(w, &args[3]);

    lik1 = liklong(w, args, &args[3], 1, 3);
    lik2 = pow(10, likelihood(w, args, 0));

    ck_assert(fabs(lik1 - lik2) < 1e-6);

    UNPROTECT(1);
    pcbr_free(w);
    pllNewickParseDestroy(&tree);
}
END_TEST

Suite *likelihood_suite(void)
{
    Suite *s;
    TCase *tc_guess, *tc_trees, *tc_lik, *tc_fuzz;

    s = suite_create("likelihood");

    tc_trees = tcase_create("Trees");
    tcase_add_test(tc_trees, test_flatten);
    tcase_add_test(tc_trees, test_scale_branches);
    suite_add_tcase(s, tc_trees);

    tc_guess = tcase_create("Guessing");
    tcase_add_test(tc_guess, test_guess_means);
    tcase_add_test(tc_guess, test_guess_parameters);
    suite_add_tcase(s, tc_guess);

    tc_lik = tcase_create("Likelihood");
    tcase_add_test(tc_lik, test_liklong);
    tcase_add_test(tc_lik, test_liklong_3);
    tcase_add_test(tc_lik, test_pi);
    tcase_add_test(tc_lik, test_pmut);
    tcase_add_test(tc_lik, test_pmut_3);
    tcase_add_test(tc_lik, test_pt_nonterminal);
    tcase_add_test(tc_lik, test_pt_nonterminal_3);
    tcase_add_test(tc_lik, test_pt_terminal);
    tcase_add_test(tc_lik, test_pt_terminal_3);
    tcase_add_test(tc_lik, test_likelihood_terminals);
    tcase_add_test(tc_lik, test_likelihood_3);
    tcase_add_test(tc_lik, test_reconstruct_terminals);
    tcase_add_test(tc_lik, test_reconstruct_3);
    tcase_add_test(tc_lik, test_kb);
    tcase_add_test(tc_lik, test_kb_large);
    tcase_set_timeout(tc_lik, 10);
    suite_add_tcase(s, tc_lik);

    tc_fuzz = tcase_create("Fuzz");
    tcase_add_test(tc_lik, test_likelihood_random);
    tcase_add_test(tc_lik, test_likelihood_random_3);
    suite_add_tcase(s, tc_fuzz);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    setenv("R_HOME", R_HOME, 0);

    s = likelihood_suite();
    sr = srunner_create(s);

    start_R();
    srunner_run_all(sr, CK_NORMAL);
    stop_R();
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return number_failed;
}
