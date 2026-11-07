// Adapter regression test: the legacy solver contract (CPLEX-era
// Jjsolver.c) treats a definitive Infeasible or Unbounded model
// status as a valid solver outcome, not a failure. JJdualopt and
// JJoptimize return 0 and the caller reads the status through
// JJgetstat; solve_lp's status switch is the designated handler for
// those cases (retry pricing, then close the node at the upper
// bound). If the solve functions turn those statuses into failures,
// the switch becomes unreachable and the exit propagation added in
// PR #43 turns a legitimate infeasible node LP into error 7002.
//
// Ground truth is the legacy CPLEX-era implementation, which returned
// the solver status directly (XPRSminim / CPXdualopt / SCIP_OKAY)
// and only signalled an error when the solver call itself failed.

#include "Jjsolver.h"

#include <cstdio>

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

// min c0  s.t.  c0 >= 2,  0 <= c0 <= 1.  Infeasible by construction.
static JJLPptr load_infeasible_lp()
{
    const int mac = 1, mar = 1;
    double o[1] = { 1.0 };
    double rhs[1] = { 2.0 };
    char sense[1] = { 'G' };
    int matbeg[1] = { 0 };
    int matcnt[1] = { 1 };
    int matind[1] = { 0 };
    double matval[1] = { 1.0 };
    double lb[1] = { 0.0 };
    double ub[1] = { 1.0 };
    return JJloadprob(
        (char*)"t", mac, mar, 0, 1,
        o, rhs, sense,
        matbeg, matcnt, matind, matval, lb, ub, (double*)0,
        (int*)0, (int*)0, (int*)0, (int*)0, (int*)0, (double*)0,
        (char*)"d", (char*)"o", (char*)"r", (char*)"n", (char*)"b",
        (char**)0, (char*)"cns", (char**)0, (char*)"rns", (char**)0,
        (char*)"nns", 0, 0, 0, 0, 0, 0, 0, 0);
}

// min -c0  s.t.  c0 >= 0,  0 <= c0 <= 1e20 (legacy INFBOUND, which
// HiGHS treats as an infinite bound).  Unbounded below.
static JJLPptr load_unbounded_lp()
{
    const int mac = 1, mar = 1;
    double o[1] = { -1.0 };
    double rhs[1] = { 0.0 };
    char sense[1] = { 'G' };
    int matbeg[1] = { 0 };
    int matcnt[1] = { 1 };
    int matind[1] = { 0 };
    double matval[1] = { 1.0 };
    double lb[1] = { 0.0 };
    double ub[1] = { 1.0E+20 };
    return JJloadprob(
        (char*)"t", mac, mar, 0, 1,
        o, rhs, sense,
        matbeg, matcnt, matind, matval, lb, ub, (double*)0,
        (int*)0, (int*)0, (int*)0, (int*)0, (int*)0, (double*)0,
        (char*)"d", (char*)"o", (char*)"r", (char*)"n", (char*)"b",
        (char**)0, (char*)"cns", (char**)0, (char*)"rns", (char**)0,
        (char*)"nns", 0, 0, 0, 0, 0, 0, 0, 0);
}

int main()
{
    // Dual path: the solve must not report a failure, and the status
    // must be readable as the definitive Infeasible result.
    JJLPptr lp = load_infeasible_lp();
    CHECK(lp != (JJLPptr)0);
    if (lp == (JJLPptr)0) return 1;
    CHECK(JJdualopt(lp) == 0);
    CHECK(JJgetstat(lp) == JJ_INFEASIBLE);
    JJfreeprob(&lp);

    // Primal path: same contract.
    JJLPptr lp2 = load_infeasible_lp();
    CHECK(lp2 != (JJLPptr)0);
    if (lp2 == (JJLPptr)0) return 1;
    CHECK(JJoptimize(lp2) == 0);
    CHECK(JJgetstat(lp2) == JJ_INFEASIBLE);
    JJfreeprob(&lp2);

    // Unbounded models are a valid outcome on both paths as well.
    JJLPptr lp3 = load_unbounded_lp();
    CHECK(lp3 != (JJLPptr)0);
    if (lp3 == (JJLPptr)0) return 1;
    CHECK(JJdualopt(lp3) == 0);
    CHECK(JJgetstat(lp3) == JJ_UNBOUNDED);
    JJfreeprob(&lp3);

    JJLPptr lp4 = load_unbounded_lp();
    CHECK(lp4 != (JJLPptr)0);
    if (lp4 == (JJLPptr)0) return 1;
    CHECK(JJoptimize(lp4) == 0);
    CHECK(JJgetstat(lp4) == JJ_UNBOUNDED);
    JJfreeprob(&lp4);

    return failures ? 1 : 0;
}
