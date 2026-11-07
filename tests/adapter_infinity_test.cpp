// Adapter regression test: the legacy B&B LP (cspsolve.c) carries a
// dummy column 0 with objective 1.0E+20, and a B&B branch can force
// that column positive. HiGHS's default infinite_cost is exactly 1e20
// and treats any cost with magnitude greater than or equal to that as
// infinite, so at the default the dummy objective is infinite: the
// objective value comes out +inf (and the fuller node LPs report
// UnboundedOrInfeasible, surfaced as error 7002). The solver must keep
// 1.0E+20 a finite coefficient so the objective stays finite, matching
// the CPLEX backend the port replaced.
//
// Ground truth is the HiGHS C API on the JJLPptr: with the dummy forced
// positive, the LP must solve Optimal with a finite objective of exactly
// 1.0E+20. The sense is minimize (objsen 1), matching cspsolve.c.

#include <highs/interfaces/highs_c_api.h>
#include "Jjsolver.h"

#include <cmath>
#include <cstdio>

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

static HighsInt hstat(JJLPptr lp)
{
    return Highs_getModelStatus(lp);
}

// The shape of the initial B&B LP in cspsolve.c: one dummy column with
// objective 1.0E+20, plus the row x0 >= 0. A B&B branch forces the dummy
// to 1 (lower bound 1), so the objective is non-trivial. Sense minimize.
static JJLPptr load_dummy_lp()
{
    const int mac = 1, mar = 1;
    double o[1] = { 1.0E+20 };
    double rhs[1] = { 0.0 };
    char sense[1] = { 'G' };
    int matbeg[1] = { 0 };
    int matcnt[1] = { 1 };
    int matind[1] = { 0 };
    double matval[1] = { 1.0 };
    double lb[1] = { 1.0 };
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

int main()
{
    JJLPptr lp = load_dummy_lp();
    CHECK(lp != (JJLPptr)0);
    if (lp == (JJLPptr)0) return 1;

    // The LP solves to its exact finite optimum: the forced-positive
    // dummy contributes 1.0E+20 to the (minimized) objective.
    CHECK(Highs_run(lp) == kHighsStatusOk);
    CHECK(hstat(lp) == kHighsModelStatusOptimal);
    const double obj = Highs_getObjectiveValue(lp);
    CHECK(std::isfinite(obj));
    CHECK(obj == 1.0E+20);

    JJfreeprob(&lp);
    return failures ? 1 : 0;
}
