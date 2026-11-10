// JJlpiterlimit test (VSCIP / HiGHS port).
//
// Solves go through JJoptimize, the library entry every CSP call
// path uses; the iteration cap is applied there (and in the dual
// path), matching where the legacy CPLEX parameter took effect.
// Ground truth is the HiGHS model status after the solve. The legacy
// CPLEX backend exposed CPXsetintparam(CPX_PARAM_ITERATIONS, val)
// through JJlpiterlimit(int), which applied to every subsequent solve;
// the heuristic's constructive phase (cspheur.c) sets the cap to 1000
// before the protection-level network LPs and resets it to 1000000
// after, and accepts the current bound as the level when the solver
// stops at the limit (JJgetstat == 10, cspheur.c:649). The port must
// forward the cap to the solver (HiGHS simplex_iteration_limit):
// without it, HiGHS drives the degenerate network LPs to full
// optimality and OPT hangs on large tables.
//
// Covers:
//   - no cap: the solve reaches kHighsModelStatusOptimal
//   - cap 1: a fresh solve stops at kHighsModelStatusIterationLimit
//     (JJgetstat 10), not optimal
//   - the cap persists across solves (not a one-shot setting)
//   - reset to 1000000: optimal again

#include <highs/interfaces/highs_c_api.h>
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

static HighsInt hstat(JJLPptr lp)
{
    return Highs_getModelStatus(lp);
}

// max x0 + x1 + x2 s.t. x0 + x1 <= 1, x1 + x2 <= 1, 0 <= x <= 1
// (objsen -1: HiGHS 1 minimizes, -1 maximizes). The optimum
// (0.5, 0.5, 0.5) is not reachable in a single simplex pivot from
// the zero starting basis, so a 1-iteration cap cannot reach
// optimality. Presolve is off (string option: HiGHS 1.15 rejects a
// bool) so the solve starts from the standard initial basis instead
// of a postsolve basis.
static JJLPptr load_lp(void)
{
    const int mac = 3, mar = 2;
    double o[3] = { 1.0, 1.0, 1.0 };
    double rhs[2] = { 1.0, 1.0 };
    char sense[2] = { 'L', 'L' };
    int matbeg[3] = { 0, 1, 3 };
    int matcnt[3] = { 1, 2, 1 };
    int matind[4] = { 0, 0, 1, 1 };
    double matval[4] = { 1.0, 1.0, 1.0, 1.0 };
    double lb[3] = { 0.0, 0.0, 0.0 };
    double ub[3] = { 1.0, 1.0, 1.0 };
    JJLPptr lp = JJloadprob(
        (char*)"t", mac, mar, 0, -1,
        o, rhs, sense,
        matbeg, matcnt, matind, matval, lb, ub, (double*)0,
        (int*)0, (int*)0, (int*)0, (int*)0, (int*)0, (double*)0,
        (char*)"d", (char*)"o", (char*)"r", (char*)"n", (char*)"b",
        (char**)0, (char*)"cns", (char**)0, (char*)"rns", (char**)0,
        (char*)"nns", 0, 0, 0, 0, 0, 0, 0, 0);
    if (lp != (JJLPptr)0)
        Highs_setStringOptionValue(lp, "presolve", "off");
    return lp;
}

int main()
{
    setvbuf(stdout, 0, _IONBF, 0);

    // Baseline: no cap set, the solve reaches optimal.
    {
        JJLPptr lp = load_lp();
        CHECK(lp != (JJLPptr)0);
        CHECK(JJoptimize(lp) == 0);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJfreeprob(&lp);
    }

    // Cap of 1 iteration: a fresh solve stops at the limit.
    CHECK(JJlpiterlimit(1) == 0);
    {
        JJLPptr lp = load_lp();
        CHECK(lp != (JJLPptr)0);
        CHECK(JJoptimize(lp) == 0);
        CHECK(hstat(lp) == kHighsModelStatusIterationLimit);
        CHECK(JJgetstat(lp) == 10);
        JJfreeprob(&lp);
    }

    // The cap persists across solves (legacy CPXsetintparam contract:
    // it is not a one-shot setting).
    {
        JJLPptr lp = load_lp();
        CHECK(lp != (JJLPptr)0);
        CHECK(JJoptimize(lp) == 0);
        CHECK(hstat(lp) == kHighsModelStatusIterationLimit);
        JJfreeprob(&lp);
    }

    // Reset to the heuristic's post-constructive value: optimal again.
    CHECK(JJlpiterlimit(1000000) == 0);
    {
        JJLPptr lp = load_lp();
        CHECK(lp != (JJLPptr)0);
        CHECK(JJoptimize(lp) == 0);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJfreeprob(&lp);
    }

    if (failures) {
        std::printf("%d failures\n", failures);
        return 1;
    }
    std::printf("ok\n");
    return 0;
}
