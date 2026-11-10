// Protection-level network LP feasibility tolerance.
//
// The per-level protection LP built by the heuristic (cspheur.c) has
// zero-rhs balance rows, so rounding can never make it inconsistent
// and it must run at the legacy CPLEX feasibility tolerance (EPRHS =
// 1e-6). JJloadprob sets 1e-4 for the data-valued models (the B&B
// node LP and the capacity-separator network LP, whose redundant
// relations are inconsistent by ~1e-6 on real tables; documented
// deviation, orchestration docs/d2.6-alignment.md). On the
// protection LP that loose tolerance is actively harmful: a forced
// protection flow below 1e-4 is accepted as a plain violation of the
// balance rows, the LP reports optimal with objective 0, no flow is
// rerouted through the other cells, and no secondary suppression is
// left to mark (the MOD D=0 regression: 0 secondary cells where
// legacy tau-ARGUS 4.1 produced 224).
//
// Contract under test: with a protection level below 1e-4, the final
// solution must protect the sensitive cell through at least one
// secondary suppression (a non-sensitive cell marked 'm' by
// CSPsolution), i.e. the forced level was rerouted through it.

#include "cspback.h"
#include "cspmain.h"
#include "resource.h"

#include <cstdio>

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

// The 4-cell, 1-total table that csp_exit_test solves cleanly, with
// the sensitive cell protected at level 1e-5 (above JJZERO = 1e-7,
// below the 1e-4 tolerance the buggy path ran the protection LP at).
static void load_level_problem(void)
{
    SCIPv::CSPSetFileNames("");
    SCIPv::CSPSetDoubleConstant(JJMAXTIME, 3600.0);
    SCIPv::CSPSetDoubleConstant(JJZERO, 1e-7);
    SCIPv::CSPSetDoubleConstant(JJINF, 2.14e9);
    SCIPv::CSPSetDoubleConstant(JJMINVIOLA, 1e-4);
    SCIPv::CSPSetDoubleConstant(JJMAXSLACK, 1e-4);
    SCIPv::CSPSetDoubleConstant(JJFEASTOL, 1e-6);
    SCIPv::CSPSetDoubleConstant(JJOPTTOL, 1e-9);
    SCIPv::CSPSetIntegerConstant(JJMAXCOLSLP, 50000);
    SCIPv::CSPSetIntegerConstant(JJMAXROWSLP, 15000);
    SCIPv::CSPSetIntegerConstant(JJMAXCUTSPOOL, 500000);
    SCIPv::CSPSetIntegerConstant(JJMAXCUTSITER, 50);

    const int ncells = 4;
    const int nsums = 1;
    static char names[4][4] = { "x0", "x1", "x2", "x3" };
    char *namelist[4] = { names[0], names[1], names[2], names[3] };
    double data[4]    = { 10, 5, 7, 3 };
    int    weight[4]  = { 1, 1, 1, 1 };
    char   status[4]  = { 'u', 's', 's', 's' };
    double lpl[4]     = { 1e-5, 0, 0, 0 };
    double upl[4]     = { 1e-5, 0, 0, 0 };
    double lb[4]      = { 5, 3, 5, 1 };
    double ub[4]      = { 15, 7, 9, 5 };
    double rhs[1]     = { 25 };
    int    ncard[1]   = { 4 };
    int    list[4]    = { 0, 1, 2, 3 };
    signed char coef[4] = { 1, 1, 1, 1 };

    CHECK(SCIPv::CSPloadprob(nsums, rhs, ncells, data, weight, status,
                        lpl, upl, lb, ub, namelist, ncard, list, coef) == 0);
}

int main(void)
{
    load_level_problem();
    int rc = SCIPv::CSPoptimize(NULL);
    CHECK(rc == 0);
    CHECK(CSPerror() == 0);
    int lcost = 0, ucost = 0;
    char st[4] = { 'u', 's', 's', 's' };
    CHECK(SCIPv::CSPsolution(&lcost, &ucost, st) == 0);
    CHECK(lcost <= ucost);
    int secondary = (st[1] == 'm') + (st[2] == 'm') + (st[3] == 'm');
    /* The forced 1e-5 flow must be rerouted through the other cells
       of the total; the rerouted cell(s) are the secondary
       suppressions. With the 1e-4 tolerance the flow vanishes, no
       cell is rerouted, and no secondary exists. */
    CHECK(secondary >= 1);
    SCIPv::CSPfreeprob();
    SCIPv::CSPFreeFileNames();
    if (failures == 0) {
        std::printf("protection_flow_test: all checks passed\n");
        return 0;
    }
    std::printf("protection_flow_test: %d failures\n", failures);
    return 1;
}
