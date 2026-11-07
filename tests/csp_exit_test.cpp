// CSPexit error propagation.
//
// Legacy tau-ARGUS called exit(1) when the LP solver or the
// branch-and-bound hit a fatal error. The engine replaced exit(1) with
// CSPexit(), but no host ever registered an exit handler
// (CSPdefineexit), so the call was a no-op: the optimizer kept retrying
// a failed run forever, re-dumping the LP model on every attempt.
//
// The contract under test: CSPexit records the first error flag,
// PPCSPoptimize checks the flag at each node of the search and returns
// it (bounded, nonzero) instead of looping, and a clean run still
// returns 0 with no error recorded. The host (hitas) maps a nonzero
// CSPoptimize result to its failure codes.
//
// The third case fails the next DUAL LP solve via a test hook in
// Jjsolver.c. Dual solves only happen inside branch-and-bound (the
// initial heuristic uses primal solves), so the failure lands in the
// same place as the field bug. The test carries a ctest timeout: a
// run that retries the failed solve forever must be reported as a
// failure, not a pass.

#include "cspback.h"
#include "cspmain.h"
#include "resource.h"

#include <cstdio>

int JJtest_failNextDualLp(void);

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

// The 4-cell, 1-total table that FullJJ solves cleanly (status 'u' is
// the sensitive cell), with the HiTaSCtrl production defaults.
static void load_small_problem(void)
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
    double lpl[4]     = { 1, 0, 0, 0 };
    double upl[4]     = { 1, 0, 0, 0 };
    double lb[4]      = { 5, 3, 5, 1 };
    double ub[4]      = { 15, 7, 9, 5 };
    double rhs[1]     = { 25 };
    int    ncard[1]   = { 4 };
    int    list[4]    = { 0, 1, 2, 3 };
    signed char coef[4] = { 1, 1, 1, 1 };

    CHECK(SCIPv::CSPloadprob(nsums, rhs, ncells, data, weight, status,
                        lpl, upl, lb, ub, namelist, ncard, list, coef) == 0);
}

static void test_flag(void)
{
    CSPclearerror();
    CHECK(CSPerror() == 0);
    CSPexit(EXIT_ERROR);
    CHECK(CSPerror() == EXIT_ERROR);
    CSPexit(EXIT_LPSOLVER);
    CHECK(CSPerror() == EXIT_ERROR); // first error wins
    CSPclearerror();
    CHECK(CSPerror() == 0);
}

static void test_clean_run(void)
{
    load_small_problem();
    int rc = SCIPv::CSPoptimize(NULL);
    CHECK(rc == 0);
    CHECK(CSPerror() == 0);
    int lcost = 0, ucost = 0;
    /* CSPsolution fills the per-cell status array of the whole table
       (the production callers pass &(*states) / Tab.status). */
    char st[4] = { 'u', 's', 's', 's' };
    CHECK(SCIPv::CSPsolution(&lcost, &ucost, st) == 0);
    CHECK(lcost <= ucost);
    SCIPv::CSPfreeprob();
    SCIPv::CSPFreeFileNames();
}

static void test_forced_dual_lp_failure(void)
{
    load_small_problem();
    JJtest_failNextDualLp();
    int rc = SCIPv::CSPoptimize(NULL);
    CHECK(rc != 0);
    CHECK(rc == CSPerror());
    SCIPv::CSPfreeprob();
    SCIPv::CSPFreeFileNames();
}

int main(void)
{
    test_flag();
    test_clean_run();
    test_forced_dual_lp_failure();
    if (failures == 0) {
        std::printf("csp_exit_test: all checks passed\n");
        return 0;
    }
    std::printf("csp_exit_test: %d failures\n", failures);
    return 1;
}
