// Solver log lifecycle: the CSP log file must not outlive the solve.
//
// PPCSPloadprob opens the global filestr on fCSPlog and redirects the
// process std::cout into it. Neither was ever closed nor restored, so
//   1) the file stays locked on Windows for the process lifetime
//      (front-end cleanup of the solver temp dir failed with
//      PermissionError), and
//   2) after the first MIP solve, every module's std::cout output
//      silently went into the log file instead of the process stdout.
//
// Contract under test: after a full
// loadprob / optimize / solution / freeprob cycle, std::cout is back on
// its original streambuf, and a second solve in the same process behaves
// identically: the restore must be idempotent, not "restore the log's
// own buffer" (which the pre-fix code did on its second backup capture).

#include <cstdio>
#include <fstream>

#include "cspback.h"
#include "cspmain.h"
#include "resource.h"

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

// The 4-cell, 1-total table that solves cleanly (status 'u' is the
// sensitive cell), with the HiTaSCtrl production defaults. Same problem
// as csp_exit_test.
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

// One full solve through the production call sequence.
static void solve_once(void)
{
    int rc = SCIPv::CSPoptimize(NULL);
    CHECK(rc == 0);
    int lcost = 0, ucost = 0;
    char st[4] = { 'u', 's', 's', 's' };
    CHECK(SCIPv::CSPsolution(&lcost, &ucost, st) == 0);
    CHECK(lcost <= ucost);
    SCIPv::CSPfreeprob();
    SCIPv::CSPFreeFileNames();
}

int main(void)
{
    std::streambuf* orig = std::cout.rdbuf();

    load_small_problem();
    solve_once();

    // The solve must have gone through the log path (CSPSetFileNames("")
    // puts fCSPlog in the working directory, as in the legacy contract).
    CHECK(std::ifstream("CSPlogfile.txt").good());
    // The redirect is undone when the solve ends.
    CHECK(std::cout.rdbuf() == orig);

    // A second solve in the same process: the restore above must be
    // idempotent, so the backup captured here is still the real stdout.
    load_small_problem();
    solve_once();
    CHECK(std::cout.rdbuf() == orig);

    if (failures == 0) {
        std::printf("csp_log_lifecycle_test: all checks passed\n");
        return 0;
    }
    std::printf("csp_log_lifecycle_test: %d failures\n", failures);
    return 1;
}
