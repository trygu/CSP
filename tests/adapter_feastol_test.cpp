// Feasibility tolerance contract for the network LPs.
//
// The per-level network LPs built by the capacity separator contain
// redundant equations: the same totals appear as row right-hand sides
// that were computed in different floating-point addition orders, so
// the system can be infeasible by ~1e-6 (absolute) from rounding
// alone. The legacy CPLEX backend set no feasibility tolerance, so
// it ran at the CPLEX 7.5 default (EPRHS = 1e-6), whose row scaling
// effectively accepted these models; HiGHS at 1e-6 rejects them and
// the run fails at the first separator level, so the adapter runs at
// 1e-4 to restore the legacy effective acceptance (documented
// deviation, orchestration docs/d2.6-alignment.md).
//
// The values below are a real failing dump (sdcnet.lp from a modular
// batch on HiGHS 1.15.1): the c1 implied by rows r0 and r3 differs by
// ~1.05e-6. Also checked: the network-solve path reports a failed
// solve as a nonzero netopt result (1101).

#include <highs/interfaces/highs_c_api.h>
#include "Jjsolver.h"

#include <cstdio>

int JJtest_failNextNetLp(void);

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

int main(void)
{
    const int mac = 3, mar = 4;
    double o[3]    = { 0.0, 1.0, 0.0 }; // max c1
    double rhs[4]  = { 6940067.14999956,
                       11388498.2749995,
                       6160132.41999991,
                       12168433.0050002 };
    char sense[4]  = { '=', '=', '=', '=' };
    int matbeg[3]  = { 0, 2, 4 };
    int matcnt[3]  = { 2, 2, 2 };
    int matind[6]  = { 0, 2, 0, 3, 1, 3 }; // c0: r0,r2; c1: r0,r3; c2: r1,r3
    double matval[6] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    double lb[3]    = { 0.0, 0.0, 0.0 };
    double ub[3]    = { 13565258231.0, 13565258231.0, 13565258231.0 };

    JJLPptr lp = JJloadprob(
        (char*)"t", mac, mar, 0, 1,
        o, rhs, sense,
        matbeg, matcnt, matind, matval, lb, ub, (double*)0,
        (int*)0, (int*)0, (int*)0, (int*)0, (int*)0, (double*)0,
        (char*)"d", (char*)"o", (char*)"r", (char*)"n", (char*)"b",
        (char**)0, (char*)"cns", (char**)0, (char*)"rns", (char**)0,
        (char*)"nns", 0, 0, 0, 0, 0, 0, 0, 0);
    CHECK(lp != (JJLPptr)0);

    int rc = JJdualopt(lp);
    CHECK(rc == 0);
    CHECK(JJgetstat(lp) == JJ_OPTIMAL);

    /* The protection-level network path must surface a failed solve
       as a nonzero netopt result (the host records CSPexit
       EXIT_LPSOLVER on that; Cspnet.c). */
    JJtest_failNextNetLp();
    int netstatus = 0, nnodes = 0, narcs = 0, iters = 0;
    rc = JJnetopt(lp, &netstatus, &nnodes, &narcs, &iters);
    CHECK(rc != 0);
    CHECK(netstatus == 1101);

    JJfreeprob(&lp);

    if (failures == 0) {
        std::printf("adapter_feastol_test: all checks passed\n");
        return 0;
    }
    std::printf("adapter_feastol_test: %d failures\n", failures);
    return 1;
}
