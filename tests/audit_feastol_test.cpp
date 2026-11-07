// Feasibility tolerance contract for the audit LP.
//
// Real tables carry relations whose right-hand sides are the same
// totals computed in different floating-point addition orders, so the
// table is mutually inconsistent by ~1e-6 (absolute) from rounding
// alone. The legacy audit (intervalle.exe) accepted such tables;
// HiGHS at its 1e-6 default declares the audit LP infeasible and the
// whole audit aborts (rc -1) on tables the legacy audit accepted, so
// the audit runs at 1e-4 to restore the legacy effective acceptance
// (documented deviation, orchestration docs/d2.6-alignment.md).
//
// The numbers below mirror a real failing dump (sdcnet.lp from a
// modular batch on HiGHS 1.15.1): the implied value of the suppressed
// cell differs by ~1.05e-6 between two overlapping relations.

#include "cspaudit.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

static int write_jj(const char* path, const char* body)
{
    FILE* f = fopen(path, "w");
    if (!f) {
        return -1;
    }
    fputs(body, f);
    fclose(f);
    return 0;
}

int main(void)
{
    const double v0 = 1000000.123456;   // suppressed cell value
    const double v1 = 5000000.54321;    // safe cell value
    const double a = v0 + v1;           // relation 0 rhs
    const double b = v0 + 1.05e-6;      // relation 1 rhs: implied v0 + 1.05e-6

    char path[256];
    char body[512];
    snprintf(body, sizeof(body),
             "0\n"
             "2\n"
             "0 %.17g 1 u %.17g %.17g 0.0 0.0 0.0\n"
             "1 %.17g 1 s %.17g %.17g 0.0 0.0 0.0\n"
             "2\n"
             "%.17g 2 : 0 (1) 1 (1)\n"
             "%.17g 1 : 0 (1)\n",
             v0, v0 - 1000.0, v0 + 1000.0,
             v1, v1 - 1000.0, v1 + 1000.0,
             a, b);
    snprintf(path, sizeof(path), "audit_feastol.jj");
    CHECK(write_jj(path, body) == 0);

    int n = 0;
    long* idx = 0;
    double* mn = 0;
    double* mx = 0;
    double* val = 0;
    char* st = 0;
    int* unsafe = 0;

    int rc = TauAuditJj(path, &n, &idx, &mn, &mx, &val, &st, &unsafe);
    remove(path);
    CHECK(rc == 0);
    if (rc == 0) {
        CHECK(n == 1);
        if (n == 1) {
            CHECK(idx[0] == 0);
            CHECK(st[0] == 'u');
            CHECK(std::fabs(val[0] - v0) < 1e-9);
            // The LP may violate rows by up to the 1e-4 tolerance.
            CHECK(std::fabs(mn[0] - v0) < 1e-3);
            CHECK(std::fabs(mx[0] - b) < 1e-3);
            CHECK(unsafe[0] == 0);
        }
        TauAuditJjFree(n, idx, mn, mx, val, st, unsafe);
    }

    if (failures == 0) {
        std::printf("audit_feastol_test: all checks passed\n");
        return 0;
    }
    std::printf("audit_feastol_test: %d failures\n", failures);
    return 1;
}
