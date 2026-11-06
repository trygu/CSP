// TauAuditJj on fractional tables.
//
// The legacy intervalle program Int() cast the JJ bounds. On integer
// tables the cast is a no-op (bounds are integer-valued), but on
// fractional tables it clips a bound below the value an equality
// margin demands, so the audit LP is infeasible and TauAuditJj fails
// with rc -1. The audit must use the JJ bounds as-is, which keeps the
// integer behaviour exactly the same.

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

// Margin identity -x0 + x1 = 0 with a suppressed x0 whose upper bound
// is its own value. "frac" uses 10.5, "int" uses 10.
static void run_case(const char* tag, const char* cell_value, double expected)
{
    char path[256];
    char body[512];
    int n = 0;
    long* idx = 0;
    double* mn = 0;
    double* mx = 0;
    double* val = 0;
    char* st = 0;
    int* unsafe = 0;

    snprintf(body, sizeof(body),
             "0\n"
             "2\n"
             "0 %s 1 u 0.0 %s 0.0 0.0 0.0\n"
             "1 %s 1 s 0.0 %s 0.0 0.0 0.0\n"
             "1\n"
             "0 2 : 0 (-1) 1 (1)\n",
             cell_value, cell_value, cell_value, cell_value);
    snprintf(path, sizeof(path), "audit_%s.jj", tag);
    if (write_jj(path, body) != 0) {
        CHECK(false && "write_jj");
        return;
    }

    int rc = TauAuditJj(path, &n, &idx, &mn, &mx, &val, &st, &unsafe);
    CHECK(rc == 0);
    if (rc == 0) {
        CHECK(n == 1);
        if (n == 1) {
            CHECK(idx[0] == 0);
            CHECK(std::fabs(val[0] - expected) < 1e-9);
            CHECK(std::fabs(mn[0] - expected) < 1e-6);
            CHECK(std::fabs(mx[0] - expected) < 1e-6);
            CHECK(unsafe[0] == 0);
        }
        TauAuditJjFree(n, idx, mn, mx, val, st, unsafe);
    }
    remove(path);
}

int main()
{
    // Fractional: x0 = 10.5, hi = 10.5. The legacy Int() cast clips
    // the bound to 10 and the LP is infeasible (rc -1).
    run_case("frac", "10.5", 10.5);

    // Integer regression guard: x0 = 10, hi = 10. Behaviour must be
    // identical to the legacy Int() path.
    run_case("int", "10", 10.0);

    if (failures == 0) {
        std::printf("PASS audit_fractional\n");
    }
    return failures ? 1 : 0;
}
