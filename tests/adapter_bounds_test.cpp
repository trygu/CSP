// Adapter tests for the no-op bound skip (column and row bounds).
//
// The bound values are read back from the HiGHS C API as ground truth.
// Whether a change was forwarded to HiGHS or skipped as a no-op is
// observed through Highs_getModelStatus: after Highs_run the status is
// kHighsModelStatusOptimal; a forwarded change mutates the model and
// resets the status, a skipped no-op leaves it untouched.

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

// 3-column LP: max obj.x s.t. x0 + x1 + x2 <= 30, 0 <= xi <= 10.
// Feasible and bounded, so Highs_run reaches kHighsModelStatusOptimal.
static JJLPptr load_lp(const double obj[3])
{
    const int mac = 3, mar = 1;
    double o[3] = { obj[0], obj[1], obj[2] };
    double rhs[1] = { 30.0 };
    char sense[1] = { 'L' };
    int matbeg[3] = { 0, 1, 2 };
    int matcnt[3] = { 1, 1, 1 };
    int matind[3] = { 0, 0, 0 };
    double matval[3] = { 1.0, 1.0, 1.0 };
    double lb[3] = { 0.0, 0.0, 0.0 };
    double ub[3] = { 10.0, 10.0, 10.0 };
    return JJloadprob(
        (char*)"t", mac, mar, 0, 1,
        o, rhs, sense,
        matbeg, matcnt, matind, matval, lb, ub, (double*)0,
        (int*)0, (int*)0, (int*)0, (int*)0, (int*)0, (double*)0,
        (char*)"d", (char*)"o", (char*)"r", (char*)"n", (char*)"b",
        (char**)0, (char*)"cns", (char**)0, (char*)"rns", (char**)0,
        (char*)"nns", 0, 0, 0, 0, 0, 0, 0, 0);
}

static void read_bounds(JJLPptr lp, int first, int last, double* lo, double* hi)
{
    HighsInt num_col = 0, num_nz = 0;
    Highs_getColsByRange(lp, first, last, &num_col, (double*)0, lo, hi,
                         &num_nz, (HighsInt*)0, (int*)0, (double*)0);
}

int main()
{
    JJLPptr lp = load_lp((const double[]){ 1.0, 2.0, 3.0 });
    CHECK(lp != (JJLPptr)0);
    int idx[1] = { 0 };
    char lu;
    double bd;
    double lo[3], hi[3];

    CHECK(Highs_run(lp) == kHighsStatusOk);
    CHECK(hstat(lp) == kHighsModelStatusOptimal);

    // Column lower bound: real change, then no-op.
    idx[0] = 0; lu = 'L'; bd = 1.0;
    JJchgbds(lp, 1, idx, &lu, &bd);
    CHECK(hstat(lp) != kHighsModelStatusOptimal);
    read_bounds(lp, 0, 2, lo, hi);
    CHECK(lo[0] == 1.0);
    CHECK(Highs_run(lp) == kHighsStatusOk);
    CHECK(hstat(lp) == kHighsModelStatusOptimal);
    JJchgbds(lp, 1, idx, &lu, &bd);   // no-op
    CHECK(hstat(lp) == kHighsModelStatusOptimal);

    // Column upper bound: real change, then no-op.
    idx[0] = 1; lu = 'U'; bd = 9.0;
    JJchgbds(lp, 1, idx, &lu, &bd);
    CHECK(hstat(lp) != kHighsModelStatusOptimal);
    read_bounds(lp, 0, 2, lo, hi);
    CHECK(hi[1] == 9.0);
    CHECK(Highs_run(lp) == kHighsStatusOk);
    CHECK(hstat(lp) == kHighsModelStatusOptimal);
    JJchgbds(lp, 1, idx, &lu, &bd);   // no-op
    CHECK(hstat(lp) == kHighsModelStatusOptimal);

    // Fixed column: real change, then no-op.
    idx[0] = 2; lu = 'B'; bd = 4.0;
    JJchgbds(lp, 1, idx, &lu, &bd);
    CHECK(hstat(lp) != kHighsModelStatusOptimal);
    read_bounds(lp, 0, 2, lo, hi);
    CHECK(lo[2] == 4.0 && hi[2] == 4.0);
    CHECK(Highs_run(lp) == kHighsStatusOk);
    CHECK(hstat(lp) == kHighsModelStatusOptimal);
    JJchgbds(lp, 1, idx, &lu, &bd);   // no-op
    CHECK(hstat(lp) == kHighsModelStatusOptimal);

    // Row bound through JJchgcoef(row, -1, value). The row is 'L' with
    // rhs 30, i.e. lo = -inf, hi = 30: real change, then no-op.
    JJchgcoef(lp, 0, -1, 25.0);
    CHECK(hstat(lp) != kHighsModelStatusOptimal);
    {
        double rlo[1], rhi[1];
        HighsInt nr = 0, nnz = 0;
        Highs_getRowsByRange(lp, 0, 0, &nr, rlo, rhi, &nnz,
                             (HighsInt*)0, (int*)0, (double*)0);
        CHECK(rhi[0] == 25.0);
    }
    CHECK(Highs_run(lp) == kHighsStatusOk);
    CHECK(hstat(lp) == kHighsModelStatusOptimal);
    JJchgcoef(lp, 0, -1, 25.0);       // no-op
    CHECK(hstat(lp) == kHighsModelStatusOptimal);

    JJfreeprob(&lp);

    if (failures) {
        std::printf("adapter_bounds: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("adapter_bounds: all checks passed\n");
    return 0;
}
