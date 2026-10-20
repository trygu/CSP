// Adapter lifecycle tests for the Jjsolver objective cost shadow
// (VSCIP / HiGHS port).
//
// Ground truth is the HiGHS C API on the JJLPptr. Whether a change was
// forwarded to HiGHS or skipped as a no-op is observed through
// Highs_getModelStatus: after Highs_run the status is
// kHighsModelStatusOptimal; a forwarded change mutates the model and
// resets the status, a skipped no-op leaves it untouched.
//
// Covers:
//   - load: the shadow is installed from the passed objective
//   - no-op cost change: skipped (solved state preserved)
//   - real cost change: forwarded (solved state reset)
//   - JJaddcols: the shadow is extended with the new column costs
//   - JJaddrows(ccnt > 0): the shadow is extended with zero costs
//   - JJdelcols: the shadow is compacted (column index shift)
//   - drift: a desynchronized shadow is dropped on delete, and
//     changes are then forwarded (never mis-skipped)
//   - free/recreate: the shadow is erased on free and reinstalled on
//     load (looped to raise the chance of pointer reuse)
//   - JJchgobj: the batch path skips/applies per element

#include <highs/interfaces/highs_c_api.h>
#include "Jjsolver.h"

#include <cstdio>

// Test-only hooks, compiled into the library under TAUARGUS_CSP_UNIT_TEST.
void JJtest_corruptCostShadow(JJLPptr lp);
long JJtest_shadowCount(void);

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

static void read_costs(JJLPptr lp, int first, int last, double* out)
{
    HighsInt num_col = 0, num_nz = 0;
    Highs_getColsByRange(lp, first, last, &num_col, out,
                         (double*)0, (double*)0,
                         &num_nz, (HighsInt*)0, (int*)0, (double*)0);
}

int main()
{
    CHECK(JJtest_shadowCount() == 0);

    // T1-T6: one LP through the full edit lifecycle.
    {
        JJLPptr lp = load_lp((const double[]){ 1.0, 2.0, 3.0 });
        CHECK(lp != (JJLPptr)0);
        CHECK(JJtest_shadowCount() == 1);
        double c[8];
        read_costs(lp, 0, 2, c);
        CHECK(c[0] == 1.0 && c[1] == 2.0 && c[2] == 3.0);
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);

        // T2: no-op cost change is skipped: solved state preserved.
        JJchgcoef(lp, -1, 1, 2.0);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);

        // T3: a real change is forwarded: solved state reset.
        JJchgcoef(lp, -1, 1, 2.5);
        CHECK(hstat(lp) != kHighsModelStatusOptimal);
        read_costs(lp, 0, 2, c);
        CHECK(c[0] == 1.0 && c[1] == 2.5 && c[2] == 3.0);

        // T4: JJaddcols extends the shadow with the new column costs.
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        double nc[2] = { 10.0, 11.0 };
        double nlb[2] = { 0.0, 0.0 }, nub[2] = { 10.0, 10.0 };
        int nbeg[3] = { 0, 0, 0 };
        CHECK(JJaddcols(lp, 2, 0, nc, nbeg, (int*)0, (double*)0,
                        nlb, nub, (char**)0) == 0);
        read_costs(lp, 0, 4, c);
        CHECK(c[0] == 1.0 && c[1] == 2.5 && c[2] == 3.0
              && c[3] == 10.0 && c[4] == 11.0);
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgcoef(lp, -1, 4, 11.0);   // no-op on the appended column
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgcoef(lp, -1, 4, 12.0);   // real change
        CHECK(hstat(lp) != kHighsModelStatusOptimal);
        read_costs(lp, 4, 4, c);
        CHECK(c[0] == 12.0);

        // T4b: JJaddrows(ccnt > 0) appends zero-cost columns.
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        {
            double rhs2[1] = { 100.0 };
            char sense2[1] = { 'L' };
            int rbeg[2] = { 0, 2 };
            int rind[2] = { 5, 6 };   // the two appended columns
            double rval[2] = { 1.0, 1.0 };
            CHECK(JJaddrows(lp, 2, 1, 2, rhs2, sense2,
                            rbeg, rind, rval, (char**)0, (char**)0) == 0);
        }
        read_costs(lp, 5, 6, c);
        CHECK(c[0] == 0.0 && c[1] == 0.0);
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgcoef(lp, -1, 5, 0.0);    // no-op on a zero-cost column
        CHECK(hstat(lp) == kHighsModelStatusOptimal);

        // T5: JJdelcols compacts the shadow (column index shift).
        // Columns are [1, 2.5, 3, 10, 12, 0, 0].
        CHECK(JJdelcols(lp, 2, 2) == 0);   // drop the cost-3 column
        read_costs(lp, 0, 5, c);
        CHECK(c[0] == 1.0 && c[1] == 2.5 && c[2] == 10.0 && c[3] == 12.0
              && c[4] == 0.0 && c[5] == 0.0);
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgcoef(lp, -1, 3, 12.0);   // no-op at the shifted index
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgcoef(lp, -1, 1, 7.0);    // real change
        CHECK(hstat(lp) != kHighsModelStatusOptimal);
        read_costs(lp, 1, 1, c);
        CHECK(c[0] == 7.0);

        // T6: drift fallback. A desynchronized shadow is dropped on the
        // next delete; a later change with the model's current value must
        // be forwarded (a stale shadow would mis-skip it).
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJtest_corruptCostShadow(lp);          // shadow size 7 != mac 6
        CHECK(JJdelcols(lp, 5, 5) == 0);       // drops the shadow
        read_costs(lp, 0, 4, c);
        CHECK(c[0] == 1.0 && c[1] == 7.0 && c[2] == 10.0
              && c[3] == 12.0 && c[4] == 0.0);
        JJchgcoef(lp, -1, 1, 7.0);
        CHECK(hstat(lp) != kHighsModelStatusOptimal);

        JJfreeprob(&lp);
        CHECK(JJtest_shadowCount() == 0);
    }

    // T7: free/recreate. The shadow is erased on free and reinstalled on
    // load; looped so the allocator can reuse the previous address.
    for (int it = 0; it < 64; it++) {
        JJLPptr lp = load_lp((const double[]){ 5.0, 5.0, 5.0 });
        CHECK(JJtest_shadowCount() == 1);
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgcoef(lp, -1, 0, 5.0);   // no-op against the fresh shadow
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgcoef(lp, -1, 0, 9.9);   // real change; shadow now 9.9
        CHECK(hstat(lp) != kHighsModelStatusOptimal);
        JJfreeprob(&lp);
        CHECK(JJtest_shadowCount() == 0);
    }
    {
        JJLPptr lp = load_lp((const double[]){ 5.0, 5.0, 5.0 });
        CHECK(Highs_run(lp) == kHighsStatusOk);
        double c[1];
        read_costs(lp, 0, 0, c);
        CHECK(c[0] == 5.0);
        JJchgcoef(lp, -1, 0, 9.9);   // applied, not mis-skipped
        read_costs(lp, 0, 0, c);
        CHECK(c[0] == 9.9);
        JJfreeprob(&lp);
    }

    // T8: JJchgobj batch path: skip and apply are per element.
    {
        JJLPptr lp = load_lp((const double[]){ 1.0, 2.0, 3.0 });
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        int idx[2] = { 0, 1 };
        double obj[2] = { 1.0, 9.0 };    // idx 0 no-op, idx 1 real
        JJchgobj(lp, 2, idx, obj);
        CHECK(hstat(lp) != kHighsModelStatusOptimal);
        double c[3];
        read_costs(lp, 0, 2, c);
        CHECK(c[0] == 1.0 && c[1] == 9.0);
        CHECK(Highs_run(lp) == kHighsStatusOk);
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJchgobj(lp, 2, idx, obj);       // now a full no-op batch
        CHECK(hstat(lp) == kHighsModelStatusOptimal);
        JJfreeprob(&lp);
    }

    if (failures) {
        std::printf("adapter_lifecycle: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("adapter_lifecycle: all checks passed\n");
    return 0;
}
