// free_col/unfree_col on a cell that is not in the network LP.
//
// cell2net maps a cell to its NLP column (0-based); -1 means the cell
// is not part of this network. The legacy -1 guard was inside
// #ifdef STAMP only, so release builds passed -1 as a column index to
// the JJ shim and from there to HiGHS (out of range). free_col and
// unfree_col must be a clean no-op for out-of-network cells: no
// out-of-range call, no bound change, no clobbering of the saved
// bound buffer.
//
// The test-only hook JJtest_netFreeColState points the network state
// at a small LP and a cell mapping, so the pair can be exercised
// without a full SDC problem.

#include <highs/interfaces/highs_c_api.h>
#include "Jjsolver.h"

#include <cstdio>

// Plain C++ declarations: the library is compiled as C++.
void free_col(int index, double *bd);
void unfree_col(int index, double *bd);
int JJtest_netFreeColState(JJLPptr lp, int ncells, const int *mapping);

static int failures = 0;

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            ++failures; \
        } \
    } while (0)

static void read_bounds(JJLPptr lp, int first, int last, double* lo, double* hi)
{
    HighsInt num_col = 0, num_nz = 0;
    Highs_getColsByRange(lp, first, last, &num_col, (double*)0, lo, hi,
                         &num_nz, (HighsInt*)0, (int*)0, (double*)0);
}

// 3-column LP: max obj.x s.t. x0 + x1 + x2 <= 30, 0 <= xi <= 10.
static JJLPptr load_lp()
{
    const int mac = 3, mar = 1;
    double obj[3] = { 1.0, 2.0, 3.0 };
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
        obj, rhs, sense,
        matbeg, matcnt, matind, matval, lb, ub, (double*)0,
        (int*)0, (int*)0, (int*)0, (int*)0, (int*)0, (double*)0,
        (char*)"d", (char*)"o", (char*)"r", (char*)"n", (char*)"b",
        (char**)0, (char*)"cns", (char**)0, (char*)"rns", (char**)0,
        (char*)"nns", 0, 0, 0, 0, 0, 0, 0, 0);
}

int main()
{
    JJLPptr lp = load_lp();
    CHECK(lp != (JJLPptr)0);
    if (lp == (JJLPptr)0) return 1;

    // Cell 0 -> NLP column 1 (in the network), cell 1 -> not in it.
    const int mapping[2] = { 1, -1 };
    CHECK(JJtest_netFreeColState(lp, 2, mapping) == 0);

    double lo[3], hi[3];
    read_bounds(lp, 0, 2, lo, hi);
    CHECK(lo[0] == 0.0 && hi[0] == 10.0);
    CHECK(lo[1] == 0.0 && hi[1] == 10.0);
    CHECK(lo[2] == 0.0 && hi[2] == 10.0);

    // In-network cell: free, then restore. This is the normal path
    // and must keep working.
    double bd[2];
    free_col(0, bd);
    CHECK(bd[0] == 0.0 && bd[1] == 10.0);      // saved bounds
    read_bounds(lp, 0, 2, lo, hi);
    CHECK(lo[1] <= -1e19 && hi[1] >= 1e19);    // effectively unbounded
    unfree_col(0, bd);
    read_bounds(lp, 0, 2, lo, hi);
    CHECK(lo[1] == 0.0 && hi[1] == 10.0);      // restored

    // Out-of-network cell: must be a clean no-op.
    double bd2[2] = { -7.0, -8.0 };
    free_col(1, bd2);
    unfree_col(1, bd2);
    read_bounds(lp, 0, 2, lo, hi);
    CHECK(lo[0] == 0.0 && hi[0] == 10.0);
    CHECK(lo[1] == 0.0 && hi[1] == 10.0);
    CHECK(lo[2] == 0.0 && hi[2] == 10.0);
    CHECK(bd2[0] == -7.0 && bd2[1] == -8.0);   // buffer untouched

    JJfreeprob(&lp);

    if (failures) {
        std::printf("net_freecol: %d failure(s)\n", failures);
        return 1;
    }
    std::printf("net_freecol: all checks passed\n");
    return 0;
}
