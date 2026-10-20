/*******************************************************************/
/*  Interface to the LP solver: HiGHS (VSCIP)                     */
/*                                                                 */
/*  Rewritten on top of the HiGHS C API to replace the old         */
/*  CPLEX/XPRESS/SCIPlpi back-ends.                               */
/*                                                                 */
/*  The LP is handed to HiGHS in its original form: the objective */
/*  keeps its sense (max stays max) and rows keep their L/G/E     */
/*  bounds. HiGHS therefore returns row duals and reduced costs    */
/*  in the same convention as CPLEX 5, which the rest of the CSP  */
/*  code expects.                                                  */
/*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <unordered_map>
#include <vector>
#include "Jjsolver.h"

/* Shadow objective for skipping no-op cost changes.
 * The LP pointer is the key; what makes it a valid identity is the
 * lifecycle pairing below (every HiGHS object is created in
 * JJloadprob and destroyed in JJfreeprob, nothing else):
 *   JJloadprob                     install from the passed objective
 *   JJaddcols, JJaddrows(ccnt>0)   extend with the new column costs
 *   JJdelcols, JJdelsetcols        compact with the delete mask
 *   JJfreeprob                     erase before Highs_destroy
 * A cost change equal to the shadowed value (exact double equality)
 * is not forwarded to HiGHS: the model is left untouched, so any
 * solver state HiGHS has cached stays valid. The shadow is updated
 * only after a successful HiGHS call. */
struct JjCostShadow {
    std::vector<double> cost;
};
static std::unordered_map<const void*, JjCostShadow> g_costshadow;

/* Compacts the shadow after Highs_deleteColsByMask. HiGHS reuses the
 * mask array in place: on return, mask[i] is the new index of kept
 * column i (>= 0) or -1 for a deleted column. Kept columns keep their
 * relative order, so filtering the old order by mask[i] >= 0 yields
 * the shadow in the new column order. */
static void jj_shadow_delcols(const void* lp, HighsInt mac, const HighsInt* mask)
{
    auto sit = g_costshadow.find(lp);
    if (sit == g_costshadow.end()) return;
    if (sit->second.cost.size() != (size_t)mac) {
        /* drift: stop skipping for this LP rather than risk a
         * wrong skip on a stale index mapping */
        g_costshadow.erase(sit);
        return;
    }
    std::vector<double> keep;
    keep.reserve((size_t)mac);
    for (HighsInt i = 0; i < mac; i++)
        if (mask[i] >= 0) keep.push_back(sit->second.cost[i]);
    sit->second.cost.swap(keep);
}

#ifdef TAUARGUS_CSP_UNIT_TEST
/* Test-only: desynchronize the shadow from the model so the next column
 * delete exercises the drift fallback (the shadow must be dropped, and
 * changes must then be forwarded, never mis-skipped). */
void JJtest_corruptCostShadow(JJLPptr lp)
{
    auto sit = g_costshadow.find(lp);
    if (sit != g_costshadow.end())
        sit->second.cost.resize(sit->second.cost.size() + 1);
}

/* Test-only: number of live shadows (verify install/free pairing). */
long JJtest_shadowCount(void)
{
    return (long)g_costshadow.size();
}
#endif

extern char* fproblemlp;
extern double FEAS_TOL;
extern double OPT_TOL;

#ifdef VSCIP
#define JJ_HIGHS_INF  1.0e20
static unsigned int lpJJ = 0;
#endif

#ifdef VSCIP
/* A bound is treated as infinite once its magnitude reaches the
 * HiGHS infinity threshold (1e20). The CSP uses 1e20 for its own
 * infinite bounds, so this matches both sides. */
static int JJ_isinf(double v) { return v >=  JJ_HIGHS_INF || v <= -JJ_HIGHS_INF; }

/* Map a HiGHS model status to the legacy JJ_* status codes. */
static int JJ_status(void* lp)
{
    int ms = (int)Highs_getModelStatus(lp);
    switch (ms) {
        case kHighsModelStatusOptimal:                    return JJ_OPTIMAL;
        case kHighsModelStatusUnbounded:                  return JJ_UNBOUNDED;
        case kHighsModelStatusInfeasible:                 return JJ_INFEASIBLE;
        case kHighsModelStatusTimeLimit:                  return 11;
        case kHighsModelStatusIterationLimit:             return 10;
        case kHighsModelStatusObjectiveBound:
        case kHighsModelStatusObjectiveTarget:            return 12;
        case kHighsModelStatusUnboundedOrInfeasible:      return JJ_INForUNB;
        default:                                          return JJ_INForUNB;
    }
}

/* Solve the LP with the dual simplex method. Returns 0 on success. */
static int JJ_dual_solve(void* lp)
{
    Highs_setStringOptionValue(lp, "dual_feasibility_tolerance", "1e-7");
    HighsInt st = Highs_run(lp);
    if (st == kHighsStatusError) return 1;
    int ms = (int)Highs_getModelStatus(lp);
    if (ms == kHighsModelStatusUnbounded || ms == kHighsModelStatusInfeasible)
        return 1;
    return 0;
}

/* Solve the LP with the primal simplex method. Returns 0 on success. */
static int JJ_primal_solve(void* lp)
{
    HighsInt st = Highs_run(lp);
    if (st == kHighsStatusError) return 1;
    int ms = (int)Highs_getModelStatus(lp);
    if (ms == kHighsModelStatusUnbounded || ms == kHighsModelStatusInfeasible)
        return 1;
    return 0;
}
#endif

#ifdef VSCIP

/* ------------------------------- load / free ------------------------------- */

JJLPptr JJloadprob (
char *probname,
int numcols, int numrows, int numrims, int objsen,
double *obj, double *rhs,
char *sense,
int *matbeg, int *matcnt, int *matind,
double *matval, double *lb, double *ub, double *rngval,
int *freerowind, int *rimtype, int *rimbeg, int *rimcnt, int *rimind,
double *rimval,
char *dataname, char *objname, char *rhsname, char *rngname, char *bndname,
char **colname,
char *colnamestore,
char **rowname,
char *rownamestore,
char **rimname,
char *rimnamestore,
int colspace, int rowspace, int nzspace,
int rimspace, int rimnzspace,
unsigned colnamespace, unsigned rownamespace, unsigned rimnamespace
)
{
    (void)probname; (void)numrims; (void)rngval;
    (void)freerowind; (void)rimtype; (void)rimbeg; (void)rimcnt; (void)rimind;
    (void)rimval; (void)dataname; (void)objname; (void)rhsname; (void)rngname;
    (void)bndname; (void)colname; (void)colnamestore; (void)rowname;
    (void)rownamestore; (void)rimname; (void)rimnamestore; (void)colspace;
    (void)rowspace; (void)nzspace; (void)rimspace; (void)rimnzspace;
    (void)colnamespace; (void)rownamespace; (void)rimnamespace;

    void* H = Highs_create();
    if (!H) return NULL;
    Highs_setBoolOptionValue(H, "output_flag", 0);

    HighsInt mac = numcols, mar = numrows;
    /* The CSP stores the matrix in column-wise sparse form where matbeg[k]
     * is the offset of the first element of column k. HiGHS wants a_start
     * of length mac+1 with a_start[mac] = total nz. */
    HighsInt num_nz = (mac > 0) ? (matbeg[mac - 1] + matcnt[mac - 1]) : 0;
    HighsInt* a_start = (HighsInt*)malloc((mac + 1) * sizeof(HighsInt));
    for (int k = 0; k < mac; k++) a_start[k] = matbeg[k];
    a_start[mac] = num_nz;

    double* row_lower = (double*)malloc(mar * sizeof(double));
    double* row_upper = (double*)malloc(mar * sizeof(double));
    for (int i = 0; i < mar; i++) {
        if (sense[i] == 'L') { row_lower[i] = -JJ_HIGHS_INF; row_upper[i] = rhs[i]; }
        else if (sense[i] == 'G') { row_lower[i] = rhs[i]; row_upper[i] = JJ_HIGHS_INF; }
        else { row_lower[i] = rhs[i]; row_upper[i] = rhs[i]; }
    }

    HighsInt st = Highs_passLp(H, mac, mar, num_nz,
        kHighsMatrixFormatColwise, objsen, 0.0,
        obj, lb, ub, row_lower, row_upper,
        a_start, matind, matval);

    free(a_start);
    free(row_lower);
    free(row_upper);

    if (st != kHighsStatusOk) {
        Highs_destroy(H);
        return NULL;
    }

    {
        JjCostShadow& sh = g_costshadow[H];
        sh.cost.assign(obj, obj + mac);
    }

    lpJJ++;
    return H;
}

void JJfreeprob (JJLPptr *lp)
{
    if (lp == NULL || *lp == NULL) return;
    g_costshadow.erase(*lp);
    Highs_destroy(*lp);
    *lp = NULL;
    if (lpJJ > 0) lpJJ--;
}

int JJcheckprob (
char *probname, int numcols, int numrows, int numrims, int objsen,
double *obj, double *rhs, char *sense,
int *matbeg, int *matcnt, int *matind, double *matval,
double *lb, double *ub, double *rngval,
int *freerowind, int *rimtype, int *rimbeg, int *rimcnt, int *rimind,
double *rimval, char *dataname, char *objname, char *rhsname,
char *rngname, char *bndname,
char **colname, char *colnamestore, char **rowname, char *rownamestore,
char **rimname, char *rimnamestore,
int colspace, int rowspace, int nzspace,
int rimspace, int rimnzspace,
unsigned colnamespace, unsigned rownamespace, unsigned rimnamespace
)
{
    (void)probname;(void)numcols;(void)numrows;(void)numrims;(void)objsen;
    (void)obj;(void)rhs;(void)sense;(void)matbeg;(void)matcnt;(void)matind;
    (void)matval;(void)lb;(void)ub;(void)rngval;(void)freerowind;(void)rimtype;
    (void)rimbeg;(void)rimcnt;(void)rimind;(void)rimval;(void)dataname;
    (void)objname;(void)rhsname;(void)rngname;(void)bndname;(void)colname;
    (void)colnamestore;(void)rowname;(void)rownamestore;(void)rimname;
    (void)rimnamestore;(void)colspace;(void)rowspace;(void)nzspace;
    (void)rimspace;(void)rimnzspace;(void)colnamespace;(void)rownamespace;
    (void)rimnamespace;
    return 0;
}

/* ------------------------------- solve ---------------------------------- */

int JJoptimize (JJLPptr lp)
{
    return JJ_primal_solve(lp);
}

int JJdualopt (JJLPptr lp)
{
    return JJ_dual_solve(lp);
}

int JJnetopt (JJLPptr lp, int *netstatus_p, int *numnodes_p, int *numarcs_p, int *itcnt_p)
{
    *numnodes_p = 0;
    *numarcs_p  = 0;
    *itcnt_p    = 0;
    int rc = JJ_dual_solve(lp);
    *netstatus_p = (rc == 0) ? 0 : 1101;
    return *netstatus_p;
}

int JJmipopt (JJLPptr lp)
{
    return JJ_primal_solve(lp);
}

int JJmipoptimize (JJLPptr lp)
{
    JJoptimize(lp);
    return 0;
}

/* ------------------------------- dimension / info ----------------------- */

int JJgetmac (JJLPptr lp)  { return (int)Highs_getNumCol(lp); }
int JJgetmar (JJLPptr lp)  { return (int)Highs_getNumRow(lp); }
int JJgetmat (JJLPptr lp)  { return (int)Highs_getNumNz(lp); }
int JJgetobjsen (JJLPptr lp)
{
    HighsInt sense = 1;
    Highs_getObjectiveSense(lp, &sense);
    return (int)sense;
}

int JJgetitc (JJLPptr lp)
{
    HighsInt it = Highs_getSimplexIterationCount(lp);
    return (it >= 0) ? (int)it : 0;
}

int JJgetitci (JJLPptr lp)
{
    int status = JJgetstat(lp);
    if (status != JJ_OPTIMAL && status != JJ_UNBOUNDED && status != JJ_INFEASIBLE)
        return 0;
    return JJgetitc(lp);
}

int JJlpiterlimit (int val)
{
    (void)val;
    return 0;
}

int JJsetscr_ind (JJLPptr lp, int scr_ind)
{
    Highs_setBoolOptionValue(lp, "output_flag", scr_ind ? 1 : 0);
    return 0;
}

int JJlpwrite (JJLPptr lp, std::string f)
{
    Highs_writeModel(lp, f.c_str());
    return 0;
}

int JJmpswrite (JJLPptr lp, char *f)
{
    Highs_writeModel(lp, f);
    return 0;
}

/* ------------------------------- data getters --------------------------- */

int JJgetobj (JJLPptr lp, double *obj, int begin, int end)
{
    HighsInt mac = Highs_getNumCol(lp);
    double* cost = (double*)malloc(mac * sizeof(double));
    HighsInt num_col = mac, num_row = 0, num_nz = 0, sense = 1;
    double offset = 0.0;
    double* cl = (double*)malloc(mac * sizeof(double));
    double* cu = (double*)malloc(mac * sizeof(double));
    Highs_getLp(lp, kHighsMatrixFormatColwise, &num_col, &num_row, &num_nz,
                &sense, &offset, cost, cl, cu, NULL, NULL, NULL, NULL, NULL, NULL);
    for (int i = begin; i <= end; i++)
        obj[i - begin] = (sense == kHighsObjSenseMaximize) ? -cost[i] : cost[i];
    free(cost); free(cl); free(cu);
    return 0;
}

int JJgetrhs (JJLPptr lp, double *rhs, int begin, int end)
{
    int n = end - begin + 1;
    double* lo = (double*)malloc(n * sizeof(double));
    double* hi = (double*)malloc(n * sizeof(double));
    HighsInt num_row = n, num_nz = 0;
    Highs_getRowsByRange(lp, begin, end, &num_row, lo, hi, &num_nz, NULL, NULL, NULL);
    for (int i = 0; i < n; i++)
        rhs[i] = JJ_isinf(hi[i]) ? lo[i] : hi[i];
    free(lo); free(hi);
    return 0;
}

int JJgetsense (JJLPptr lp, char *sense, int begin, int end)
{
    int n = end - begin + 1;
    double* lo = (double*)malloc(n * sizeof(double));
    double* hi = (double*)malloc(n * sizeof(double));
    HighsInt num_row = n, num_nz = 0;
    Highs_getRowsByRange(lp, begin, end, &num_row, lo, hi, &num_nz, NULL, NULL, NULL);
    for (int i = 0; i < n; i++)
        if (JJ_isinf(lo[i])) sense[i] = 'L';
        else if (JJ_isinf(hi[i])) sense[i] = 'G';
        else sense[i] = 'E';
    free(lo); free(hi);
    return 0;
}

int JJgetcols (JJLPptr lp, int *nzcnt, int *cmatbeg, int *cmatind,
              double *cmatval, int cmatspace, int *surplus, int begin, int end)
{
    (void)cmatspace;
    HighsInt num_col, num_nz;
    Highs_getColsByRange(lp, begin, end, &num_col, NULL, NULL, NULL, &num_nz,
                         cmatbeg, cmatind, cmatval);
    *nzcnt = (int)num_nz;
    *surplus = cmatspace - *nzcnt;
    return 0;
}

int JJgetrows (JJLPptr lp, int *nzcnt, int *rmatbeg, int *rmatind, double *rmatval,
              int rmatspace, int *surplus, int begin, int end)
{
    (void)rmatspace;
    HighsInt num_row, num_nz;
    Highs_getRowsByRange(lp, begin, end, &num_row, NULL, NULL, &num_nz,
                         rmatbeg, rmatind, rmatval);
    *nzcnt = (int)num_nz;
    *surplus = rmatspace - *nzcnt;
    return 0;
}

int JJgetbdl (JJLPptr lp, double *bd, int begin, int end)
{
    int n = end - begin + 1;
    double* lo = (double*)malloc(n * sizeof(double));
    double* hi = (double*)malloc(n * sizeof(double));
    HighsInt num_col = n, num_nz = 0;
    Highs_getColsByRange(lp, begin, end, &num_col, NULL, lo, hi, &num_nz, NULL, NULL, NULL);
    for (int i = 0; i < n; i++) bd[i] = lo[i];
    free(lo); free(hi);
    return 0;
}

int JJgetbdu (JJLPptr lp, double *bd, int begin, int end)
{
    int n = end - begin + 1;
    double* lo = (double*)malloc(n * sizeof(double));
    double* hi = (double*)malloc(n * sizeof(double));
    HighsInt num_col = n, num_nz = 0;
    Highs_getColsByRange(lp, begin, end, &num_col, NULL, lo, hi, &num_nz, NULL, NULL, NULL);
    for (int i = 0; i < n; i++) bd[i] = hi[i];
    free(lo); free(hi);
    return 0;
}

/* ------------------------------- solution ------------------------------- */

int JJgetstat (JJLPptr lp)
{
    return JJ_status(lp);
}

int JJgetobjval (JJLPptr lp, double *objval)
{
    *objval = Highs_getObjectiveValue(lp);
    return 0;
}

int JJgetx (JJLPptr lp, double *x, int begin, int end)
{
    HighsInt mac = Highs_getNumCol(lp);
    double* colv = (double*)malloc(mac * sizeof(double));
    Highs_getSolution(lp, colv, NULL, NULL, NULL);
    for (int i = begin; i <= end; i++) x[i - begin] = colv[i];
    free(colv);
    return 0;
}

int JJgetpi (JJLPptr lp, double *pi, int begin, int end)
{
    HighsInt mar = Highs_getNumRow(lp);
    double* rowd = (double*)malloc(mar * sizeof(double));
    Highs_getSolution(lp, NULL, NULL, NULL, rowd);
    for (int i = begin; i <= end; i++) pi[i - begin] = rowd[i];
    free(rowd);
    return 0;
}

int JJgetslack (JJLPptr lp, double *slack, int begin, int end)
{
    HighsInt mar = Highs_getNumRow(lp);
    double* rowv = (double*)malloc(mar * sizeof(double));
    double* lo   = (double*)malloc(mar * sizeof(double));
    double* hi   = (double*)malloc(mar * sizeof(double));
    Highs_getSolution(lp, NULL, NULL, rowv, NULL);
    HighsInt num_row = mar, num_nz = 0;
    Highs_getRowsByRange(lp, 0, mar - 1, &num_row, lo, hi, &num_nz, NULL, NULL, NULL);
    for (int i = begin; i <= end; i++)
        slack[i - begin] = JJ_isinf(hi[i]) ? lo[i] - rowv[i] : hi[i] - rowv[i];
    free(rowv); free(lo); free(hi);
    return 0;
}

int JJgetdj (JJLPptr lp, double *dj, int begin, int end)
{
    HighsInt mac = Highs_getNumCol(lp);
    double* cold = (double*)malloc(mac * sizeof(double));
    Highs_getSolution(lp, NULL, cold, NULL, NULL);
    for (int i = begin; i <= end; i++) dj[i - begin] = cold[i];
    free(cold);
    return 0;
}

int JJsolution (JJLPptr lp, int *lpstat_p, double *objval_p, double *x,
                double *pi, double *slack, double *dj)
{
    if (lpstat_p)  *lpstat_p  = JJgetstat(lp);
    if (objval_p)  *objval_p  = Highs_getObjectiveValue(lp);
    if (x)     JJgetx(lp, x, 0, JJgetmac(lp) - 1);
    if (pi)    JJgetpi(lp, pi, 0, JJgetmar(lp) - 1);
    if (slack) JJgetslack(lp, slack, 0, JJgetmar(lp) - 1);
    if (dj)    JJgetdj(lp, dj, 0, JJgetmac(lp) - 1);
    return 0;
}

int JJgetmx (JJLPptr lp, double *x, int begin, int end)
{
    return JJgetx(lp, x, begin, end);
}

/* ------------------------------- basis ---------------------------------- */

int JJgetbase (JJLPptr lp, int *cstat, int *rstat)
{
    HighsInt mac = Highs_getNumCol(lp);
    HighsInt mar = Highs_getNumRow(lp);
    int* tc = cstat ? cstat : (int*)malloc(mac * sizeof(int));
    int* tr = rstat ? rstat : (int*)malloc(mar * sizeof(int));
    Highs_getBasis(lp, tc, tr);
    if (!cstat) free(tc);
    if (!rstat) free(tr);
    return 0;
}

int JJloadbase (JJLPptr lp, int *cstat, int *rstat)
{
    if (!cstat || !rstat) return 0;
    Highs_setBasis(lp, cstat, rstat);
    return 0;
}

/* b = B^{-1} (effective rhs). The effective rhs of a row is its lower
 * bound when the upper is infinite ('G'/'>='), its upper bound when the
 * lower is infinite ('L'/'<='), or either for equality rows. */
static void JJ_effective_rhs(void* H, double* b)
{
    HighsInt mar = Highs_getNumRow(H);
    double* lo = (double*)malloc(mar * sizeof(double));
    double* hi = (double*)malloc(mar * sizeof(double));
    HighsInt num_row = mar, num_nz = 0;
    Highs_getRowsByRange(H, 0, mar - 1, &num_row, lo, hi, &num_nz, NULL, NULL, NULL);
    for (int i = 0; i < mar; i++)
        b[i] = JJ_isinf(hi[i]) ? lo[i] : hi[i];
    free(lo); free(hi);
}

/* Values of the current basic variables: x_B = B^{-1} b. */
int JJgetbhead (JJLPptr lp, int *head, double *x)
{
    HighsInt mar = Highs_getNumRow(lp);
    HighsInt mac = Highs_getNumCol(lp);

    HighsInt* basic = (HighsInt*)malloc(mar * sizeof(HighsInt));
    Highs_getBasicVariables(lp, basic);

    double* colv = (double*)malloc(mac * sizeof(double));
    double* rowv = (double*)malloc(mar * sizeof(double));
    Highs_getSolution(lp, colv, NULL, rowv, NULL);

    for (int i = 0; i < mar; i++) {
        HighsInt v = basic[i];
        if (v >= 0 && v < mac) { head[i] = (int)v;       x[i] = colv[v]; }
        else if (v < -1)       { head[i] = -(int)(v + 1); x[i] = -rowv[-(int)(v + 1)]; }
        else                   { head[i] = (int)mac;     x[i] = 0.0; }
    }
    free(basic); free(colv); free(rowv);
    return 0;
}

/* Row i of the inverse basis matrix B^{-1} (dense, length mar). */
int JJbinvrow (JJLPptr lp, int i, double *y)
{
    HighsInt mar = Highs_getNumRow(lp);
    double* vec = (double*)calloc(mar, sizeof(double));
    HighsInt nnz = 0;
    HighsInt* idx = (HighsInt*)malloc(mar * sizeof(HighsInt));
    Highs_getBasisInverseRow(lp, i, vec, &nnz, idx);
    memcpy(y, vec, mar * sizeof(double));
    free(vec); free(idx);
    return 0;
}

/* Row i of B^{-1}A (dense, length mac). */
int JJbinvarow (JJLPptr lp, int i, double *z)
{
    HighsInt mac = Highs_getNumCol(lp);
    double* vec = (double*)calloc(mac, sizeof(double));
    HighsInt nnz = 0;
    HighsInt* idx = (HighsInt*)malloc(mac * sizeof(HighsInt));
    Highs_getReducedRow(lp, i, vec, &nnz, idx);
    memcpy(z, vec, mac * sizeof(double));
    free(vec); free(idx);
    return 0;
}

int JJbinvcol (JJLPptr lp, int j, double *y)
{
    HighsInt mar = Highs_getNumRow(lp);
    double* vec = (double*)calloc(mar, sizeof(double));
    HighsInt nnz = 0;
    HighsInt* idx = (HighsInt*)malloc(mar * sizeof(HighsInt));
    Highs_getBasisInverseCol(lp, j, vec, &nnz, idx);
    memcpy(y, vec, mar * sizeof(double));
    free(vec); free(idx);
    return 0;
}

/* ------------------------------- model edits ---------------------------- */

/* The HiGHS LP back-end solves continuous relaxations; variable integrality
 * types are not tracked. This is a no-op stub kept for API compatibility. */
int JJcopyctype (JJLPptr lp, char *xctype)
{
    (void)lp; (void)xctype;
    return 1;
}

int JJaddrows (JJLPptr lp, int ccnt, int rcnt, int nzcnt, double *rhs, char *sense,
              int *rmatbeg, int *rmatind, double *rmatval, char **colname, char **rowname)
{
    (void)colname; (void)rowname;

    if (ccnt > 0) {
        HighsInt oldmac = Highs_getNumCol(lp);
        HighsInt* starts = (HighsInt*)malloc(ccnt * sizeof(HighsInt));
        double* costs = (double*)malloc(ccnt * sizeof(double));
        double* lo = (double*)malloc(ccnt * sizeof(double));
        double* hi = (double*)malloc(ccnt * sizeof(double));
        for (int j = 0; j < ccnt; j++) {
            starts[j] = 0;
            costs[j]  = 0.0;
            lo[j]     = 0.0;
            hi[j]     = JJ_HIGHS_INF;
        }
        HighsInt stc = Highs_addCols(lp, ccnt, costs, lo, hi, 0, starts, NULL, NULL);
        if (stc == kHighsStatusOk) {
            auto sit = g_costshadow.find(lp);
            if (sit != g_costshadow.end())
                sit->second.cost.resize(sit->second.cost.size() + ccnt, 0.0);
        }
        free(starts); free(costs); free(lo); free(hi);
        (void)oldmac;
    }

    double* lo = (double*)malloc(rcnt * sizeof(double));
    double* hi = (double*)malloc(rcnt * sizeof(double));
    for (int i = 0; i < rcnt; i++) {
        if (sense[i] == 'L') { lo[i] = -JJ_HIGHS_INF; hi[i] = rhs[i]; }
        else if (sense[i] == 'G') { lo[i] = rhs[i]; hi[i] = JJ_HIGHS_INF; }
        else { lo[i] = rhs[i]; hi[i] = rhs[i]; }
    }
    HighsInt st = Highs_addRows(lp, rcnt, lo, hi, nzcnt, rmatbeg, rmatind, rmatval);
    free(lo); free(hi);
    return (st != kHighsStatusOk) ? 1 : 0;
}

int JJaddcols (JJLPptr lp, int ccnt, int nzcnt, double *cobj, int *cmatbeg,
              int *cmatind, double *cmatval, double *cbdl, double *cbdu, char **colname)
{
    (void)colname;
    HighsInt st = Highs_addCols(lp, ccnt, cobj, cbdl, cbdu, nzcnt, cmatbeg, cmatind, cmatval);
    if (st == kHighsStatusOk) {
        auto sit = g_costshadow.find(lp);
        if (sit != g_costshadow.end())
            sit->second.cost.insert(sit->second.cost.end(), cobj, cobj + ccnt);
    }
    return (st != kHighsStatusOk) ? 1 : 0;
}

int JJdelrows (JJLPptr lp, int begin, int end)
{
    HighsInt mar = Highs_getNumRow(lp);
    HighsInt* mask = (HighsInt*)calloc(mar, sizeof(HighsInt));
    for (HighsInt i = begin; i <= end; i++) mask[i] = 1;
    HighsInt st = Highs_deleteRowsByMask(lp, mask);
    free(mask);
    return (st != kHighsStatusOk) ? 1 : 0;
}

int JJdelsetrows (JJLPptr lp, int *delstat)
{
    HighsInt mar = Highs_getNumRow(lp);
    HighsInt* mask = (HighsInt*)calloc(mar, sizeof(HighsInt));
    for (HighsInt i = 0; i < mar; i++) mask[i] = delstat[i] ? 1 : 0;
    HighsInt st = Highs_deleteRowsByMask(lp, mask);
    free(mask);
    return (st != kHighsStatusOk) ? 1 : 0;
}

int JJdelcols (JJLPptr lp, int begin, int end)
{
    HighsInt mac = Highs_getNumCol(lp);
    HighsInt* mask = (HighsInt*)calloc(mac, sizeof(HighsInt));
    for (HighsInt i = begin; i <= end; i++) mask[i] = 1;
    HighsInt st = Highs_deleteColsByMask(lp, mask);
    if (st == kHighsStatusOk) jj_shadow_delcols(lp, mac, mask);
    free(mask);
    return (st != kHighsStatusOk) ? 1 : 0;
}

int JJdelsetcols (JJLPptr lp, int *delstat)
{
    HighsInt mac = Highs_getNumCol(lp);
    HighsInt* mask = (HighsInt*)calloc(mac, sizeof(HighsInt));
    for (HighsInt i = 0; i < mac; i++) mask[i] = delstat[i] ? 1 : 0;
    HighsInt st = Highs_deleteColsByMask(lp, mask);
    if (st == kHighsStatusOk) jj_shadow_delcols(lp, mac, mask);
    free(mask);
    return (st != kHighsStatusOk) ? 1 : 0;
}

int JJchgcoef (JJLPptr lp, int i, int j, double newvalue)
{
    if (i == -1) {
        /* change objective cost of column j; skip when the value is
         * already in effect (exact equality against the shadow) */
        auto sit = g_costshadow.find(lp);
        if (sit != g_costshadow.end() && j >= 0
            && (size_t)j < sit->second.cost.size()
            && sit->second.cost[j] == newvalue)
            return 0;
        HighsInt st = Highs_changeColCost(lp, j, newvalue);
        if (st == kHighsStatusOk && sit != g_costshadow.end() && j >= 0
            && (size_t)j < sit->second.cost.size())
            sit->second.cost[j] = newvalue;
        return 0;
    }
    if (j == -1) {
        /* change the effective rhs of row i */
        double* lo = (double*)malloc(sizeof(double));
        double* hi = (double*)malloc(sizeof(double));
        HighsInt num_row = 1, num_nz = 0;
        Highs_getRowsByRange(lp, i, i, &num_row, lo, hi, &num_nz, NULL, NULL, NULL);
        /* skip when the bound is already at the requested value
         * (exact equality; the lookup above is needed for the
         * opposite bound anyway) */
        if (JJ_isinf(lo[0])) {
            if (hi[0] != newvalue) Highs_changeRowBounds(lp, i, lo[0], newvalue);
        } else if (JJ_isinf(hi[0])) {
            if (lo[0] != newvalue) Highs_changeRowBounds(lp, i, newvalue, hi[0]);
        } else {
            if (lo[0] != newvalue || hi[0] != newvalue)
                Highs_changeRowBounds(lp, i, newvalue, newvalue);
        }
        free(lo); free(hi);
        return 0;
    }
    Highs_changeCoeff(lp, i, j, newvalue);
    return 0;
}

int JJchgbds (JJLPptr lp, int cnt, int *index, char *lu, double *bd)
{
    for (int i = 0; i < cnt; i++) {
        int k = index[i];
        double* lo = (double*)malloc(sizeof(double));
        double* hi = (double*)malloc(sizeof(double));
        HighsInt num_col = 1, num_nz = 0;
        Highs_getColsByRange(lp, k, k, &num_col, NULL, lo, hi, &num_nz, NULL, NULL, NULL);
        /* skip when the bound is already at the requested value
         * (exact equality; the lookup above is needed for the
         * opposite bound anyway) */
        if (lu[i] == 'L') {
            if (lo[0] != bd[i]) Highs_changeColBounds(lp, k, bd[i], hi[0]);
        } else if (lu[i] == 'U') {
            if (hi[0] != bd[i]) Highs_changeColBounds(lp, k, lo[0], bd[i]);
        } else {
            if (lo[0] != bd[i] || hi[0] != bd[i])
                Highs_changeColBounds(lp, k, bd[i], bd[i]);
        }
        free(lo); free(hi);
    }
    return 0;
}

int JJchgobj (JJLPptr lp, int cnt, int *index, double *obj)
{
    auto sit = g_costshadow.find(lp);
    for (int i = 0; i < cnt; i++) {
        if (sit != g_costshadow.end() && index[i] >= 0
            && (size_t)index[i] < sit->second.cost.size()
            && sit->second.cost[index[i]] == obj[i])
            continue;
        HighsInt st = Highs_changeColCost(lp, index[i], obj[i]);
        if (st == kHighsStatusOk && sit != g_costshadow.end() && index[i] >= 0
            && (size_t)index[i] < sit->second.cost.size())
            sit->second.cost[index[i]] = obj[i];
    }
    return 0;
}

int JJchgrhs (JJLPptr lp, int cnt, int *index, double *rhs)
{
    for (int i = 0; i < cnt; i++)
        JJchgcoef(lp, index[i], -1, rhs[i]);
    return 0;
}

#endif /* VSCIP */
