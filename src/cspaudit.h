/*************************************************************************/
/* Audit (realized feasibility intervals) for cell-suppressed tables    */
/* in Statistical Disclosure Control (SDC)                              */
/*                                                                      */
/* CSPAUDIT.h                                                           */
/* Native port of the standalone `intervalle.exe` audit program         */
/* (reference/intervalle, Pascal/FPC) into the tauargus_csp module.     */
/*                                                                      */
/* Contract (legacy intervalle defaults: dateirechnen=primsec,          */
/* p0=0, p1=50, p2=50, p3=0):                                           */
/*   - Variables: every JJ cell with status 'u' (primary) or 'm'        */
/*     (secondary), in ascending cell-index order.                      */
/*   - Bounds: [Int(LBound), Int(UBound)] per cell, taken straight      */
/*     from the JJ cell line (p0=0 branch of the legacy bounds code).   */
/*   - Constraints: the JJ restrictions that contain at least one       */
/*     suppressed cell; safe-cell contributions are moved to the        */
/*     right-hand side.  All rows are equalities.                       */
/*   - Per cell: minimise the cell's variable -> realized minimum,      */
/*     maximise -> realized maximum.                                    */
/*   - Unsafe flag: 1 if (max - value) < UPL or (max - min) < LPL,      */
/*     mirroring the legacy protection check.                           */
/*                                                                      */
/* The JJ input is exactly what CJJFormat::WriteJJFormat writes with    */
/* AsPerc=false (the same call the legacy SaveTable.writeJJ made).      */
/*************************************************************************/
#ifndef CSPAUDIT_H
#define CSPAUDIT_H

#ifdef __cplusplus
extern "C" {
#endif

/* Audit the suppressed cells of the table described by the JJ file.
 *
 * On success (*n_out > 0 or 0) six malloc'd arrays are returned, one
 * entry per audited cell, in ascending cell-index order:
 *   cell_idx  0-based JJ cell index
 *   a_min     realized minimum of the cell value
 *   a_max     realized maximum of the cell value
 *   a_value   response value of the cell (from the JJ line)
 *   a_status  suppression status, 'u' or 'm'
 *   a_unsafe  1 if the legacy protection check fails, else 0
 *
 * Return codes:
 *    0  success (also when the table has no suppressed cells: n_out=0)
 *   -1  JJ file could not be opened
 *   -2  JJ file parse error
 *   -3  HiGHS model setup failed
 *   -4  LP infeasible / unbounded for some cell (legacy aborts the
 *       whole audit on the first non-optimal solve; *n_out is set to 0)
 *
 * On rc != 0 no arrays are returned and TauAuditJjFree is not needed.
 * On rc == 0 with *n_out > 0 the caller must free via TauAuditJjFree.
 */
int TauAuditJj(const char *jjfile,
               int *n_out,
               long **cell_idx,
               double **a_min,
               double **a_max,
               double **a_value,
               char **a_status,
               int **a_unsafe);

/* Release the arrays returned by TauAuditJj. */
void TauAuditJjFree(int n, long *cell_idx, double *a_min, double *a_max,
                    double *a_value, char *a_status, int *a_unsafe);

#ifdef __cplusplus
}
#endif

#endif /* CSPAUDIT_H */
