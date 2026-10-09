/*************************************************************************/
/* Audit (realized feasibility intervals) for cell-suppressed tables    */
/* in Statistical Disclosure Control (SDC)                              */
/*                                                                      */
/*                                                                      */
/* CSPAUDIT.c  (HiGHS backend)                                          */
/* Native port of the standalone `intervalle.exe` audit program         */
/* (reference/intervalle, Pascal/FPC) into the tauargus_csp module.     */
/*                                                                      */
/* The audit subproblem is a pure LP that is loaded once and then       */
/* re-solved 2 x nvar times with a single column's objective cost       */
/* toggled between +1 / -1 (mirroring native/crp/src/crpSaudit.c).      */
/*                                                                      */
/* Legacy defaults (dateirechnen=primsec, p0=0, p1=50, p2=50, p3=0):    */
/*   variables = JJ cells with status 'u' or 'm'                        */
/*   bounds    = [Int(LBound), Int(UBound)] per cell from the JJ line   */
/*   rows      = JJ restrictions containing >= 1 suppressed cell,       */
/*               safe-cell terms moved to the RHS, sense = E            */
/*   objective = x_j for the cell under audit, minimise and maximise    */
/*                                                                      */
/* The legacy program aborts the whole audit (raises) on the first      */
/* non-optimal solve; we mirror that by returning -4.                   */
/*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <highs/interfaces/highs_c_api.h>

#include "cspaudit.h"

#define AUDIT_MAXLINE 4096

typedef struct {
    long   cell;    /* 0-based JJ cell index (== line order) */
    char   status;  /* 'u', 'm', 's', 'z'                    */
    double value;   /* response value                        */
    double lo, hi;  /* per-cell bounds from the JJ line      */
    double lpl, upl;/* protection levels                     */
} AuditCell;

typedef struct {
    long   cell;
    double coef;
} AuditTerm;

typedef struct {
    double rhs;
    int    nterms;
    AuditTerm *terms;
} AuditRestriction;

/* Split `line` (in place) on spaces/tabs; returns token count. */
static int tokenize(char *line, char **tok, int maxtok)
{
    int n = 0;
    char *p = line;
    while (n < maxtok) {
        while (*p == ' ' || *p == '\t') p++;
        if (!*p) break;
        tok[n++] = p;
        while (*p && *p != ' ' && *p != '\t') p++;
        if (*p) { *p = '\0'; p++; }
    }
    return n;
}

/* Parse the JJ file. Returns 0 on success, -1 open error, -2 parse error. */
static int read_jj(const char *jjfile, int *n_cell_out, AuditCell **cells_out,
                   int *n_rel_out, AuditRestriction **rels_out)
{
    FILE *f;
    char line[AUDIT_MAXLINE];
    int ncell, nrel, i;
    AuditCell *cells;
    AuditRestriction *rels;

    f = fopen(jjfile, "r");
    if (!f) {
        return -1;
    }

    /* line 0: fixed value ("0") — skip */
    if (!fgets(line, sizeof(line), f)) {
        fclose(f);
        return -2;
    }

    /* line 1: number of cells */
    if (!fgets(line, sizeof(line), f) ||
        (ncell = atoi(line)) <= 0 || ncell > 1000000) {
        fclose(f);
        return -2;
    }

    cells = (AuditCell *)calloc((size_t)ncell, sizeof(AuditCell));
    if (!cells) {
        fclose(f);
        return -2;
    }

    for (i = 0; i < ncell; i++) {
        char *tok[16];
        int nt;
        if (!fgets(line, sizeof(line), f)) {
            free(cells);
            fclose(f);
            return -2;
        }
        nt = tokenize(line, tok, 16);
        /* idx value cost status lo hi lpl upl sliding */
        if (nt < 9) {
            free(cells);
            fclose(f);
            return -2;
        }
        cells[i].cell = (long)i;
        cells[i].value = atof(tok[1]);
        cells[i].status = tok[3][0];
        cells[i].lo = atof(tok[4]);
        cells[i].hi = atof(tok[5]);
        cells[i].lpl = atof(tok[6]);
        cells[i].upl = atof(tok[7]);
    }

    /* restriction count */
    if (!fgets(line, sizeof(line), f) || (nrel = atoi(line)) < 0 || nrel > 1000000) {
        free(cells);
        fclose(f);
        return -2;
    }

    rels = (AuditRestriction *)calloc((size_t)(nrel > 0 ? nrel : 1),
                                      sizeof(AuditRestriction));
    if (!rels) {
        free(cells);
        fclose(f);
        return -2;
    }

    for (i = 0; i < nrel; i++) {
        char *tok[AUDIT_MAXLINE / 2];
        int nt, n;
        if (!fgets(line, sizeof(line), f)) {
            free(rels);
            free(cells);
            fclose(f);
            return -2;
        }
        nt = tokenize(line, tok, AUDIT_MAXLINE / 2);
        /* "0 n : c0 (k0) c1 (k1) ..." */
        if (nt < 3 + 2 * atoi(tok[1])) {
            free(rels);
            free(cells);
            fclose(f);
            return -2;
        }
        n = atoi(tok[1]);
        rels[i].rhs = atof(tok[0]);
        rels[i].nterms = n;
        rels[i].terms = (AuditTerm *)malloc((size_t)n * sizeof(AuditTerm));
        if (!rels[i].terms) {
            free(rels);
            free(cells);
            fclose(f);
            return -2;
        }
        for (int j = 0; j < n; j++) {
            /* coefficient token is "(1)" or "(-1)"; skip the '(' */
            rels[i].terms[j].cell = atol(tok[3 + 2 * j]);
            rels[i].terms[j].coef = (double)atoi(tok[3 + 2 * j + 1] + 1);
        }
    }

    fclose(f);
    *n_cell_out = ncell;
    *cells_out = cells;
    *n_rel_out = nrel;
    *rels_out = rels;
    return 0;
}

static int is_suppressed(char status)
{
    return status == 'u' || status == 'm';
}

int TauAuditJj(const char *jjfile,
               int *n_out,
               long **cell_idx,
               double **a_min,
               double **a_max,
               double **a_value,
               char **a_status,
               int **a_unsafe)
{
    int ncell, nrel, i, j, nvar = 0, nrow = 0;
    long *var_of = NULL;       /* cell -> var index (-1 = safe) */
    long *var_cell = NULL;     /* var index -> cell index */
    double *clb = NULL, *cub = NULL, *cost = NULL;
    int *active = NULL;        /* row -> active flag */
    int *active_row = NULL;    /* active row -> original restriction index */
    double *rlo = NULL, *rhi = NULL;
    HighsInt *mstart = NULL, *mindex = NULL, *colcount = NULL;
    double *mvalue = NULL;
    HighsInt nnz = 0;
    void *h = NULL;
    HighsInt st;
    int fail = 0;

    *n_out = 0;
    *cell_idx = NULL; *a_min = NULL; *a_max = NULL;
    *a_value = NULL; *a_status = NULL; *a_unsafe = NULL;

    AuditCell *cells = NULL;
    AuditRestriction *rels = NULL;
    if (read_jj(jjfile, &ncell, &cells, &nrel, &rels) != 0) {
        return -2;
    }

    /* variable map: every 'u'/'m' cell, ascending cell order */
    var_of = (long *)malloc((size_t)ncell * sizeof(long));
    var_cell = (long *)malloc((size_t)ncell * sizeof(long));
    clb = (double *)malloc((size_t)ncell * sizeof(double));
    cub = (double *)malloc((size_t)ncell * sizeof(double));
    cost = (double *)calloc((size_t)(ncell > 0 ? ncell : 1), sizeof(double));
    if (!var_of || !var_cell || !clb || !cub || !cost) {
        goto error;
    }
    for (i = 0; i < ncell; i++) {
        if (is_suppressed(cells[i].status)) {
            var_of[i] = nvar;
            var_cell[nvar] = i;
            /* legacy Int() truncates toward zero */
            clb[nvar] = (double)(long)cells[i].lo;
            cub[nvar] = (double)(long)cells[i].hi;
            nvar++;
        }
        else {
            var_of[i] = -1;
        }
    }

    /* active restrictions: those containing >= 1 suppressed cell */
    active = (int *)calloc((size_t)(nrel > 0 ? nrel : 1), sizeof(int));
    active_row = (int *)malloc((size_t)(nrel > 0 ? nrel : 1) * sizeof(int));
    rlo = (double *)malloc((size_t)(nrel > 0 ? nrel : 1) * sizeof(double));
    rhi = (double *)malloc((size_t)(nrel > 0 ? nrel : 1) * sizeof(double));
    if (!active || !active_row || !rlo || !rhi) {
        goto error;
    }
    for (i = 0; i < nrel; i++) {
        AuditRestriction *r = &rels[i];
        double rhs = r->rhs;
        int hasvar = 0;
        for (j = 0; j < r->nterms; j++) {
            if (var_of[r->terms[j].cell] >= 0) {
                hasvar = 1;
            }
            else {
                /* move the safe-cell contribution to the RHS */
                rhs -= r->terms[j].coef *
                       cells[r->terms[j].cell].value;
            }
        }
        if (hasvar) {
            active[i] = 1;
            active_row[nrow] = i;
            rlo[nrow] = rhs;
            rhi[nrow] = rhs; /* equality */
            nrow++;
        }
    }

    /* constraint matrix in column-major (CSC) form:
       column = variable, entries from the active rows */
    colcount = (HighsInt *)calloc((size_t)(nvar > 0 ? nvar : 1), sizeof(HighsInt));
    if (!colcount) {
        goto error;
    }
    for (i = 0; i < nrow; i++) {
        AuditRestriction *r = &rels[active_row[i]];
        for (j = 0; j < r->nterms; j++) {
            if (var_of[r->terms[j].cell] >= 0) {
                nnz++;
                colcount[var_of[r->terms[j].cell]]++;
            }
        }
    }
    mstart = (HighsInt *)malloc((size_t)(nvar + 1) * sizeof(HighsInt));
    mindex = (HighsInt *)malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(HighsInt));
    mvalue = (double *)malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(double));
    if (!mstart || !mindex || !mvalue) {
        goto error;
    }
    mstart[0] = 0;
    for (i = 1; i <= nvar; i++) {
        mstart[i] = mstart[i - 1] + colcount[i - 1];
    }
    {
        HighsInt *pos = (HighsInt *)malloc((size_t)(nvar > 0 ? nvar : 1) *
                                           sizeof(HighsInt));
        if (!pos) {
            goto error;
        }
        for (i = 0; i < nvar; i++) {
            pos[i] = mstart[i];
        }
        for (i = 0; i < nrow; i++) {
            AuditRestriction *r = &rels[active_row[i]];
            for (j = 0; j < r->nterms; j++) {
                long v = var_of[r->terms[j].cell];
                if (v >= 0) {
                    mindex[pos[v]] = (HighsInt)i;
                    mvalue[pos[v]] = r->terms[j].coef;
                    pos[v]++;
                }
            }
        }
        free(pos);
    }

    /* output arrays, one entry per audited cell */
    *cell_idx = (long *)calloc((size_t)(nvar > 0 ? nvar : 1), sizeof(long));
    *a_min = (double *)calloc((size_t)(nvar > 0 ? nvar : 1), sizeof(double));
    *a_max = (double *)calloc((size_t)(nvar > 0 ? nvar : 1), sizeof(double));
    *a_value = (double *)calloc((size_t)(nvar > 0 ? nvar : 1), sizeof(double));
    *a_status = (char *)calloc((size_t)(nvar > 0 ? nvar : 1), sizeof(char));
    *a_unsafe = (int *)calloc((size_t)(nvar > 0 ? nvar : 1), sizeof(int));
    if (nvar > 0 && (!*cell_idx || !*a_min || !*a_max || !*a_value ||
                     !*a_status || !*a_unsafe)) {
        TauAuditJjFree(nvar, *cell_idx, *a_min, *a_max, *a_value,
                       *a_status, *a_unsafe);
        *n_out = 0;
        goto error;
    }
    for (i = 0; i < nvar; i++) {
        (*cell_idx)[i] = var_cell[i];
        (*a_value)[i] = cells[var_cell[i]].value;
        (*a_status)[i] = cells[var_cell[i]].status;
    }

    /* load the LP: zero objective, maximise (toggled per cell below) */
    h = Highs_create();
    Highs_setBoolOptionValue(h, "output_flag", 0);
    Highs_setBoolOptionValue(h, "log_to_console", 0);

    st = Highs_passLp(h,
                      (HighsInt)nvar, (HighsInt)nrow, nnz,
                      kHighsMatrixFormatColwise, kHighsObjSenseMaximize,
                      0.0,
                      cost, clb, cub, rlo, rhi,
                      mstart, mindex, mvalue);
    if (st != kHighsStatusOk) {
        TauAuditJjFree(nvar, *cell_idx, *a_min, *a_max, *a_value,
                       *a_status, *a_unsafe);
        *n_out = 0;
        goto error;
    }

    /* per cell: minimise (cost -1 under MAX sense) and maximise */
    for (j = 0; j < nvar && !fail; j++) {
        st = Highs_changeColCost(h, (HighsInt)j, -1.0);
        if (st == kHighsStatusOk) {
            st = Highs_run(h);
        }
        if (st == kHighsStatusOk &&
            Highs_getModelStatus(h) == kHighsModelStatusOptimal) {
            (*a_min)[j] = -Highs_getObjectiveValue(h);
        }
        else {
            fail = 1;
        }

        if (!fail) {
            st = Highs_changeColCost(h, (HighsInt)j, 1.0);
            if (st == kHighsStatusOk) {
                st = Highs_run(h);
            }
            if (st == kHighsStatusOk &&
                Highs_getModelStatus(h) == kHighsModelStatusOptimal) {
                (*a_max)[j] = Highs_getObjectiveValue(h);
            }
            else {
                fail = 1;
            }
        }
        (void)Highs_changeColCost(h, (HighsInt)j, 0.0);
    }
    Highs_clear(h);

    if (fail) {
        /* legacy raises on the first non-optimal solve: abort the
           whole audit and release the result arrays */
        TauAuditJjFree(nvar, *cell_idx, *a_min, *a_max, *a_value,
                       *a_status, *a_unsafe);
        *n_out = 0;
        goto error;
    }

    /* legacy protection check:
         unsafe = (max - value) < UPL  or  (max - min) < LPL */
    for (i = 0; i < nvar; i++) {
        AuditCell *c = &cells[var_cell[i]];
        if ((*a_max)[i] - (*a_value)[i] < c->upl ||
            (*a_max)[i] - (*a_min)[i] < c->lpl) {
            (*a_unsafe)[i] = 1;
        }
    }

    *n_out = nvar;
    return 0;

error:
    free(colcount);
    free(mstart);
    free(mindex);
    free(mvalue);
    free(rlo);
    free(rhi);
    free(active_row);
    free(active);
    free(cost);
    free(cub);
    free(clb);
    free(var_cell);
    free(var_of);
    if (cells) {
        if (rels) {
            for (i = 0; i < nrel; i++) {
                free(rels[i].terms);
            }
            free(rels);
        }
        free(cells);
    }
    return -1;
}

void TauAuditJjFree(int n, long *cell_idx, double *a_min, double *a_max,
                    double *a_value, char *a_status, int *a_unsafe)
{
    (void)n;
    free(cell_idx);
    free(a_min);
    free(a_max);
    free(a_value);
    free(a_status);
    free(a_unsafe);
}
