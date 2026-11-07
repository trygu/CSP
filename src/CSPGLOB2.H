/****************************************************/
/* The Cell Suppression Problem (CSP)               */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* SDCGLOB2.h                                       */
/* Version 0.0.1                                    */
/*                                                  */
/* Matteo FISCHETTI & Juan Jose SALAZAR GONZALEZ    */
/*                                                  */
/* Last modified 15 July 1995                       */
/****************************************************/
#include <string>
/* Global variables for UGFLP */

/* On Windows the DLL must export the four data globals the cut-hash
 * test references through this header: /EXPORT:*
 * (WINDOWS_EXPORT_ALL_SYMBOLS) exports the functions but skips the
 * data symbols, so the test link fails with LNK2019. The export
 * comes from the matching __declspec(dllexport) on the definitions
 * in Cspmain.c; TAUARGUS_CSP_BUILDING_DLL (set by CMake) is defined
 * only for the DLL's own translation units, never for the tests. */
#if defined(_WIN32)
  #if defined(TAUARGUS_CSP_BUILDING_DLL)
    #define CSP_GLOB_EXPORT __declspec(dllexport)
  #else
    #define CSP_GLOB_EXPORT
  #endif
#else
  #define CSP_GLOB_EXPORT
#endif

extern int            Rncells;    /* number of cells after reduction      */
extern int            Rnsums;     /* number of sums after reduction       */
extern double         *Rrhs;      /* RHS for sums after reduction         */

extern int            ncols;      /* number of cells in the table         */
extern VARIABLE       *columns;   /* potential variables                  */

CSP_GLOB_EXPORT extern int            nrows;      /* number of current constraints        */
CSP_GLOB_EXPORT extern CONSTRAINT     **rows;     /* potential constraints                */

extern int            mac;        /* number of current columns in the LP  */
extern VARIABLE       **cind;     /* real variables (in LP)               */

extern int            mar;        /* number of current rows in the LP     */
extern CONSTRAINT     **rind;     /* real constraints (in LP)             */

extern int            nsensitive; /* number of sensitive cells            */
extern SENSITIVE      *sensitive; /* list of sensitive cells              */

extern int            nprot_level;/* number of protection levels          */
extern PROT_LEVEL     *prot_level;/* list of protection levels            */


/* Globals data statistics */

extern int     iterations; /* number of iterations                 */
extern int     branchs;    /* number of studied branch-tree nodes  */
extern int     ntail;      /* number of iterations for tailing-off */

extern int     nsupport;   /* number of variables in 'support' (LB)*/
extern VARIABLE **support; /* variables with non-zero variable     */

extern struct PRICE *list_pricing;/* array to save WAITING columns */

extern int     nbetter;    /* number of variables in 'better'  (UB)*/
extern VARIABLE **better;  /* variables with non-zero variable     */

extern double  lowerb;     /* current lower bound                  */
extern int     upperb;     /* current upper bound (heuristic)      */

extern float   t0;         /* inicial clock-ticks CPU time         */
extern float   theur;      /* CPU time to for the heuristic        */
extern float   topti;      /* CPU time to find the current lowerb  */
extern char    ubtype;     /* UB type: find by Lp / find by Heuris */

extern int     cpool;      /* total number of POOL constraints     */
extern int     cbend;      /* total number of BENDERS constraints  */
extern int     ccapa;      /* total number of CAPACITY constraints */
extern int     cbrid;      /* total number of CAPACITY constraints */
extern int     ccove;      /* total number of COVER constraints    */
extern int     ccove2;     /* total number of COVER2 constraints   */
extern int     cgomo;      /* total number of GOMORY constraints   */
extern int     cmatc;      /* total number of COMB constraints     */
extern int     cpath;      /* total number of PATH constraints     */
extern int     cpath2;     /* total number of PATH2 constraints    */

extern int     bad_heuristic; /* flag =1 when no heuristic was applied */

// PWOF: 11-04-2013 added to control outputfilenames and variables
extern char* fsolution;
extern char* fheuristi;
extern char* fout;
extern char* fbranch;
extern char* fproblemlp;
extern char* fsdcnetlp;
extern char* fsdclp;
extern char* fpartial;
extern char* fmpsnet;
CSP_GLOB_EXPORT extern double ZERO;
extern double INF;
extern double MAX_TIME;
extern int MAX_COLS_LP;
extern int MAX_ROWS_LP;
CSP_GLOB_EXPORT extern int MAX_CUTS_POOL;
extern int MAX_CUTS_ITER;
extern double MIN_VIOLA;
extern double MAX_SLACK;
extern double FEAS_TOL;
extern double OPT_TOL;
extern char* logfile;
