/****************************************************/
/* The Cell Suppression Problem (CSP)               */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* SDCGLOB.h                                        */
/* Version 0.0.1                                    */
/*                                                  */
/* Matteo FISCHETTI & Juan Jose SALAZAR GONZALEZ    */
/*                                                  */
/* Last modified 15 July 1995                       */
/****************************************************/

/* Global variables for CSP */

/* On Windows the DLL must export the data symbols the tests reference
 * through CSPGLOB2.H; WINDOWS_EXPORT_ALL_SYMBOLS does not pick up BSS
 * globals on MSVC. Cspmain.c is the only includer of this header. */
#if defined(_WIN32)
#define CSP_GLOB_EXPORT __declspec(dllexport)
#else
#define CSP_GLOB_EXPORT
#endif

int            Rncells;    /* number of cells after reduction      */
int            Rnsums;     /* number of sums after reduction       */
double         *Rrhs;      /* RHS for sums after reduction         */

int            nsensitive; /* number of sensitive cells            */
SENSITIVE      *sensitive; /* list of sensitive cells              */

int            nprot_level;/* number of protection levels          */
PROT_LEVEL     *prot_level;/* list of protection levels            */

int            ncols;      /* number of cells in the table         */
VARIABLE       *columns;   /* potential variables                  */

/* extern: the definitions live in Cspmain.c. The plain declarations
 * used to create C common symbols, which MSVC does not export from
 * the DLL (uninitialized data is skipped by /EXPORT:*). */
CSP_GLOB_EXPORT extern int            nrows;      /* number of current constraints        */
CSP_GLOB_EXPORT extern CONSTRAINT     **rows;     /* potential constraints                */

int            mac;        /* number of current columns in the LP  */
VARIABLE       **cind;     /* real variables (in LP)               */

int            mar;        /* number of current rows in the LP     */
CONSTRAINT     **rind;     /* real constraints (in LP)             */



/* Globals data statistics */

int     iterations; /* number of iterations                 */
int     branchs;    /* number of studied branch-tree nodes  */
int     ntail;      /* number of iterations for tailing-off */

int     nsupport;   /* number of variables in 'support' (LB)*/
VARIABLE **support; /* variables with non-zero variable     */

struct PRICE *list_pricing; /* array to save WAITING columns     */

int     nbetter;    /* number of variables in 'better'  (UB)*/
VARIABLE **better;  /* variables with non-zero variable     */

double  lowerb;     /* current lower bound                  */
int     upperb;     /* current upper bound (heuristic)      */

float   t0;         /* inicial clock-ticks CPU time         */
float   theur;      /* CPU time to for the heuristic        */
float   topti;      /* CPU time to find the current lowerb  */
char    ubtype;     /* UB type: find by Lp / find by Heuris */

int     cpool;      /* total number of POOL constraints     */
int     cbend;      /* total number of BENDERS constraints  */
int     ccapa;      /* total number of CAPACITY constraints */
int     cbrid;      /* total number of CAPACITY constraints */
int     ccove;      /* total number of COVER constraints    */
int     ccove2;     /* total number of COVER2 constraints   */
int     cgomo;      /* total number of GOMORY constraints   */
int     cmatc;      /* total number of COMB constraints     */
int     cpath;      /* total number of PATH constraints     */
int     cpath2;     /* total number of PATH2 constraints    */

int     bad_heuristic; /* flag =1 when no heuristic was applied */

// PWOF: 11-04-2013 added to control outputfilenames and variables
char* fsolution;
char* fheuristi;
char* fout;
char* fbranch;
char* fproblemlp;
char* fsdcnetlp;
char* fCSPlog;
char* fsdclp;
char* fpartial;
char* fmpsnet;
CSP_GLOB_EXPORT extern double ZERO;
double INF;
double MAX_TIME;
int MAX_COLS_LP;
int MAX_ROWS_LP;
CSP_GLOB_EXPORT extern int MAX_CUTS_POOL;
int MAX_CUTS_ITER;
double MIN_VIOLA;
double MAX_SLACK;
double FEAS_TOL;
double OPT_TOL;
char* logfile;
