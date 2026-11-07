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

/* nrows, rows, ZERO and MAX_CUTS_POOL are declared in CSPGLOB2.H
 * (with the Windows export macro) and defined in Cspmain.c; they are
 * not declared here so Cspmain.c sees exactly one declaration of each. */

int            Rncells;    /* number of cells after reduction      */
int            Rnsums;     /* number of sums after reduction       */
double         *Rrhs;      /* RHS for sums after reduction         */

int            nsensitive; /* number of sensitive cells            */
SENSITIVE      *sensitive; /* list of sensitive cells              */

int            nprot_level;/* number of protection levels          */
PROT_LEVEL     *prot_level;/* list of protection levels            */

int            ncols;      /* number of cells in the table         */
VARIABLE       *columns;   /* potential variables                  */

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
double INF;
double MAX_TIME;
int MAX_COLS_LP;
int MAX_ROWS_LP;
int MAX_CUTS_ITER;
double MIN_VIOLA;
double MAX_SLACK;
double FEAS_TOL;
double OPT_TOL;
char* logfile;
