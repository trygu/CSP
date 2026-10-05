/****************************************************/
/* The Cell Suppression Problem (CSP) for 3-dim     */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* CSPDEFNS.c                                       */
/* Version 0.0.1                                    */
/*                                                  */
/* Matteo FISCHETTI & Juan Jose SALAZAR GONZALEZ    */
/*                                                  */
/* Last modified May 20, 1997                       */
/****************************************************/

/* Definitions for SDC */

#define   VERSION      "1.0.0.4"  /* software current version                */

#define   CAPACITY  1
#define   BRIDGE    2
#define   COVER     3
#define   BRANCHCUT 4
#define   GOMORY    5



struct CELDA{
    int          index;
    signed char  coef;
    struct CELDA *next;
};

typedef struct{
    int    label;             /* label of the cell (for output purppose)   */
    struct CELDA *next;       /* pointer to a constraints where it is      */
    double nominal;           /* nominal value                             */
    char*  name;              /* cell name                                 */
    int    weight;            /* information value                         */
    int    index;             /* name of the cell                          */
    double lvalue;            /* lower possible value                      */
    double uvalue;            /* upper possible value                      */
    int    sensitive;         /* =0 iff corresponds to a NON sensitive cell*/
    double val;               /* value in the current fractional solution  */
    int    stat;              /* status of the variable                    */
    int    lp;                /* positive number of col if it is in the LP
                                 if stat=FIX then 0 iff fixed at root node */
} VARIABLE;

#define LP_LB   0
#define LP_UB   2
#define LP_BA   1
#define FIX_LB  -2
#define FIX_UB  -1
#define WAITING -3            /* variables outsite the LP but not fixed   */

typedef struct CON{
    int      index;           /* RHS cell for bridge constraint            */
    int      card;            /* number of non-zero elements               */
    VARIABLE **stack;         /* non-zero variables                        */
    double   *coef;           /* non-zero coefficients                     */
    double   rhs;             /* right-hand-side of the constraint         */
    char     sense;           /* 'G' , 'L'  or   'E'                       */
    int      type;            /* type of constraint                        */
    int      stat;            /* status of the constraint                  */
    int      lp;              /* positive number of row if it is in the LP */
    long int hash;            /* hash number for the row identification    */
    struct CON *con;          /* for saving original constraint in a COVER */
} CONSTRAINT;

struct BRANCH{
    double         val;       /* value of an initial bound                 */
    long           file;      /* position in '.bra' where it is described  */
    struct BRANCH  *next;     /* pointer to the next node in the tree      */
};

struct PRICE{
    VARIABLE       *col;      /* current column in PRICING structure       */
    struct PRICE   *next;     /* pointer to the next column in the struct. */
};

typedef struct{
    VARIABLE *var;            /* pointer to the cell                       */
    double   lpl;             /* lower protection level                    */
    double   upl;             /* upper protection level                    */
} SENSITIVE;

typedef struct {
    SENSITIVE *sen;           /* pointer to the cell                       */
    int       sense;          /* lower / upper                             */
    double    level;          /* protection level value                    */
    int       study;          /* 1=not known autoprotection                */
} PROT_LEVEL;


/* Setting output-file names */
//#define fsolution "cspSCIP.sol"   /* name of the file with the fractional sol. */
//#define fheuristi "cspSCIP.heu"   /* name of the file with the feasible sol.   */
//#define fout      "cspSCIP.sta"   /* name of the file for results              */
//#define fbranch   "cspSCIP.bra"   /* name of the file for saving bases         */

/* Setting parameters        */
//#define ZERO       1.0E-7   /* zero-epsilon                             */
//#define INF        1.0E+9   //2140000000   /* infinity                                 */
//#define MAX_TIME  18000000.0   /* maximum total CPU time                   */
//#define MAX_COLS_LP   10110  /* maximum number of columns in the LP      */
//#define MAX_ROWS_LP    4000  /* maximum number of cuts in the LP         */
//#define MAX_CUTS_POOL 500000  /* maximum number of cuts in the POOL       */
//#define MAX_CUTS_ITER   50   /* maximum number of new cuts per iteration */
//#define MIN_VIOLA    0.001   /* minimum violation for valid cuts         */
//#define MAX_SLACK    0.01    /* maximum slack for cuts in the LP         */



/* Mathematical abreviations */
#ifndef MIN
        #define   MIN(x,y)   ( (x) < (y) ? (x) : (y) )
#endif
#ifndef MAX
        #define   MAX(x,y)   ( (x) > (y) ? (x) : (y) )
#endif

/* Time function             */
float  seconds(void);