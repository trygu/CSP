/*******************************************************************/
/*  Interface to LP-solvers:                                       */
/*                            CPLEX 3.0, CPLEX 5.0, XPRESS 10.0    */
/*                                                                 */
/*  J.J.Salazar                       Rotterdam, August 10, 1997   */
/*                                                                 */
/*  Interface to LP-solver: SCIP 3.0.1    (VSCIP)                  */
/*  M-S Hern�ndez-Garc�a              Spain, 2013                  */
/*******************************************************************/

#include "cspmain.h"
#include "cspback.h"

#ifndef CPLEX7
#define INFBOUND  1.0E+20
#endif

#define JJ_MIN                          1
#define JJ_MAX                         -1

#define JJ_OPTIMAL                      1
#define JJ_OPTIMAL_INFEAS               5  // 11  :  JJ changed on 25-oct-2010 because cplex 12 changed!
#define JJ_INFEASIBLE                   3  // 2 : Salome changed on 08-feb-2012
#define JJ_UNBOUNDED                    2  // 3 : Salome changed on 08-feb-2012
#define JJ_INForUNB                  1101

#define JJ_NETOPTIMAL                  -1
#define JJ_NETINFEASIBLE               -2
#define JJ_NETUNBOUNDED                -3

#ifdef VSCIP
#include <highs/interfaces/highs_c_api.h>

/* The HiGHS instance is an opaque void* in the C API. */
typedef void* JJLPptr;
#endif

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

int
   JJdualopt (JJLPptr),
   JJmipopt (JJLPptr),
   JJoptimize (JJLPptr);

int
   JJnetopt    (JJLPptr, int *, int *, int *, int *),
   JJhybnetopt (JJLPptr, char),
   JJembwrite  (JJLPptr, char *);

int
   JJcopyctype(JJLPptr, char *);

int
   JJpivotin (JJLPptr, int *, int);

int
   JJtwiddle (JJLPptr, int, int *);

int JJgetbhead     (JJLPptr, int *, double *),
   JJbinvcol      (JJLPptr, int, double *),
   JJbinvacol     (JJLPptr, int, double *),
   JJbinvrow      (JJLPptr, int, double *),
   JJbinvarow     (JJLPptr, int, double *),
   JJftran        (JJLPptr, double *),
   JJbtran        (JJLPptr, double *),
   JJgetijrow     (JJLPptr, int, int, int *),
   JJgetweight    (JJLPptr, int, int *, int *, double *, double *, int),
   JJstrongbranch (JJLPptr, int *, int, double *, double *, int);

void
   JJfreeprob     (JJLPptr*);

int
   JJreallocprob  (JJLPptr, double **, double **, char **, int **,
                 int **, double **, int **, double **, double **,
                 double **, char ***, char **, char ***, char **,
                 char **, int, int, int, unsigned, unsigned);

JJLPptr
   JJloadprob     (char *, int, int, int, int, double *, double *,
                 char *, int *, int *, int *, double *, double *,
                 double *, double *, int *, int *, int *, int *,
                 int *, double *, char *, char *, char *, char *,
                 char *, char **, char *, char **, char *, char **,
                 char *, int, int, int, int, int, unsigned, unsigned,
                 unsigned),
   JJloadlp       (char *,int, int, int, double *, double *, char *,
                 int *, int *, int *, double *, double *,double *,
                 double *, int, int, int),
   JJloadlpwnames (char *,int, int, int, double *, double *, char *,
                 int *, int *, int *, double *, double *,double *,
                 double *, char **, char *, char **, char *, int,
                 int, int, unsigned, unsigned);



void JJfreedata (char*,double*,double*,char*,
				   int*,int*,int*,double*,double*,double*,double*,
				   int*,int*,int*,int*,int*,double*,
				   char*,char*,char*,char*,char*,
				   char**,char*,char**,char*,char**,char*);


int
   JJlpread  (char *, int *, int *, int *, double **, double **,
            char **, int **, int **, int **, double **, double **,
            double **, char **, char **, char ***, char **,
            char ***, char **, int *, int *, int *, unsigned *,
            unsigned *),
   JJlpmread (char *, int *, int *, int *, double **, double **,
            char **, int **, int **, int **, double **, double **,
            double **, char **, char **, char ***, char **,
            char ***, char **, int *, int *, int *, unsigned *,
            unsigned *, char **);

int
   JJlpwrite   (JJLPptr, std::string),//char *),
   JJlprewrite (JJLPptr, char *);

int
   JJmbaseread  (char *, int, int, char **, char **, int *, int *),
   JJmsbaseread (char *, int, int, char **, char **, int *, int *,
               double **, double **),
   JJmbasewrite (JJLPptr, char *),
   JJvecread    (JJLPptr, char *),
   JJvecwrite   (JJLPptr, char *);

int
   JJmpsread  (char *, int *, int *, int *, int *, double **,
             double **, char **, int **, int **, int **,
             double **, double **, double **, double **, int **,
             int **, int **, int **, int **, double **,
             char **, char **, char **, char **, char **,
             char ***, char **, char ***, char **, char ***,
             char **, int *, int *, int *, int *, int *,
             unsigned *, unsigned *, unsigned *),
   JJmpsmread (char *, int *, int *, int *, int *, double **,
             double **, char **, int **, int **, int **,
             double **, double **, double **, double **, int **,
             int **, int **, int **, int **, double **,
             char **, char **, char **, char **, char **,
             char ***, char **, char ***, char **, char ***,
             char **, int *, int *, int *, int *, int *,
             unsigned *, unsigned *, unsigned *, char **,
             int *, int *, char **, int **, int **, double **);


int
   JJmpswrite   (JJLPptr, char *),
   JJmpsrewrite (JJLPptr, char *),
   JJdualwrite  (JJLPptr, char *, double *);

int
   JJgetspace   (JJLPptr, int *, int *, int *, unsigned *, unsigned *),
   JJloadbase   (JJLPptr, int *, int *),
   JJloadstart  (JJLPptr, int *, int *, double *, double *, double *, double *),
   JJloaddnorms (JJLPptr, double *, int *, int),
   JJloadpnorms (JJLPptr, double *, double *, int);

int
   JJgetprobname (JJLPptr, char *),
   JJgetmac      (JJLPptr),
   JJgetmar      (JJLPptr),
   JJgetmat      (JJLPptr),
   JJgetmae      (JJLPptr),
   JJgetnr       (JJLPptr),
   JJgetnrnz     (JJLPptr),
   JJgetenz      (JJLPptr),
   JJgetobjsen   (JJLPptr),
   JJgetobj      (JJLPptr, double *, int, int),
   JJgetrhs      (JJLPptr, double *, int, int),
   JJgetsense    (JJLPptr, char *, int, int),
   JJgetcols     (JJLPptr, int *, int *, int *, double *, int,
                int *, int, int),
   JJgetnewrhs (JJLPptr,double *,int,  int),
   JJgetrows     (JJLPptr, int *, int *, int *, double *, int,
                int *, int, int),
   JJgetnrows    (JJLPptr, int *, int *, int *, double *, int,
                int *, int, int),
   JJgetextra    (JJLPptr, int *, int *, int *, double *, int,
                int *, int, int),
   JJgetbdl      (JJLPptr, double *, int, int),
   JJgetbdu      (JJLPptr, double *, int, int),
   JJgetrngval   (JJLPptr, double *, int, int),
   JJgetnrind    (JJLPptr, int *, int, int),
   JJgetetype    (JJLPptr, int *, int, int),
   JJgetdataname (JJLPptr, char *),
   JJgetobjname  (JJLPptr, char *),
   JJgetrhsname  (JJLPptr, char *),
   JJgetrngname  (JJLPptr, char *),
   JJgetbndname  (JJLPptr, char *),
   JJgetcname    (JJLPptr, char **, char *, int, int *, int, int),
   JJgetrname    (JJLPptr, char **, char *, int, int *, int, int),
   JJgetename    (JJLPptr, char **, char *, int, int *, int, int),
   JJgetmacsz    (JJLPptr),
   JJgetmarsz    (JJLPptr),
   JJgetmatsz    (JJLPptr),
   JJgetmaesz    (JJLPptr),
   JJgetenzsz    (JJLPptr),
   JJgetbase     (JJLPptr, int *, int *),
   JJgetdnorms   (JJLPptr, double *, int *, int *),
   JJgetpnorms   (JJLPptr, double *, double *, int *),
   JJgetitc      (JJLPptr),
   JJgetitci     (JJLPptr),
   JJgetgrad     (JJLPptr, int, int *, double *),
   JJgetcoef     (JJLPptr, int, int, double *),
   JJgetsbcnt    (JJLPptr),
   JJgetijdiv    (JJLPptr, int *, int *);

unsigned
   JJgetcstorsz  (JJLPptr),
   JJgetrstorsz  (JJLPptr),
   JJgetestorsz  (JJLPptr);

double
   JJgetkappa    (JJLPptr);

int
   JJsolution  (JJLPptr, int *, double *, double *, double *,
                double *, double *),
   JJgetmethod (JJLPptr),
   JJgetobjval (JJLPptr, double *),
   JJgetstat   (JJLPptr),
   JJgetx      (JJLPptr, double *, int, int),
   JJgetpi     (JJLPptr, double *, int, int),
   JJgetslack  (JJLPptr, double *, int, int),
   JJgetdj     (JJLPptr, double *, int, int);

void
   JJchgobjsen  (JJLPptr, int);

int
   JJaddrows    (JJLPptr, int, int, int, double *, char *,
               int *, int *, double *, char **, char **),
   JJfaddrows   (JJLPptr, int, int, int, double *, char *,
               int *, int *, double *, char **, char **),
   JJdelrows    (JJLPptr, int, int),
   JJdelsetrows (JJLPptr, int *),
   JJaddcols    (JJLPptr, int, int, double *, int *, int *,
               double *, double *, double *, char **),
   JJdelcols    (JJLPptr, int, int),
   JJdelsetcols (JJLPptr, int *),
   JJchgcoef    (JJLPptr, int, int, double),
   JJchgbds     (JJLPptr, int, int *, char *, double *),
   JJchgobj     (JJLPptr, int, int *, double *),
   JJchgrhs     (JJLPptr, int, int *, double *),
   JJchgsense   (JJLPptr, int, int *, char *);

int
   JJsetaggfill  (int, int *, int *),
   JJsetitlim    (int, int *, int *),
   JJsetreinv    (int, int *, int *),
   JJsetmlim     (int, int *, int *),
   JJsetnlim     (int, int *, int *),
   JJsetnzlim    (int, int *, int *),
   JJsetelim     (int, int *, int *),
   JJsetenzlim   (int, int *, int *),
   JJsetedlimu   (int, int *, int *),
   JJsetepopt    (double, double *, double *),
   JJsetepper    (double, double *, double *),
   JJsettilim    (double, double *, double *),
   JJsetperind   (int),
   JJsetscr_ind  (JJLPptr, int),
   JJsetlogfile  (FILE *),
   JJsetitfoind  (int, int *, int *),
   JJsetppriind  (int, int *, int *),
   JJsetdpriind  (int, int *, int *),
   JJsetadvind   (int),
   JJsetobjulim  (double, double *, double *),
   JJsetobjllim  (double, double *, double *),
   JJsetsinglim  (int, int *, int *),
   JJsetsingtol  (double, double *, double *),
   JJseteprhs    (double, double *, double *),
   JJsetepmrk    (double, double *, double *),
   JJsetcraind   (int, int *, int *),
   JJsetscaind   (int, int *, int *),
   JJsetbasintv  (int, int *, int *),
   JJsetpreind   (int),
   JJsetdepind   (int),
   JJsetaggind   (int),
   JJsetxxxind   (int),
   JJsetfastmip  (int),
   JJsetrfilemul (double, double *, double *),
   JJsetcfilemul (double, double *, double *),
   JJsetdefaults (void),
   JJsetnetfind  (int, int *, int *),
   JJsetlpcallbackfunc (int (*)(JJLPptr, int));

int JJloadctype (JJLPptr, char*);
int JJmipoptimize (JJLPptr);
int JJgetmx (JJLPptr, double *, int, int);
/*
JJLPptr
   JJloadmprob   (char *, int, int, int, int, double *, double *,
                 char *, int *, int *, int *, double *, double *,
                 double *, double *, int *, int *, int *, int *,
                 int *, double *, char *, char *, char *, char *,
                 char *, char **, char *, char **, char *, char **,
                 char *, int, int, int, int, int, unsigned, unsigned,
                 unsigned, char *);
*/
int JJlpiterlimit(/*JJLPptr,*/int);
