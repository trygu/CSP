/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPMAIN.c                                        */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/
//#define STAMP
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "Cspdefns.h"
#include "CSPGLOB.H"
#include "cspmain.h"
#include "cspprep.h"
#include "cspnet.h"
#include "cspbranc.h"
#include "cspheur.h"
#include "cspsep.h"
#include "cspcapa.h"
#include "cspsolve.h"
#include "cspprice.h"
#include "cspcard.h"
#include "cspback.h"
#include "cspdebug.h"
#include "CSPGLOB2.H"

/*#ifdef WIN32 //|| WIN64
    #include <direct.h>
#else
    #include <unistd.h>
#endif*/

std::streambuf *psbuf, *backup;
std::ofstream filestr;

/*  PROTOTYPES OF FUNCTIONS  */


extern int    l1u1,l1u0,l0u1,l0u0;
extern int    pricing_done;
extern struct BRANCH *tree;
extern float  tpricing,toptimize;


static int    change_status(char*);
static void   init_variable_type(void);
static void   init_constraint_pool(void);
static void   putcellsum(int,int,signed char *);
static void   delcell(int,double,double*);
static void   delsumcell(int,int,double,double*);
static void   libera(struct CELDA *);
static char*  copyname(char *);


static int    *newcellname; /* new number of an original cell     */
static int    *oldcellname; /* original number of an new cell     */
static int    *newsumname;  /* new number of an original sum      */
static int    *oldsumname;  /* original number of an new sum      */
static int    CSPcells;     /* original total number of cells     */
static int    CSPnsums;     /* original total number of sums      */

static double bestLB;       /* for saving the best lower bound    */
/*  FUNCTIONS  */
void    PPCSPSetFileNames(const char*);
void    PPCSPFreeFileNames();
void    PPCSPSetDoubleConstant(const int, double);
double  PPCSPGetDoubleConstant(const int);
void    PPCSPSetIntegerConstant(const int, int);
int     PPCSPGetIntegerConstant(const int);
int	PPCSPoptimize(IProgressListener* ProgressListener);
int     PPCSPloadprob(int,double*,int,double*,int*,char*,double*,double*,double*,double*,char**,int*,int*,signed char*);
int     PPCSPfreeprob();
int     PPCSPsolution(int*,int*,char*);
int     PPCSPrelbounds(int,int*,double*,double*,char);

// namespace SCIPv
#ifdef VSCIP
namespace SCIPv{
void *_scip;
void    CSPSetFileNames(const char* dir){ PPCSPSetFileNames(dir); }
void    CSPFreeFileNames(){ PPCSPFreeFileNames(); }
void    CSPSetDoubleConstant(const int VarNumber, double VarValue){ PPCSPSetDoubleConstant(VarNumber, VarValue); }
double  CSPGetDoubleConstant(const int VarNumber){ return PPCSPGetDoubleConstant(VarNumber); }
void    CSPSetIntegerConstant(const int VarNumer, int VarValue){ PPCSPSetIntegerConstant(VarNumer, VarValue); }
int     CSPGetIntegerConstant(const int VarNumber){ return PPCSPGetIntegerConstant(VarNumber); }
int	CSPoptimize(IProgressListener* ProgressListener){ return PPCSPoptimize(ProgressListener); }
int     CSPloadprob(int nsums_,double *rhs_,int ncells_,double *data_,int  *weight_,char *status_,double *lpl_,double *upl_,double *lb_,double *ub_,char **names_,int  *nlist_,int  *listcell_,signed char *listcoef_)
			{ return PPCSPloadprob(nsums_,rhs_,ncells_,data_,weight_,status_,lpl_,upl_,lb_,ub_,names_,nlist_,listcell_,listcoef_); }
int     CSPfreeprob(){ return PPCSPfreeprob(); }
int     CSPsolution(int *lowerb_,int *upperb_,char *status_){ return PPCSPsolution(lowerb_,upperb_,status_); }
int     CSPrelbounds(int nlist,int *list,double *ub,double *lb,char type){ return PPCSPrelbounds(nlist,list,ub,lb,type); }
}
#endif

// namespace CPLEX7
#ifdef CPLEX7
namespace CPLEXv{
EXPORTFUNC CPXENVptr Env;
void    CSPSetFileNames(const char* dir){ PPCSPSetFileNames(dir); }
void    CSPFreeFileNames(){ PPCSPFreeFileNames(); }
void    CSPSetDoubleConstant(const int VarNumber, double VarValue){ PPCSPSetDoubleConstant(VarNumber, VarValue); }
double  CSPGetDoubleConstant(const int VarNumber){ return PPCSPGetDoubleConstant(VarNumber); }
void    CSPSetIntegerConstant(const int VarNumer, int VarValue){ PPCSPSetIntegerConstant(VarNumer, VarValue); }
int     CSPGetIntegerConstant(const int VarNumber){ return PPCSPGetIntegerConstant(VarNumber); }
int	CSPoptimize(IProgressListener* ProgressListener){ return PPCSPoptimize(ProgressListener); }
int     CSPloadprob(int nsums_,double *rhs_,int ncells_,double *data_,int  *weight_,char *status_,double *lpl_,double *upl_,double *lb_,double *ub_,char **names_,int  *nlist_,int  *listcell_,signed char *listcoef_)
			{ return PPCSPloadprob(nsums_,rhs_,ncells_,data_,weight_,status_,lpl_,upl_,lb_,ub_,names_,nlist_,listcell_,listcoef_); }
int     CSPfreeprob(){ return PPCSPfreeprob(); }
int     CSPsolution(int *lowerb_,int *upperb_,char *status_){ return PPCSPsolution(lowerb_,upperb_,status_); }
int     CSPrelbounds(int nlist,int *list,double *ub,double *lb,char type){ return PPCSPrelbounds(nlist,list,ub,lb,type); }
}
#endif

// namespace XPRESS
#ifdef XPRESS_13
namespace XPRESSv{
void    CSPSetFileNames(const char* dir){ PPCSPSetFileNames(dir); }
void    CSPFreeFileNames(){ PPCSPFreeFileNames(); }
void    CSPSetDoubleConstant(const int VarNumber, double VarValue){ PPCSPSetDoubleConstant(VarNumber, VarValue); }
double  CSPGetDoubleConstant(const int VarNumber){ return PPCSPGetDoubleConstant(VarNumber); }
void    CSPSetIntegerConstant(const int VarNumer, int VarValue){ PPCSPSetIntegerConstant(VarNumer, VarValue); }
int     CSPGetIntegerConstant(const int VarNumber){ return PPCSPGetIntegerConstant(VarNumber); }
int	CSPoptimize(IProgressListener* ProgressListener){ return PPCSPoptimize(ProgressListener); }
int     CSPloadprob(int nsums_,double *rhs_,int ncells_,double *data_,int  *weight_,char *status_,double *lpl_,double *upl_,double *lb_,double *ub_,char **names_,int  *nlist_,int  *listcell_,signed char *listcoef_)
			{ return PPCSPloadprob(nsums_,rhs_,ncells_,data_,weight_,status_,lpl_,upl_,lb_,ub_,names_,nlist_,listcell_,listcoef_); }
int     CSPfreeprob(){ return PPCSPfreeprob(); }
int     CSPsolution(int *lowerb_,int *upperb_,char *status_){ return PPCSPsolution(lowerb_,upperb_,status_); }
int     CSPrelbounds(int nlist,int *list,double *ub,double *lb,char type){ return PPCSPrelbounds(nlist,list,ub,lb,type); }
}
#endif

/*  FUNCTIONS  */
// PWOF: 11-04-2013 added to control outputfilenames and variables

#define JJZERO          101
#define JJZERO1         102
#define JJZERO2         103
#define JJINF           104
#define JJMAXCOLSLP     105
#define JJMAXROWSLP     106
#define JJMAXCUTSPOOL   107
#define JJMAXCUTSITER   108
#define JJMINVIOLA      109
#define JJMAXSLACK      110
#define JJFEASTOL       111
#define JJOPTTOL        112
#define JJMAXTIME	113

void PPCSPSetFileNames(const char* dir)
{
        int maxfilenamesize = strlen(dir) + 15;
        int dirnamesize = strlen(dir)+1;
        
        fsolution = (char*) malloc(maxfilenamesize*sizeof(char));
	fheuristi = (char*) malloc(maxfilenamesize*sizeof(char));
	fout = (char*) malloc(maxfilenamesize*sizeof(char));
	fbranch = (char*) malloc(maxfilenamesize*sizeof(char));
        fproblemlp = (char*) malloc(maxfilenamesize*sizeof(char));
        fsdcnetlp = (char*) malloc(maxfilenamesize*sizeof(char));
        fsdclp = (char*) malloc(maxfilenamesize*sizeof(char));
        fCSPlog = (char*) malloc(maxfilenamesize*sizeof(char));
        fpartial = (char*) malloc(maxfilenamesize*sizeof(char));
        fmpsnet = (char*) malloc(maxfilenamesize*sizeof(char));
	strncpy(fsolution,dir,dirnamesize);
	strncat(fsolution,"cspSCIP.sol",15);
	strncpy(fheuristi,dir,dirnamesize);
	strncat(fheuristi,"cspSCIP.heu",15);
	strncpy(fout,dir,dirnamesize);
	strncat(fout,"cspSCIP.sta",15);
	strncpy(fbranch,dir,dirnamesize);
	strncat(fbranch,"cspSCIP.bra",15);
        strncpy(fproblemlp,dir,dirnamesize);
	strncat(fproblemlp,"problem.lp",15);
        strncpy(fsdcnetlp,dir,dirnamesize);
	strncat(fsdcnetlp,"sdcnet.lp",15);
        strncpy(fsdclp,dir,dirnamesize);
	strncat(fsdclp,"sdc.lp",15);
        strncpy(fpartial,dir,dirnamesize);
	strncat(fpartial,"partial.sol",15); 
        strncpy(fmpsnet,dir,dirnamesize);
	strncat(fmpsnet,"sdcnet.mps",15); 
        strncpy(fCSPlog,dir,dirnamesize);
        strncat(fCSPlog,"CSPlogfile.txt",15);
}

void PPCSPFreeFileNames()
{
	free(fsolution);
        fsolution = NULL; /*PWOF*/
	free(fheuristi);
        fheuristi = NULL; /*PWOF*/
	free(fout);
        fout = NULL; /*PWOF*/
	free(fbranch);
        fbranch = NULL; /*PWOF*/
        free(fproblemlp);
        fproblemlp = NULL; /*PWOF*/
        free(fsdcnetlp);
        fsdcnetlp = NULL; /*PWOF*/
        free(fsdclp);
        fsdclp = NULL; /*PWOF*/        
        free(fCSPlog);
        fCSPlog = NULL; /*PWOF*/
        free(fmpsnet);
        fmpsnet = NULL;
        free(fpartial);
        fpartial = NULL;
        
}

void PPCSPSetDoubleConstant(const int ConstName, double ConstValue)
{
	switch (ConstName){
		case JJZERO:
			ZERO = ConstValue;
			break;
		case JJINF:
			INF = ConstValue;
			break;
		case JJMAXTIME:
			MAX_TIME = ConstValue;
			break;
		case JJMINVIOLA:
			MIN_VIOLA = ConstValue;
			break;
		case JJMAXSLACK:
			MAX_SLACK = ConstValue;
			break;
                case JJFEASTOL:
                        FEAS_TOL = ConstValue;
                        break;
                case JJOPTTOL:
                        OPT_TOL = ConstValue;
                        break;
		default:
			std::cout << "Unknown constant " << ConstName << std::endl;
			break;
	}
}

double PPCSPGetDoubleConstant(const int ConstName)
{
	switch (ConstName){
		case JJZERO:
			return ZERO;
		case JJINF:
			return INF;
		case JJMAXTIME:
			return MAX_TIME;
		case JJMINVIOLA:
			return MIN_VIOLA;
		case JJMAXSLACK:
			return MAX_SLACK;
                case JJFEASTOL:
                        return FEAS_TOL;
                case JJOPTTOL:
                        return OPT_TOL;
		default:
			std::cout << "Unknown constant " << ConstName << std::endl;
			return -9;
	}
}

void PPCSPSetIntegerConstant(const int ConstName, int ConstValue)
{
	switch (ConstName){
        case JJMAXCOLSLP:
            MAX_COLS_LP = ConstValue;
            break;
        case JJMAXROWSLP:
            MAX_ROWS_LP = ConstValue;
            break;
        case JJMAXCUTSPOOL:
            MAX_CUTS_POOL = ConstValue;
            break;
        case JJMAXCUTSITER:
            MAX_CUTS_ITER = ConstValue;
            break;
		default:
			std::cout << "Unknown constant " << ConstName << std::endl;
			break;
	}
}

int PPCSPGetIntegerConstant(const int ConstName)
{
	switch (ConstName){
        case JJMAXCOLSLP:
            return MAX_COLS_LP;
        case JJMAXROWSLP:
            return MAX_ROWS_LP;
        case JJMAXCUTSPOOL:
            return MAX_CUTS_POOL;
        case JJMAXCUTSITER:
            return MAX_CUTS_ITER;
		default:
			std::cout << "Unknown constant " << ConstName << std::endl;
			return -9;
	}
}
// PWOF: added end

int  PPCSPoptimize(IProgressListener* ProgressListener)
{
    int        k,lcuts,upperb_root=0,upperb_init=0,nbr,old_ub,root,lprows,rowsinit,initpl;
    int        res_prep;
    CONSTRAINT *stack[MAX_CUTS_ITER];

    float      t1,tprep,theur0,theur1,troot,tseparation,tstart;
    double     lowerb_root;
    int        NEWcells,NEWnsums;
//#ifdef STAMP
    FILE       *fich;
    float      cpu;
//#endif
    t0=seconds();
    tstart=t0;
    
    root        = 1;
    nbr         = 0;
    tprep       = 0;
    tseparation = 0;
    theur0      = 0;
    theur1      = 0;
    lowerb_root = 0;
    lowerb      = 0;
    troot       = 0;
    bestLB      = 0;
    lprows      = 0;
    rowsinit    = 0;
    tree        = NULL;
    initpl      = 0;
    iterations  = 0;
    branchs     = 0;
    ntail       = 0;
    for(k=0;k<nsensitive;k++){
        if( sensitive[k].lpl>ZERO ) initpl++;
        if( sensitive[k].upl>ZERO ) initpl++;
    }

    // Added PWOF initializing global variables!
    // Needed when repeatedly calling PCSPoptimize() from modular approach
    nsupport = 0;
    nbetter = 0;
    upperb = 0;
    ubtype = ' ';
    cpool = 0;
    cbend = 0;
    ccapa = 0;
    cbrid = 0;
    ccove = 0;
    ccove2 = 0;
    cgomo = 0;
    cmatc = 0;
    cpath = 0;
    cpath2 = 0;
    // End Added PWOF
    /**** preprocessing ****/

    // k = preprocessing();
    // k replaced by res_prep: was confusing to use k as result of preprocessing as well as counter
    res_prep = preprocessing();

    tprep = seconds()-t0;
    NEWcells = Rncells;
    NEWnsums = Rnsums;
#ifdef STAMP
        std::cout << " New number of cells=" << Rncells << "     New number of links=" << Rnsums << std::endl;
#endif
        
    // k replaced by res_prep: was confusing to use k as result of preprocessing as well as counter
    if( res_prep==0 ){
        upperb_root = upperb_init = upperb = 0;
        nbetter = nsensitive;
        for(k=0;k<nbetter;k++){
            better[k] = sensitive[k].var;
            upperb += better[k]->weight;
        }
#ifdef STAMP
        std::cout << " All the sensitive cells are auto-protected " << std::endl;
#endif
        branchs = -1;
        ubtype = 'P';
        goto OUT;
    }
    if( res_prep<0 ){
        upperb_root = upperb_init = upperb = INF;
        nbetter = 1;
        better[0] = prot_level[-res_prep-1].sen->var;
#ifdef STAMP
        std::cout << " Not feasible solution exists because cell " << better[0]->name << std::endl;
#endif
        tprep = troot = seconds()-t0;
        branchs = -1;
        ubtype = 'P';
        goto OUT;
    }

#ifdef STAMP
    std::cout << " Time = " << seconds()-t0 << " sec.  " << std::endl;
#endif

    /**** initial feasible solution ****/
    for(k=0;k<Rncells;k++)
        if( columns[k].val < ZERO ) columns[k].val = 2*ZERO;

    if( !upperb ){
        upperb = INF;
        if( Rncells<1000 && Rnsums<1000) heuristic(1);
        else                             heuristic(0);
        theur0 = theur;
    }

    if( upperb == INF ){

#ifdef STAMP
        std::cout << "WARNING: not initial feasible solution found!!!!!!!!" << std::endl; //system("pause");
#endif
/**        nbetter = Rncells;
        upperb = 0;
        for(k=0;k<nbetter;k++){
            upperb += columns[k].weight;
            better[k] = columns+k;
        }
**/
        bad_heuristic = 1;
        goto OUT;
    }
    for(k=0;k<Rncells;k++) columns[k].val = 0.0;
    for(k=0;k<nsensitive;k++) sensitive[k].var->val = 1.0;


    upperb_root = upperb_init = upperb;
    ncols = Rncells;
    cpool = ccapa = cbend = cbrid = ccove = ccove2 = cmatc = cgomo = cpath = cpath2 = 0;
    rowsinit = nrows;

    load_branch_tree();

    while( read_prob() ){
        
        if (ProgressListener != NULL)
        {
                ProgressListener->UpdateUB(upperb);
                ProgressListener->UpdateLB(lowerb);
                ProgressListener->UpdateDiscrepancy(100*(upperb-lowerb)/lowerb); // percentages
                ProgressListener->UpdateTime((int) (seconds()-tstart));
        }
        
        bestLB = ceil(lowerb-ZERO);
        nbr++;
        ntail = 0;
        while( ceil(lowerb-ZERO) < upperb-ZERO ) {

            if( CSPstopcondition() ){   // Check whether an external stop
                unload_lp();
                goto OUT;
            }

            if( iterations%100==0 ) heuristic(0);

            t1 = seconds();
            lcuts = 0;
            separa(&lcuts,stack);
            tseparation += seconds()-t1;

            if ( lcuts==0 ) {
                 old_ub = upperb;
/**
                 if( branchs ) heuristic(0);
                 else          heuristic(1);
**/
                 if( Rncells<1000 && Rnsums<1000) heuristic(1);
                 else                             heuristic(0);



                 if( pricing_done==0 || upperb<old_ub )
                     if( pricing(1) ) goto JUMP;
                 if( integer_solution() ) break;
/****
                 if( branchs==0 )
                     insert_cutcard();
****/

/****
                 k = con_branching();
                 if((k==1) || (k==-1 && var_branching() ))
                     break;
****/
                 if ( var_branching() ) break;

            } else
                 add_rows(lcuts,stack);
JUMP:
            if(mar>lprows) lprows=mar;

            solve_lp(1);
            get_solution();
            del_cols();
            del_rows();

            t1 = seconds() - t0;
            iterations++;
#ifdef STAMP
            std::cout << "***  Iter = " << iterations << " ; Nrows = " << nrows << " ;";
            std::cout << " Time = " << t1 << " sec." << std::endl;
#endif
            if ( t1 > MAX_TIME ) 
            {
                int res = 0;

                if (ProgressListener != NULL) // Optimal approach (in case Modular approach, ProgressListener==NULL)
                {
                        float ttmp;
                        ttmp = seconds();
                        
                        ProgressListener->UpdateNSuppressed(nbetter);

                        res = CSPstoptime();
                        tstart = tstart + (seconds() - ttmp); // exclude time spent to ask for more time
                }
                
                if (res > 0)
                {
                    t0 = seconds();
                    MAX_TIME = res*60.0;  // res in minutes (integer), MAX_TIME in seconds (double))
                }
                else
                {
                    unload_lp();
                    goto OUT;
                }
            }
        }
        if(root){
            troot       = seconds()-t0;
            theur1      = theur;
            lowerb_root = lowerb;
            upperb_root = upperb;
            if( lowerb_root>upperb_root ) lowerb_root = upperb_root;
            root = 0;
        }
        unload_lp();
    }

OUT:

    if( tree ) bestLB = ceil( tree->val - ZERO );
//    else if(lowerb-ZERO<upperb) bestLB = ceil( lowerb - ZERO );
//    else bestLB = upperb;
    unload_branch_tree();



//#ifdef STAMP

    cpu = seconds() - t0;
    if ( cpu>MAX_TIME ){
        std::cout << "TIME LIMIT EXCEEDED: SDC best solution = ";
        cpu = -MAX_TIME;
    }else{
        std::cout << "END: SDC optimal solution = ";
        bestLB = upperb;
    }
    std::cout << upperb << " with " << nbetter << " suppressions " << std::endl;
    std::cout << " Total time = " << cpu << std::endl;

    fich=fopen(fout,"r");
    if(fich==NULL){
       fich=fopen(fout,"w");
       if(fich==NULL){
           std::cout << "ERROR: it not possible to write in disk" << std::endl;
           CSPexit(EXIT_ERROR); //exit(1);
       }
       fprintf(fich,"cells sums cells sums |P| pl0 l1u1 l1u0 l0u1 l0u0  tprep  theur0   UB0  theur1  theur*    UB1     LB1     time   opt   SUP   time  BestLB  nodes    Ttime   Tsepar  Tprice   Toptim  maxmar cuts  pool  bender capa  bridge cover cover2 gomo comb  path path2 iterations  UBtype\n\n");
    }
    fclose(fich);

    fich=fopen(fout,"a");
    fprintf(fich,"%3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %8.1f %8.1f %8d %8.1f %8.1f %8d %8.1f %8.1f %8d %4d %8.1f %5d %4d %8.2f %8.1f %8.1f %8.1f %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d   %5d    %c\n",
       CSPcells,CSPnsums,NEWcells,NEWnsums,nsensitive,initpl,l1u1,l1u0,l0u1,l0u0,tprep,theur0,upperb_init,theur1,theur,upperb_root,(float)ceil((double)lowerb_root-ZERO),troot,upperb,nbetter-nsensitive,topti,
       (int)bestLB,nbr,cpu,tseparation,tpricing,toptimize,lprows,nrows,cpool,rowsinit,ccapa,cbrid,ccove,ccove2,cgomo,cmatc,cpath,cpath2,iterations,ubtype);
    fclose(fich);
//#endif 
    
    //return(0);
    return(bad_heuristic);
}

/*******************************************************************/

static struct CELDA **listcell,**listsum;
static int          *nlistcell,*nlistsum;
std::string filenameA;   // Salome added 31/01/2014

int PPCSPloadprob(int nsums_,double *rhs_,int ncells_,double *data_,int  *weight_,char *status_,double *lpl_,double *upl_,double *lb_,double *ub_,char **names_,int  *nlist_,int  *listcell_,signed char *listcoef_)

{
    int      h,k,l,nz,nl,np;
    VARIABLE *col;
    struct   CELDA *c;

    //filenameA = _getcwd(NULL, 0);
    //filenameA.append("\\logfile.txt");
    //filestr.open (filenameA.c_str());//logfile);     // Salome modified 31/01/2014
    filestr.open (fCSPlog);     /*PWOF to write logfile in temp-directory 03-02-2014*/

    backup = std::cout.rdbuf();     // back up cout's streambuf

    psbuf = filestr.rdbuf();   // get file's streambuf
    std::cout.rdbuf(psbuf);         // assign streambuf to cout

/**********************
#ifdef PARTIAL
    double   max2;

    printf("changing input data for partial suppression results\n");
    max2 = 0.0;
    for(k=0;k<ncells_;k++)
        if( data_[k]>max2 ) max2 = data_[k];

//    srand(12345);
    for(k=0;k<ncells_;k++){
//        range = (int)floor( data_[k]/2.0 );

        lb_[k] = data_[k]/2.0;
        ub_[k] = 2.0*data_[k];
//        ub_[k] = ( range>0 ? data_[k]+range+ (rand() % range) : 2 );


        weight_[k] = ceil( ub_[k]-lb_[k] );
        if(max2>1000000) weight_[k] = ceil( (double)weight_[k]/1000.0 );

        if( upl_[k]>ZERO || lpl_[k]>ZERO )
            upl_[k] = lpl_[k] = data_[k]/4.0 ;
    }
#endif
**********************/


    nl=l=0;
    for(k=0;k<nsums_;k++){
        nz = 0;
        h = nlist_[k];
        while(h--){
            if( status_[ listcell_[l] ]!='z' ){
                nz++;
                np = nl;
                nl = listcell_[l];
            }
            l++;
        }
        if( nz==1 ){
            if( status_[nl]=='u' ){
#ifdef STAMP
                std::cout << "Unfeasible problem because sensitive cell " << names_[nl] << " needs to be published according to constraint " << k << std::endl;
#endif
                return 1;
            } else {
                status_[nl]='z';
            }
        } else if ( nz==2 && status_[nl] != status_[np] ){
//			status_[nl] = status_[np] = 'u';
#ifdef STAMP
            if(status_[nl]=='s')
                std::cout << "New sensitive cell: " << names_[nl] << std::endl;
            else
                std::cout << "New sensitive cell: " << names_[np] << std::endl;
#endif
        }

    }




#ifdef STAMP
    if( CSPtestprob(nsums_,rhs_,ncells_,data_,weight_,status_,lpl_,upl_,lb_,ub_,names_,nlist_,listcell_,listcoef_) )
        return 1;
#endif

    CSPcells = ncells_;
    CSPnsums = nsums_;

    listcell = (struct CELDA **)calloc( ncells_ , sizeof(struct CELDA *) );
    if( listcell==NULL )return(1);
    listsum  = (struct CELDA **)calloc( nsums_ , sizeof(struct CELDA *) );
    if( listsum==NULL )return(1);
    nlistcell = (int *)calloc( ncells_ , sizeof(int) );
    if( nlistcell==NULL )return(1);
    nlistsum = (int *)calloc( nsums_ , sizeof(int) );
    if( nlistsum==NULL )return(1);

    Rrhs = (double *)malloc( nsums_ * sizeof(double) );
    if( Rrhs==NULL )return(1);



    l = 0;
    for(k=0;k<nsums_;k++){
        h = nlist_[k];
        while(h--){
            putcellsum( k , listcell_[l] , &listcoef_[l] );
            l++;
        }
    }
    for(k=0;k<ncells_;k++)
        if( status_[k]=='z' )
            delcell(k,data_[k],rhs_);




#ifdef STAMP
    l = 0;
    for(k=0;k<nsums_;k++)
        if( nlistsum[k]==2 )
            l++;
    std::cout << "Warning: " << l << " links of TWO variables" << std::endl;


//            printf("WARNING: condition %d has %d variables\n",k,nlistsum[k]);
//            h = nlist_[k];
//            while(h--){
//                printf(" %d*A%d(%c)",listcoef_[l],listcell_[l],status_[listcell_[l]]);
//                l++;
//            }
//            printf("=%f\n",rhs_[k]);
//        } else
//            l += nlist_[k];
//    }
#endif




    /******  names ******/

    newcellname = (int *)malloc( ncells_ * sizeof(int) );
    if( newcellname==NULL )return(1);
    oldcellname = (int *)malloc( ncells_ * sizeof(int) );
    if( oldcellname==NULL )return(1);
    newsumname = (int *)malloc( nsums_ * sizeof(int) );
    if( newsumname==NULL )return(1);
    oldsumname = (int *)malloc( nsums_ * sizeof(int) );
    if( oldsumname==NULL )return(1);

    Rncells = 0;
    for(k=0;k<ncells_;k++)
       if( nlistcell[k] ){
           newcellname[k] = Rncells;
           oldcellname[Rncells] = k;
           Rncells++;
       }else
           newcellname[k] = -1;

    Rnsums = 0;
    for(k=0;k<nsums_;k++)
       if( nlistsum[k] ){
           newsumname[k] = Rnsums;
           oldsumname[Rnsums] = k;
           Rrhs[Rnsums] = rhs_[k];
           Rnsums++;
       }else
           newsumname[k] = -1;


     /* allocating memory */

     cind    = (VARIABLE **)malloc( (1+MAX_COLS_LP) * sizeof(VARIABLE *) );
     if( cind==NULL ){
         CSPexit(EXIT_MEMO);
         //exit(1);
     }
     rind    = (CONSTRAINT **)malloc( (1+MAX_ROWS_LP) * sizeof(CONSTRAINT *) );
     if( rind==NULL ){
         CSPexit(EXIT_MEMO);
         //exit(1);
     }
     columns = (VARIABLE *)malloc( Rncells * sizeof(VARIABLE) );
     if( columns==NULL ){
         CSPexit(EXIT_MEMO);
         //exit(1);
     }
     rows    = (CONSTRAINT **)malloc( MAX_CUTS_POOL * sizeof(CONSTRAINT *) );
     if( rows==NULL ){
         CSPexit(EXIT_MEMO);
         //exit(1);
     }
     support = (VARIABLE **)malloc( Rncells * sizeof(VARIABLE *) );
     if( support==NULL ){
         CSPexit(EXIT_MEMO);
         //exit(1);
     }

     nbetter = 0;
     better = (VARIABLE **)malloc( Rncells * sizeof(VARIABLE *) );
     if( better==NULL )return(1);


/************  loading table in CBS input format *******************/

     for(k=0;k<Rncells;k++){
         l = oldcellname[k];
         col = columns + k;
         col->index     = k;
         col->nominal   = data_[l];
         col->name      = (names_ ? copyname( names_[l] ) : NULL );
         col->sensitive = change_status( &status_[l] );
         col->weight    = weight_[l];
         col->uvalue    = ub_[l] - data_[l];
         col->lvalue    = data_[l] - lb_[l];
         col->index     = k;
         col->next      = NULL;
         col->val       = 0.0;
         col->lp        = 0;
         col->next      = listcell[l];
         for( c=listcell[l] ; c ; c=c->next )
             c->index = newsumname[ c->index ];
     }

     /* other cell details */

     for(k=0;k<Rncells;k++){
            col = columns + k;
            switch( col->sensitive ){
                case 0:
                          col->sensitive = 0;
                          col->stat      = LP_LB;
                          break;
                case 2:
                          col->sensitive = 1;
#ifdef PARTIAL
                          col->stat      = LP_LB;
                          col->val       = 0.5;
#else
                          col->stat      = FIX_UB;
                          col->val       = 1.0;
#endif
                          nsensitive++;
                          break;
                case 9:
                          col->sensitive = 0;
                          col->stat      = FIX_LB;
                          col->lvalue = col->uvalue = 0;
                          break;
                default:
                          std::cout << "ERROR: '" << col->sensitive << "' unknown type of cell" << std::endl;
                          return(1);
            }
     }

     /*  primary cells */

     if( nsensitive==0 ) return(0);
     sensitive  = (SENSITIVE  *)malloc( nsensitive * sizeof(SENSITIVE) );
     prot_level = (PROT_LEVEL *)malloc( 2*nsensitive * sizeof(PROT_LEVEL) );    // nivel de proteccion
     if( sensitive==NULL && prot_level==NULL ) return(1);

     nprot_level = 0;
     nsensitive  = 0;
     for(k=0;k<Rncells;k++){
            l = oldcellname[k];
            col = columns + k;
            if( col->sensitive ){
                 sensitive[ nsensitive ].var = col;
                 sensitive[ nsensitive ].upl = upl_[l] /*- data_[l]*/;
                 if( upl_[l] ){
                     prot_level[ nprot_level ].sen   = sensitive + nsensitive;
                     prot_level[ nprot_level ].sense = 1;
                     prot_level[ nprot_level ].level = upl_[l];
                     prot_level[ nprot_level ].study = 1;
                     nprot_level++;
                 }
                 sensitive[ nsensitive ].lpl = /* data_[l] -*/ lpl_[l];
                 if( lpl_[l] ){
                     prot_level[ nprot_level ].sen   = sensitive + nsensitive;
                     prot_level[ nprot_level ].sense = -1;
                     prot_level[ nprot_level ].level = lpl_[l];
                     prot_level[ nprot_level ].study = 1;
                     nprot_level++;
                 }
                 nsensitive++;
                 col->sensitive = nsensitive;
            }
     }




    /* initializing variable labels and constraint pool */

    init_variable_type();
    init_constraint_pool();

    free((void *)listcell);
    listcell = NULL; /*PWOF*/
    free((void *)nlistcell );
    nlistcell = NULL; /*PWOF*/
    for(k=0;k<nsums_;k++)
        libera( listsum[k] );
    free((void *)listsum);
    listsum = NULL; /*PWOF*/
    free((void *)nlistsum );
    nlistsum = NULL; /*PWOF*/
    return(0);
}

static void putcellsum(int  sum ,int  cell ,signed char   * coef )

{
    struct CELDA *c;

    c = (struct CELDA *)malloc( sizeof(struct CELDA) );
    if( c==NULL ){        
        std::cout << "ERROR: not enought memory for CELDA" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    c->index = cell;
    c->next  = listsum[ sum ];
    c->coef  = *coef;
    listsum[ sum ] = c;
    nlistsum[ sum ]++;

    c = (struct CELDA *)malloc( sizeof(struct CELDA) );
    if( c==NULL ){        
        std::cout << "ERROR: not enought memory for CELDA" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    c->index = sum;
    c->next  = listcell[ cell ];
    c->coef  = *coef;
    listcell[ cell ] = c;
    nlistcell[ cell ]++;
}

static void delcell(int  cell ,double  data ,double * rhs )

{
    struct CELDA *c,*c0;
    c0 = NULL;
    c  = listcell[cell];
    while(c){   
        delsumcell( c->index , cell , data , rhs );
        c0 = c;
        c  = c->next;
        free( c0 );
        c0 = NULL; /*PWOF*/
    }
    listcell[ cell ] = NULL;
}

static void delsumcell(int  sum ,int  cell ,double  data ,double * rhs )

{
    struct CELDA *c,*c0;
    c0 = NULL;
    c  = listsum[sum];
    while(c){
        if( c->index == cell ){
            if ( c0 ) c0->next     = c->next;
            else      listsum[sum] = c->next;
            rhs[ sum ] -= c->coef * data;
            free( c );
            c = NULL; /*PWOF*/
            nlistcell[cell]--;
            nlistsum[sum]--;
            return;
        }
        c0 = c;
        c  = c->next;
    }    
    std::cout << "ERROR: cell not found!" << std::endl;
    CSPexit(EXIT_ERROR); //exit(1);
}



static int change_status(char * ccc )

{
     switch( *ccc ){
        case 's': return(0);
        case 'u': return(2);
        case 'z': return(9);
        default :
              std::cout << " ERROR: input cell with status " << *ccc << std::endl;
            break;
     }
     return(-1);
}


static void libera(struct CELDA *ptr)

{
     if(ptr) libera(ptr->next);
     free(ptr);
     ptr = NULL; /*PWOF*/
}
 
static char *copyname(char *name)

{
    char *ptr;

    ptr = (char *)malloc( (strlen(name)+1)*sizeof(char) );
    if( ptr==NULL ){
        std::cout << "ERROR: not memory for names" << std::endl;
        return NULL;
    }
    strcpy(ptr,name);
    return ptr;
}
 

/*********** fixing variables outside the LP *****************/

static void  init_variable_type()
{
    int      i,j,cont,cant,num;
    VARIABLE *col;
    VARIABLE **stack;
    int      *mrow;
    int      sort_cols(const void*,const void*);
    struct   CELDA *c;

    cont = 0;
    for(i=0;i<Rncells;i++)
        if ( columns[i].stat == LP_LB ) cont++;
    if( cont<MAX_COLS_LP ) return;

    for(i=0;i<Rncells;i++)
        if ( columns[i].stat == LP_LB ) columns[i].stat = WAITING;

    stack = (VARIABLE**)malloc( Rncells * sizeof(VARIABLE*) );
    if( stack==NULL ){        
        std::cout << "ERROR: not enough memory for INIT VARIABLE TYPe" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    mrow = (int*)calloc( Rnsums , sizeof(int) );
    if( mrow==NULL ){        
        std::cout << "ERROR: not enough memory for INIT VARIABLE TYPe" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    for(i=0;i<nsensitive;i++)
        if( sensitive[i].lpl || sensitive[i].upl )
            for( c=sensitive[i].var->next ; c ; c=c->next )
                mrow[ c->index ] = 1;

    num = 0;
    for(i=0;i<Rnsums;i++)
        if( mrow[i] )
            num++;

    cant = MAX_COLS_LP / num;
    num = 0;
    for(i=0;i<Rnsums;i++){
        if( mrow[i] ){
            cont = 0;
            for( c=listsum[ oldsumname[i] ] ; c ; c=c->next ){
                col = columns + newcellname[c->index];
                if( col->stat == WAITING )
                    stack[ cont++ ] = col;
            }
            if( cont>cant ){
                //qsort( (char *)stack , cont , sizeof(VARIABLE *) , sort_cols );
                qsort( (void *)stack , cont , sizeof(VARIABLE *) , sort_cols );
                cont = cant;
            }
            for(j=0;j<cont;j++)
                stack[j]->stat = LP_LB;
            num += cont;
        }
    }
    free(stack);
    //stack = NULL; /*PWOF*/
    free(mrow);
    //mrow = NULL; /*PWOF*/
}

int sort_cols(const void*p,const void*q/*VARIABLE **p,VARIABLE **q*/)

{
    int cp,cq;

    cp = (*(VARIABLE **)p)->weight;
    cq = (*(VARIABLE **)q)->weight;
    if( cp<cq ) return(-1);
    if( cp>cq ) return(1);
    return(0);
}


/*********** inserting Basic Capacity Cuts in the POOL *****************/

static void  init_constraint_pool()
{
    int      i,nvar,sense;
    double   rhs,maxlhs;
    VARIABLE *var,*senvar;
    VARIABLE **xvar;
    CONSTRAINT *con;
    double   *cvar;
    struct   CELDA *c1,*c2;

    xvar = (VARIABLE**)malloc( Rncells * sizeof(VARIABLE*) );
    if( xvar==NULL ){        
        std::cout << "ERROR: not enough memory for INIT POOL TYPe" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    cvar = (double*)malloc( Rncells * sizeof(double) );
    if( cvar==NULL ){        
        std::cout << "ERROR: not enough memory for INIT POOL TYPe" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    nrows = 0;
    for(i=0;i<nprot_level;i++){
        senvar = prot_level[i].sen->var;
        sense  = prot_level[i].sense;
        rhs    = prot_level[i].level;
        for( c1=senvar->next ; c1 ; c1=c1->next ){
            maxlhs = rhs;
            nvar = 0;
            for( c2=listsum[ oldsumname[c1->index] ] ; c2 ; c2=c2->next ){
                var = columns + newcellname[c2->index];
                if( var != senvar ){
                    xvar[nvar] = var;
                    if( c2->coef==c1->coef )
                        cvar[nvar]=(sense==1 ? var->lvalue : var->uvalue);
                    else
                        cvar[nvar]=(sense==1 ? var->uvalue : var->lvalue);
                    if( var->sensitive ) maxlhs -= cvar[nvar];
                    if( cvar[nvar]>ZERO ) nvar++;
                }
            }

            if(nvar && maxlhs>ZERO){
                con = capacity_constraint(nvar,xvar,cvar,senvar,rhs);
                if( !new_row(con) ) remove_row(con);
                else{
                    rows[nrows-1]->stat = FIX_LB;
                    rows[nrows-1]->lp   = -1;
                }
                prot_level[i].study = 0;
            }            
        }
    }
    free(xvar);
    xvar = NULL; /*PWOF*/
    free(cvar);
    cvar = NULL; /*PWOF*/
}

int PPCSPfreeprob()
{
     int k;
     for(k=0;k<Rncells;k++){
         free( columns[k].name );
         columns[k].name = NULL; /*PWOF*/
         libera( columns[k].next );
     }
     free((void *)columns);
     columns = NULL; /*PWOF*/
     for(k=0;k<nrows;k++)
         remove_row( rows[k] );
     free((void *)rows);
     rows = NULL; /*PWOF*/
     free((void *)sensitive);
     sensitive = NULL; /*PWOF*/
     free((void *)prot_level);
     prot_level = NULL; /*PWOF*/
     free((void *)better);
     better = NULL; /*PWOF*/
     free((void *)newcellname);
     newcellname = NULL; /*PWOF*/
     free((void *)oldcellname);
     oldcellname = NULL; /*PWOF*/
     free((void *)newsumname);
     newsumname = NULL; /*PWOF*/
     free((void *)oldsumname);
     oldsumname = NULL; /*PWOF*/
     free((void *)Rrhs);
     Rrhs = NULL; /*PWOF*/
     free((void *)support);
     support = NULL; /*PWOF*/
     free((void *)cind);
     cind = NULL; /*PWOF*/
     free((void *)rind);
     rind = NULL; /*PWOF*/
     return 0;
}


int PPCSPsolution(int *lowerb_,int *upperb_,char *status_)
{
    int k;

    if(lowerb_)
        *lowerb_ = (int)ceil( bestLB-ZERO );
    if(upperb_)
        *upperb_ = upperb;
    if(status_){
        for(k=0;k<nbetter;k++)
            if( !better[k]->sensitive )
                status_[ oldcellname[better[k]->index] ] = 'm';
        if( upperb==INF )
                status_[ oldcellname[better[0]->index] ] = 'x';
    }
    write_heu(fheuristi);
    return 0;
}


int PPCSPrelbounds(int nlist,int *list,double *ub,double *lb,char type)
 //type 'C' = real numbers ; 'I' = integer numbers
{
    int    k,l;
    double *status;

    if( nbetter==0 ) return 0;
    status = (double *)calloc(sizeof(double),Rncells);
    if(status==NULL){
        std::cout << "There is not enough memory for STATUS" << std::endl;
        return 1;
    }
    for(k=0;k<nbetter;k++) status[better[k]->index] = 1;

    load_network(status,type);
    free(status);
    status = NULL; /*PWOF*/

    for(k=0;k<nlist;k++){
        l = newcellname[list[k]];
        if( l != -1){
            ub[k] =  protection_level(columns+l, 1,NULL,NULL,NULL,NULL,type);
            lb[k] = -protection_level(columns+l,-1,NULL,NULL,NULL,NULL,type);
        }
    }
    unload_network();
    return 0;
}


/* to "audit" the current partial solution */
int      CSPpartialbounds()
{
    int    k;
    double ub,lb;
    double *status;

    std::cout << " AUDIT: partial missing solution:" << std::endl;

    status = (double *)calloc(sizeof(double),Rncells);
    if(status==NULL){
        std::cout << "There is not enough memory for STATUS" << std::endl;
        return 1;
    }
    for(k=0;k<Rncells;k++)
        status[k] = columns[k].val;

    load_network(status,'C');

    for(k=0;k<Rncells;k++){
        if( status[k]>ZERO ){
            ub =  protection_level(columns+k, 1,NULL,NULL,NULL,NULL,'C');
            lb = -protection_level(columns+k,-1,NULL,NULL,NULL,NULL,'C');
            std::cout << " " << columns[k].name << " : " << lb << "  " << columns[k].nominal << "  " << ub << std::endl;
        }
    }
    unload_network();
    free(status);
    status = NULL; /*PWOF*/
    return 0;
}





int CSPwrite(char *filename)

{
    int i,j;
    FILE *file;
    struct   CELDA *c,*d;
    VARIABLE *col;
	struct CELDA **listsum;  
    int *nlistsum;



    file = fopen(filename,"w");
    if( file==NULL ){
        std::cout << "ERROR: not possible to open file " << filename << std::endl;
        return 1;
    }
    fprintf(file,"0\n %d\n",Rncells);
    for(i=0;i<Rncells;i++){
        col = columns+i;
        // if( col->name )
        //    fprintf(file, " %10s ",col->name);
        // else
        fprintf(file, " %10d ",col->index);
        fprintf(file,"%10.1f %10d ",col->nominal,col->weight);
        switch( col->stat ){
            case FIX_UB: fprintf(file,"u ");
                    break;
//            case FIX_LB: fprintf(file,"z ");
//                    break;
            default : fprintf(file,"s ");
        }
        fprintf(file,"%10.1f %10.1f ",(col->nominal)-(col->lvalue),(col->nominal)+(col->uvalue));
        if( col->sensitive )
            fprintf(file,"%10.1f %10.1f 0.0\n",
                sensitive[col->sensitive-1].lpl,
                sensitive[col->sensitive-1].upl);
        else
            fprintf(file,"%10.1f %10.1f 0.0\n",0.0,0.0);
    }

    

    listsum  = (struct CELDA **)calloc( Rnsums , sizeof(struct CELDA *) );
    if( listcell==NULL )return(1);
    nlistsum = (int *)calloc( Rnsums , sizeof(int) );
    if( nlistsum==NULL )return(1);
    for(i=0;i<Rnsums;i++)
		nlistsum[i] = 0;
    
    for(i=0;i<Rncells;i++){
        col = columns+i;
        for( c=col->next ; c ; c=c->next ){
             d = (struct CELDA *)malloc( sizeof(struct CELDA) );
             if( d==NULL ){
                 std::cout << "ERROR: not enought memory for CELDA" << std::endl;
                 CSPexit(EXIT_MEMO); //exit(1);
			 }
			 j=c->index;
             d->index = i;
             d->next  = listsum[j];
             d->coef  = c->coef;
             listsum[j] = d;
             nlistsum[j]++;
		}
	}
	
    fprintf(file,"%d\n",Rnsums);
    for(i=0;i<Rnsums;i++){
        fprintf(file,"%5.1f %5d :",Rrhs[i],nlistsum[i]);
		for( d=listsum[i] ; d ; d=d->next )
			fprintf(file," %d (%d)",d->index,d->coef);
		fprintf(file,"\n");
	}

    for(i=0;i<Rnsums;i++){
		d=listsum[i];
		while(d){
			c = d->next; 
			free(d);
			d = c;
		}
	}
	free(listsum);
        listsum = NULL; /*PWOF*/
	free(nlistsum);
        nlistsum = NULL; /*PWOF*/

	

	/**
	fprintf(file,"For each cell, the equations (and coeficients) are:\n");
    for(i=0;i<Rncells;i++){
        col = columns+i;
        j = 0;
        for( c=col->next ; c ; c=c->next ) j++;
        fprintf(file,"cell %d in %d equations :",i,j);
        for( c=col->next ; c ; c=c->next )
                fprintf(file," %d(%d)",c->index,c->coef);
        fprintf(file,"\n");
    }
    **/
    fclose(file);
    return 0;
}



int CSPtestprob(int nsums_,double *rhs_,int ncells_,double *data_,int  *weight_,char *status_,double *lpl_,double *upl_,double *lb_,double *ub_,char **names_,int  *nlist_,int  *listcell_,signed char *listcoef_)

{
    int      s,k,l;
    double   val;

    std::cout << " input checking ... ";
    if( nsums_<1 ) {
        std::cout << "ERROR: number of constraints=" << nsums_ << std::endl;
        return 1;
    }
    if( ncells_<1 ) {
        std::cout << "ERROR: number of cells=" << ncells_ << std::endl;
        return 1;
    }
    for(s=k=0;k<nsums_;k++){
        val = rhs_[k];
        for(l=0;l<nlist_[k];l++,s++)
            val -= data_[listcell_[s]]*listcoef_[s];
        if( fabs(val)>ZERO ){
            std::cout << "ERROR: nominal constraint=" << k << " with difference " << val << std::endl;
            return 1;
        }
    }
    for(k=0;k<ncells_;k++){
        if( abs(weight_[k])>INF/10 ){
            std::cout << "ERROR: weight of " << k << " = " << weight_[k] << std::endl;
            return 1;
        }
        if( data_[k]+ZERO<lb_[k] || data_[k]-ZERO>ub_[k] ){
            std::cout << "ERROR: data of " << k << " = " << lb_[k] << " < " << data_[k] << " < " << ub_[k] << std::endl;
            return 1;
        }
//        if( data_[k]-lpl_[k]<lb_[k]-ZERO || data_[k]+upl_[k]>ub_[k]+ZERO ){
//            printf("Warning: prot.levels of %d = %lf < %lf < %lf\n",k,data_[k]-lpl_[k],data_[k],data_[k]+upl_[k]);
//            return 1;
//        }
        if( status_[k]!='s' && status_[k]!='z' && status_[k]!='u'){
            std::cout << "ERROR: status of " << k << " = " << status_[k] << std::endl;
            return 1;
        }
    }
    std::cout << " ok" << std::endl;
    return 0;
}
