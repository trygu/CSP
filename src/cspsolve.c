/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPSOLVE.c                                       */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/
/* PWOF 10-12-2013                                  */
/* After each free() assignment to NULL added       */

#define MAX_NZ    250000     /* maximum number of non-zero elements in LP*/


#include <stdio.h>
#include <stdlib.h>
//#include <crtdbg.h>
#include <math.h>
#include "Cspdefns.h"
#include "CSPGLOB2.H"
#include "Jjsolver.h"
#include "cspback.h"
#ifdef   CHECKLP
#include "check.c"
#endif
#include "cspsolve.h"
#include "cspprice.h"
#include "cspnet.h"
#include "cspdebug.h"
#include "cspcover.h"
#include "CSPGOMO.H"
#include "cspback.h"



/*  PROTOTYPES OF FUNCTIONS */


static int    add_row(CONSTRAINT *,double *,char *,int *,double *);
static void   deletecols(int *);
static void   deleterows(int *);
static void   adapting_rhs(int *);



/*  PRIVATE GLOBAL VARIABLES */

float  toptimize = 0.0;

static int      pricing_util;       /* 0 iff pricing must be next done    */
int             pricing_done;       /* 1 iff pricing was done in this iter*/

static int      lowerb1;            /* lower bound of 1-fixed variables   */
static int      nfixed1;            /* number of variables in 'nfixed1'   */
static VARIABLE **fixed1;           /* variables with 1-fixed variables   */

static char     probname[]="SUPPRESS";
static double   *objx;
static double   *rhsx;
static char     *senx;
static int      *matbeg;
static int      *matind;
static int      *matcnt;
static double   *matval;
static double   *bdl;
static double   *bdu;


JJLPptr lp;

int load_lp()
{
        int      k;
        int      macsz;
        int      marsz;
        int      matsz;
        VARIABLE *col;

        list_pricing = NULL;
        pricing_util = 0;
        pricing_done = 0;

        mar = 1;

        macsz   = MAX_COLS_LP + 1;
        marsz   = MAX_ROWS_LP + 1;
        matsz   = MAX_NZ + 1;

        fixed1=(VARIABLE **)malloc(sizeof(VARIABLE *)*ncols);
        if(fixed1==NULL){                
                std::cout << "There is not enough memory for FIXED1" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        objx=(double *)malloc(sizeof(double)*macsz);
        if(objx==NULL){                
                std::cout << "There is not enough memory for OBJX" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((rhsx=(double *)malloc(sizeof(double)*marsz))==NULL){                
                std::cout << "There is not enough memory for RHSX" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((senx=(char *)malloc(sizeof(char)*marsz))==NULL){                
                std::cout << "There is not enough memory for SENX" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((matbeg=(int *)malloc(sizeof(int)*macsz))==NULL){                
                std::cout << "There is not enough memory for MATBEG" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((matcnt=(int *)malloc(sizeof(int)*macsz))==NULL){                
                std::cout << "There is not enough memory for MATCNT" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((matind=(int *)malloc(sizeof(int)*matsz))==NULL){                
                std::cout << "There is not enough memory for MATIND" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((matval=(double *)malloc(sizeof(double)*matsz))==NULL){                
                std::cout << "There is not enough memory for MATVAL" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((bdl=(double *)malloc(sizeof(double)*macsz))==NULL){                
                std::cout << "There is not enough memory for BDL" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((bdu=(double *)malloc(sizeof(double)*macsz))==NULL){                
                std::cout << "There is not enough memory for BDU" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }

        /*-------------------------   initial LP   -----------------------*/

        rhsx[0]       = 0.0;
        senx[0]       = 'G';
        rind[0]       = NULL;
        cind[0]       = NULL;
        objx[0]       = 1.0E+20;
        matbeg[0]     = 0;
        matcnt[0]     = 1;
        matind[0]     = 0;
        matval[0]     = 1.0;
        bdl[0]        = 0.0;
        bdu[0]        = 1.0;

        mac=1;

        lowerb1 = 0;
        nfixed1 = 0;
        for(k=0;k<ncols;k++){
            col = columns+k;
            switch( col->stat ){
                case WAITING:
                     col->val = 0.0;
                     insert_list_pricing( col );
                     break;
                case FIX_LB:
                     col->val = 0.0;
                     break;
                case FIX_UB:
                     col->val = 1.0;
                     lowerb1 += col->weight;
                     if(col->val < 1-ZERO ){                         
                         std::cout << "ERROR: in lowerb1 " << std::endl;
                         CSPexit(EXIT_ERROR); //exit(1);
                     }
                     fixed1[ nfixed1++ ] = col;
                     break;
                default:
                     col->lp = mac;
                     if( mac==MAX_COLS_LP ){                         
                         std::cout << "ERROR: mac=" << mac << " : too many LP variables" << std::endl;
                         CSPexit(EXIT_ERROR); //exit(1);
                     }
                     cind[mac]     = col;
                     objx[mac]     = (double)col->weight;
                     matbeg[mac]   = 1;
                     matcnt[mac]   = 0;

                     bdl[mac] = 0.0;
                     bdu[mac] = 1.0;
#ifdef PARTIAL
                     if(col->sensitive) bdl[mac] = 0.5;
#endif
                     mac++;
                     break;
            }
        }

#ifdef CHECKLP
        k = JJcheckprob (probname, mac, mar, 0, 1, objx, rhsx,
                       senx, matbeg, matcnt, matind, matval,
                       bdl , bdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       macsz, marsz, matsz, 0, 0, 0, 0, 0, NULL);

        if(k)CSPexit(EXIT_LPSOLVER); //exit(1);
#endif
        lp = JJloadprob (probname, mac, mar, 0, 1, objx, rhsx,
                       senx, matbeg, matcnt, matind, matval,
                       bdl , bdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       macsz, marsz, matsz, 0, 0, 0, 0, 0);
        if ( lp == NULL ){                
                std::cout << "ERROR: loading an LP (probably because no hardware-key)" << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
        }
#ifdef STAMP    
        JJsetscr_ind(lp, 0);
#endif
        return(0);
}


int unload_lp()
{
        int k;
        struct PRICE *ptr;

#ifdef  STAMP
        JJsetscr_ind(lp, 0);
#endif
        JJfreeprob(/*(void*)*/&lp);

        for(k=1;k<mar;k++){
            rind[k]->lp   = -1;
            rind[k]->stat = FIX_LB;
            rind[k]       = NULL;
        }
        for(k=1;k<mac;k++){
            cind[k]->lp   = -1;
            cind[k]->stat = FIX_LB;
            cind[k] = NULL;
        }
        mar = 0;
        mac = 0;
/*****
		JJfreedata(NULL, objx, rhsx,
                       senx, matbeg, matcnt, matind, matval,
                       bdl , bdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL);
*****/
        free((void *)objx);
        objx = NULL;
        free((void *)rhsx);
        rhsx = NULL;
        free((void *)senx);
        senx = NULL;
        free((void *)matbeg);
        matbeg = NULL;
        free((void *)matind);
        matind = NULL;
        free((void *)matcnt);
        matcnt = NULL;
        free((void *)matval);
        matval = NULL;
        free((void *)bdl);
        bdl = NULL;
        free((void *)bdu);
        bdu = NULL;

        free((void *)fixed1);
        fixed1 = NULL;

        while( list_pricing ){
            ptr = list_pricing->next;
            free( list_pricing );
            list_pricing = ptr;
        }
        return(0);
}


void put_base()
{
        int k;
        int *cstat,*rstat;

        if((cstat=(int *)calloc( (size_t)mac,sizeof(int)))==NULL){                
                std::cout << "No memory for CSTAT" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((rstat=(int *)calloc( (size_t)mar,sizeof(int)))==NULL){                
                std::cout << "No memory for RSTAT" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        rstat[0] = 0;
        cstat[0] = 1;
        for(k=1;k<mar;k++)
                rstat[k] = rind[k]->stat;
        for(k=1;k<mac;k++)
                cstat[k] = cind[k]->stat;
        if( JJloadbase(lp,cstat,rstat) ){                
                std::cout << "ERROR: The initial basis was not load." << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
        }
        free((void *)cstat);
        cstat = NULL; /*PWOF*/
        free((void *)rstat);
        rstat = NULL; /*PWOF*/
}

void get_base()
{
    int k;
    int *cstat,*rstat;

    cstat=(int *)malloc(sizeof(int)*mac);
    rstat=(int *)malloc(sizeof(int)*mar);
    if ( JJgetbase(lp,cstat,rstat) ){
        CSPexit(EXIT_LPSOLVER);
        //exit(1);
    }

    for( k=1 ; k<mac ; k++ )
        cind[k]->stat = cstat[k];
    for( k=1 ; k<mar ; k++ )
        rind[k]->stat = rstat[k];

    free((void *)cstat);
    cstat = NULL; /*PWOF*/
    free((void *)rstat);
    rstat = NULL; /*PWOF*/
}





/*   routine to introduce particular rows: queue[0],...,queue[rcnt-1]
     in the 'rows' structure            */

int add_rows(int rcnt,CONSTRAINT **queue)

{
        //int     i,j;
        int i;

        int     ccnt=0;
        int     nzcnt=0;
        double  *rhs    ={NULL};
        char    *sense  ={NULL};
        int     *rmatbeg={NULL};
        int     *rmatind={NULL};
        double  *rmatval={NULL};

        int     largemem;

        if(rcnt==0)return(0);

#ifdef STAMP
        std::cout << " ... adding " << rcnt << " rows ... " << std::endl;
#endif
        if( mar+rcnt >= MAX_ROWS_LP ){            
            std::cout << "ERROR: too many LP rows: " << mar+rcnt << std::endl;
            CSPexit(EXIT_ERROR); //exit(1);
        }
        largemem = rcnt * mac;
        if(( rhs=(double *)malloc( sizeof(double)*rcnt ) )==NULL){                
                std::cout << "There is not enough memory for vector RHS" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if(( sense=(char *)malloc( sizeof(char)*rcnt ) )==NULL){                
                std::cout << "There is not enough memory for vector SENSE" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if(( rmatbeg=(int *)malloc( sizeof(int)*rcnt ) )==NULL){                
                std::cout << "There is not enough memory for vector RMATBEG" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if(( rmatind=(int *)malloc( largemem*sizeof(int) ) )==NULL){                
                std::cout << "There is not enough memory for vector RMATIND" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if(( rmatval=(double *)malloc( largemem*sizeof(double) ) )==NULL){
                largemem = largemem/2;
                if( largemem==0 ){                     
                     CSPexit(EXIT_MEMO); //exit(1);
                }
                if(( rmatval=(double *)malloc( largemem*sizeof(double) ) )==NULL){                     
                     std::cout << "There is not enough memory for vector RMATVAL (addrow)" << std::endl;
                     CSPexit(EXIT_MEMO); //exit(1);
                }
        }

        for(i=0;i<rcnt;i++){
            rmatbeg[i] = nzcnt;            
            nzcnt += add_row(queue[i],rhs+i,sense+i,rmatind+nzcnt,rmatval+nzcnt);
            if( nzcnt>largemem ){                
                std::cout << "ERROR: mem for add rows = " << largemem << " not enought" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
            }
        }

#ifdef STAMP    
        JJsetscr_ind(lp,0);
#endif
        if (JJaddrows(lp,ccnt,rcnt,nzcnt,rhs,sense,
                rmatbeg,rmatind,rmatval,NULL,NULL) ) {                
                std::cout << " ERROR: it was not possible to add constraints" << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
        }

        free((void *)rhs);
        rhs = NULL; /*PWOF*/
        free((void *)sense);
        sense = NULL; /*PWOF*/
        free((void *)rmatbeg);
        rmatbeg = NULL; /*PWOF*/
        free((void *)rmatind);
        rmatind = NULL; /*PWOF*/
        free((void *)rmatval);
        rmatval = NULL; /*PWOF*/
#ifdef STAMP
        if( control_ind() ) CSPexit(EXIT_ERROR); //exit(1);
#endif
        return(0);
}

static int add_row(CONSTRAINT *con,double *rhs,char   *sense,int    *rmatind,double *rmatval)

{
    //int      card,j;
    int      card;
    int      l = 0;
    VARIABLE **stack;
    double   *coef;
    
    // PWOF intialising
    coef=NULL;
    stack=NULL;
    card=0;
    //

    if( con->lp != -1 ){        
        std::cout << "ERROR: constraints already in the LP (" << con->lp << ")" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    con->lp     = mar;
    *sense      = con->sense;
    *rhs        = con->rhs;
    rind[mar++] = con;

    switch( con->type ){

    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        stack = con->stack;
        coef = con->coef;
        card = con->card;        
        break;
    case COVER:
        stack = (VARIABLE **)malloc( ncols*sizeof(VARIABLE *) );
        coef = (double *)malloc( ncols*sizeof(double) );
        card = extend_cover(con,stack,coef);
        break;
    default:        
        std::cout << "ERROR: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    while(card--){
        switch( stack[card]->stat ){
            case WAITING:
                 break;
            case FIX_LB:
                 break;
            case FIX_UB:
                 *rhs -= coef[card];
                 break;
            default:
                 rmatind[l] = stack[card]->lp;                 
                 rmatval[l] = coef[card];
                 l++;
        }
    }

#ifdef STAMP
    if(l==0){
            std::cout << "WARNING: empty left-hand-side constraint" << std::endl;
            if(( *rhs>ZERO && *sense!='L' )||( *rhs<-ZERO && *sense!='G' ))
                std::cout << "        (probably infeasible LP)"  << std::endl;
            else if ((*rhs< ZERO && *rhs>-ZERO ) ||
                     (*rhs> ZERO && *sense=='L') ||
                     (*rhs<-ZERO && *sense=='G') )
                std::cout << "        (probably not useful cut)"  << std::endl;
            else{                
                std::cout << "        ERROR: bad cut"  << std::endl;
                /* print_row(con); */
                CSPexit(EXIT_ERROR); //exit(1);
            }
    }
#endif
    rmatind[l] = 0;
    rmatval[l] = *rhs;
    l++;
#ifdef STAMP    
    control_constraint( con );
#endif

    switch( con->type ){
    case COVER:
        free(stack);
        stack = NULL; /*PWOF*/
        free(coef);
        coef = NULL; /*PWOF*/
        break;
    case GOMORY:
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
        break;
    default:        
        std::cout << "ERROR: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    return(l);
}

int  add_cols(int  ccnt,VARIABLE **stack)

{
    int      i,k;
    double   c;
    VARIABLE *col;
    
    int     nzcnt;
    double  *cobj   ={NULL};
    int     *cmatbeg={NULL};
    int     *cmatind={NULL};
    double  *cmatval={NULL};
    double  *cbdl   ={NULL};
    double  *cbdu   ={NULL};

    if (ccnt==0) return(0);
#ifdef STAMP
    std::cout << " ... adding " << ccnt << " cols ..." << std::endl;
#endif
    
    if (mac+ccnt>=MAX_COLS_LP){        
        std::cout << " ERROR: too many LP columns" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    if(( cobj=(double *)malloc( sizeof(double)*ccnt ) )==NULL){        
        std::cout << "Not enough memory for COBJ" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cbdl=(double *)malloc( sizeof(double)*ccnt ) )==NULL){        
        std::cout << "Not enough memory for CBDL" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cbdu=(double *)malloc( sizeof(double)*ccnt ) )==NULL){        
        std::cout << "Not enough memory for CBDU" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cmatbeg=(int *)malloc( sizeof(int)*ccnt ) )==NULL){        
        std::cout << "Not enough memory for CMATBEG" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cmatind=(int *)malloc( ccnt*mar*sizeof(int) ) )==NULL){        
        std::cout << "Not enough memory for CMATIND" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cmatval=(double *)malloc( ccnt*mar*sizeof(double) ) )==NULL){        
        std::cout << "Not enough memory for CMATVAL" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }

    nzcnt=0;
    for (k=0;k<ccnt;k++){
        col = stack[k];
        col->stat   = LP_LB;
        col->lp     = mac;
        cind[mac++] = col;

        cobj[k]    =col->weight;
        cbdl[k]    =0.0;
        cbdu[k]    =1.0;
        cmatbeg[k] =nzcnt;
        for (i=1;i<mar;i++){
            c = get_coeficient(col,rind[i]);
            if( fabs(c)>ZERO ){
                cmatind[nzcnt] = i;
                cmatval[nzcnt] = c;
                nzcnt++;
            }
        }
    }

    if( JJaddcols(lp,ccnt,nzcnt,cobj,cmatbeg,cmatind,cmatval,cbdl,cbdu,NULL) ){
        CSPexit(EXIT_LPSOLVER);
        //exit(1);
    }

/*  control();  */

    free((char *)cobj);
    cobj = NULL; /*PWOF*/
    free((char *)cbdl);
    cbdl = NULL; /*PWOF*/
    free((char *)cbdu);
    cbdu = NULL; /*PWOF*/
    free((char *)cmatbeg);
    cmatbeg = NULL; /*PWOF*/
    free((char *)cmatind);
    cmatind = NULL; /*PWOF*/
    free((char *)cmatval);
    cmatval = NULL; /*PWOF*/
    return(ccnt);
}

/*
int  add_cols_all()
{
    int      i,j,k;
    VARIABLE *col;
    struct   PRICE *p;
    
    int     ccnt,nzcnt;
    double  *cobj   ={NULL};
    int     *cmatbeg={NULL};
    int     *cmatind={NULL};
    double  *cmatval={NULL};
    double  *cbdl   ={NULL};
    double  *cbdu   ={NULL};

    ccnt = npricing;
    if (ccnt==0) return(0);
#ifdef STAMP
    printf(" ... adding %3d cols ...\n",ccnt);
#endif
    npricing = 0;
    
    if (mac+ccnt>=MAX_COLS_LP){
        printf(" ERROR: too many LP columns\n");
        CSPexit(EXIT_ERROR); //exit(1);
    }

    if(( cobj=(double *)malloc( sizeof(double)*ccnt ) )==NULL){
        printf("No hay memoria para vector COBJ");
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cbdl=(double *)malloc( sizeof(double)*ccnt ) )==NULL){
        printf("No hay memoria para vector CBDL");
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cbdu=(double *)malloc( sizeof(double)*ccnt ) )==NULL){
        printf("No hay memoria para vector CBDU");
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cmatbeg=(int *)malloc( sizeof(int)*ccnt ) )==NULL){
        printf("No hay memoria para vector CMATBEG");
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cmatind=(int *)malloc( ccnt*mar*sizeof(int) ) )==NULL){
        printf("No hay memoria para vector CMATIND");
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if(( cmatval=(double *)malloc( ccnt*mar*sizeof(double) ) )==NULL){
        printf("No hay memoria para vector CMATVAL");
        CSPexit(EXIT_MEMO); //exit(1);
    }

    nzcnt=0;
    for (k=0;k<ccnt;k++){
#ifdef STAMP
        if(list_pricing==NULL){
            printf("ERROR: empty pricing list\n");
            CSPexit(EXIT_ERROR); //exit(1);
        }
#endif
        col = list_pricing->col;
        col->stat   = LP_LB;
        col->lp     = mac;
        cind[mac++] = col;

        cobj[k]    =col->weight;
        cbdl[k]    =0.0;
        cbdu[k]    =1.0;
        cmatbeg[k] =nzcnt;
        for (i=1;i<mar;i++){
            j = get_coeficient(col,rind[i]);
            if( j ){
                cmatind[nzcnt] = i;
                cmatval[nzcnt] = j;
                nzcnt++;
            }
        }
        p = list_pricing;
        list_pricing = list_pricing->next;
        free(p);
    }
#ifdef STAMP
    if(list_pricing){
        printf("ERROR: non empty pricing list\n");
        CSPexit(EXIT_ERROR); //exit(1);
    }
#endif
    if( JJaddcols(lp,ccnt,nzcnt,cobj,cmatbeg,cmatind,cmatval,cbdl,cbdu,NULL) )
        CSPexit(EXIT_LPSOLVER); //exit(1);

  control(); 

    free((char *)cobj);
    free((char *)cbdl);
    free((char *)cbdu);
    free((char *)cmatbeg);
    free((char *)cmatind);
    free((char *)cmatval);
    return(ccnt);
}
*/


void rem_cols(int         card,VARIABLE    **stack)

{
    int      k;
    struct   PRICE *p,*anterior,*posterior;

    anterior = (struct PRICE*) malloc(sizeof(struct PRICE)); // PWOF added memory allocation

    if (card==0) return;
#ifdef STAMP
    std::cout << " ... fixing " << card << " cols ..." << std::endl;
#endif

    p = list_pricing;
    k = 0;
    while( p ){
        posterior = p->next;
        if( stack[k]==p->col ){
            p->col->stat = FIX_LB;
            if( p==list_pricing ) list_pricing = posterior;
            else anterior->next = posterior;
            free( p );
            p = NULL; /*PWOF*/
            if( ++k == card ) break;
        } else {
            anterior = p;
        }
        p = posterior;
    }
}






double dualcost(double *u,double *dj)

{
    double lb;
    if ( JJoptimize(lp) ) {
        CSPexit(EXIT_LPSOLVER);
        //exit(1);
    }
    if ( JJsolution (lp, NULL, &lb, NULL, u, NULL, dj) ) {
        CSPexit(EXIT_LPSOLVER);
        //exit(1);
    }
    return lb;
}


double solve_lp(int cleaning)

{
    double lb;
    float  t1;
    
    t1 = seconds();
    if ( JJdualopt(lp) ){        
        std::cout << " it was not possible to solve with DUALOPT " << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    toptimize += seconds()-t1;
    if( pricing_util%5==0 ){
PRIC:
        while( pricing(cleaning) ){
            pricing_util = -1;
            t1 = seconds();
            if ( JJoptimize(lp) ){                
                std::cout << " it was not possible to solve with OPTIMIZE " << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
            }
            toptimize += seconds()-t1;
        }
        pricing_done = 1;
    } else {
        pricing_done = 0;
    }
    pricing_util ++;


#ifdef STAMP
    if( JJgetitc(lp)==0 ){
        std::cout << "WARNING: no LP iterations!" << std::endl;
    }
#endif


    switch( JJgetstat(lp) ){
        case JJ_OPTIMAL:
        case JJ_OPTIMAL_INFEAS:
            JJgetobjval( lp , &lb );
            lb += lowerb1;
            if( ceil(lb-ZERO) > upperb-ZERO && pricing_done==0 )
                goto PRIC;            
            return(lb);
        case JJ_INFEASIBLE:
            if( pricing_done==0 )
                goto PRIC;
#ifdef STAMP
            std::cout << " WARNING: non-feasible LP problem" << std::endl;
#endif
            return( (double)upperb );
        case JJ_UNBOUNDED:
            if( pricing_done==0 )
                goto PRIC;
#ifdef STAMP
            std::cout << " WARNING: non-bounded LP problem" << std::endl;
#endif
            return( (double)upperb );
        default:            
            std::cout << " ERROR: pstat=" << JJgetstat(lp)  << std::endl;
            JJlpwrite(lp,fsdclp);
            CSPexit(EXIT_ERROR); //exit(1);
    }
    return(0.0);
}

void get_solution()
{
    double *xval;
    int    i;
#ifdef STAMP    
    int    k;
#endif

    i = JJgetstat(lp);
    switch( i ){
        case JJ_OPTIMAL:
        case JJ_OPTIMAL_INFEAS:
            JJgetobjval( lp , &lowerb );
            lowerb += lowerb1;
            break;
        case JJ_INFEASIBLE:
#ifdef STAMP
            std::cout << " WARNING: non-feasible LP problem" << std::endl;
#endif
            lowerb = (double)upperb ;
            break;
        case JJ_UNBOUNDED:
            std::cout << " WARNING: non-bounded LP problem" << std::endl;
            lowerb = (double)upperb ;
            break;
        default:            
            std::cout << " ERROR: unknown LP problem status = " << i << std::endl;
            CSPexit(EXIT_LPSOLVER); //exit(1);
    }

#ifdef STAMP
    std::cout << " ... LP: UB=" << upperb << " LB=" << (float)lowerb << " tit=" << JJgetitc(lp) << " Iit=" << JJgetitci(lp) << " mac=" << JJgetmac(lp) << " mar=" << JJgetmar(lp) << " nz=" << JJgetmat(lp);
#endif

    for( nsupport=0 ; nsupport<nfixed1 ; nsupport++ )
        support[ nsupport ] = fixed1[ nsupport ];

    xval=(double *)malloc( mac*sizeof(double) );
    JJgetx(lp, xval, 0, mac-1);
    for(i=1;i<mac;i++){
        cind[i]->val = xval[i];
        if( xval[i]>ZERO ) support[ nsupport++ ] = cind[i];
    }
    
#ifdef STAMP
    for(k=0,i=1;i<mac;i++)
        if( ZERO<xval[i] && xval[i]<1-ZERO ){
            k++;
        }
    std::cout << " frac=" << k << std::endl;
#endif
    free((void *)xval);
    xval = NULL; /*PWOF*/
}




/*********
int reduced_cost()
{
    int    k,lb;
    double *dj;

    dj     = (double *)malloc(mac*sizeof(double));
    if( JJgetdj(lp,dj,0,mac-1) ){
        puts(" No reduced cost avalaible ");
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    for(k=1;k<mac;k++){
        lb = (int)ceil(lowerb+dj[k]);
        if( lb > cind[k]->lb )
            cind[k]->lb = lb;
    }
    free((void *)dj);
    return(0);
}
********/

void del_col(VARIABLE *var,int sta)

{
    int *status;
    status = (int *)calloc( (size_t)mac,sizeof(int));
    if(var->lp<1 || var->lp>=mac){        
        std::cout << "ERROR: var->lp=" << var->lp << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    status[var->lp] = sta;
    deletecols(status);    
    free(status);
    status = NULL; /*PWOF*/
}


int del_cols()
{
    int    i,num0,num1,num2;
    double rc;
    int    *status;
    double *dj;
#ifdef STAMP
	double obj;
#endif
    if( pricing_done==0 )
        return(0);

    dj     = (double *)malloc( (size_t)mac * sizeof(double) );
    if( JJgetdj( lp, dj, 0, mac-1) ){        
        std::cout << "ERROR: not reduced-costs" << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }

    status = (int *)calloc( (size_t)mac,sizeof(int));
    num0 = num1 = num2 = 0;
    for(i=1;i<mac;i++){
        rc = dj[i];
        if( rc > ZERO ){
            if( ceil(lowerb+rc-ZERO)+ZERO > upperb ){
                status[i]=FIX_LB;
                num0++;
            } else if( list_pricing && ceil(lowerb + 2 * rc-ZERO)+ZERO > upperb ){
                status[i]=WAITING;
                insert_list_pricing( cind[i] );
                num2++;
            }
        } else if( rc < -ZERO ){
            if( ceil(lowerb-rc-ZERO)+ZERO > upperb ){
                /* puts(" could be fixed to UB "); */
                status[i]=FIX_UB; 
                num1++; 
            }
        }
    }
#ifdef STAMP
    std::cout << " ... deleted cols ... " << num0 << " to 0, " << num1 << " to 1, " << num2 << " sleep" << std::endl;
#endif
    free((void *)dj);
    dj = NULL; /*PWOF*/

    if(num0+num1+num2==0){
        free((void *)status);
        status = NULL; /*PWOF*/
        return 0;
    }
    
    deletecols(status);
    free((void *)status);    
    status = NULL; /*PWOF*/
    if ( JJoptimize(lp) ){        
        std::cout << " it was not possible to solve with OPTIMIZE " << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
#ifdef STAMP
    if( JJgetitc(lp) ){        
        JJgetobjval( lp , &obj );
        std::cout << "WARNING: LP iterations! (new obj=" << lowerb1+obj << ")" << std::endl;
        CSPexit(EXIT_LPSOLVER);
/**
        exit(1);
**/        
    }
#endif    
    return(num0+num1+num2);
}

static void deletecols(int *status)

{
    int i,j;

    adapting_rhs(status);
    j=i=1;
    while(1){
        while(i<mac && status[i]){
            cind[i]->stat = status[i];
            cind[i]->lp   = branchs;  /*** to know when at root node ***/
            cind[i]->val  = ( status[i]==FIX_UB ? 1 : 0 );
            status[i] = 1;
            i++;
        }
        if(i==mac)break;
        ( cind[j]=cind[i] )->lp = j;
        j++;
        i++;
    }
    mac = j;
    JJdelsetcols(lp,status);    
}

static void adapting_rhs(int *status)

{
    int i,j;
    int nzcnt,surplus;
    int cmatbeg[2];
    int *cmatind;
    double *rhs,*newrhs,*cmatval;

    for(i=1,j=0;i<mac;i++)
        if( status[i]==FIX_UB ) j++;
    if(j==0) return;

    cmatind = (int *)malloc( mar * sizeof(int) );
    cmatval = (double *)malloc( mar * sizeof(double) );
    rhs     = (double *)malloc( mar * sizeof(double) );
    newrhs  = (double *)malloc( mar * sizeof(double) );
    JJgetrhs(lp,rhs,0,mar-1);
    for(i=0;i<mar;i++) newrhs[i]=rhs[i];
    for(i=1;i<mac;i++)
        if( status[i]==FIX_UB ){
            JJgetcols(lp,&nzcnt,cmatbeg,cmatind,cmatval,mar,&surplus,i,i);
            for(j=0;j<nzcnt;j++)
                newrhs[ cmatind[j] ] -= cmatval[j];
            lowerb1 += cind[i]->weight;
            fixed1[ nfixed1++ ] = cind[i];
        }
    for(i=0;i<mar;i++)
        if( newrhs[i]<rhs[i]-ZERO || newrhs[i]>rhs[i]+ZERO )
            JJchgcoef(lp,i,-1,newrhs[i]);
    free( cmatind );
    cmatind = NULL; /*PWOF*/
    free( cmatval );
    cmatval = NULL; /*PWOF*/
    free( rhs );
    rhs = NULL; /*PWOF*/
    free( newrhs );
    newrhs = NULL; /*PWOF*/
}



int del_rows()
{
    int    i,num0,num1;
    int    *status;
    double *slack;

    if( pricing_done==0 ) return(0);

    slack  = (double *)malloc( (size_t)mar * sizeof(double) );
    if( JJgetslack( lp, slack, 0, mar-1) ){        
        std::cout << "ERROR: not slacks" << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    status = (int *)calloc( (size_t)mar,sizeof(int));
    for(num0=num1=0,i=1;i<mar;i++)
        if( rind[i]->type != BRANCHCUT )
        switch( rind[i]->sense ){
            case 'G':
                if( slack[i] < - MAX_SLACK){
                    status[i]=1;
                    num0++;
                }
                break;
            case 'L':
                if( slack[i] >   MAX_SLACK){
                    status[i]=1;
                    num1++;
                }
                break;
            default:
                break;
        }
#ifdef STAMP
    std::cout << " ... deleted rows ... " << num0 << " + " << num1 << std::endl;
#endif
    free( (void *)slack);
    slack = NULL; /*PWOF*/
    
    if(num0+num1==0){
        free((void *)status);
        status = NULL; /*PWOF*/
        return 0;
    }

    deleterows(status);
    free((void *)status);
    status = NULL; /*PWOF*/
    if ( JJoptimize(lp) ){        
        std::cout << " it was not possible to solve with OPTIMIZE " << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
#ifdef STAMP
    if( JJgetitc(lp) ){        
        std::cout << "WARNING: LP iterations!" << std::endl;
        CSPexit(EXIT_LPSOLVER);
/**        
        exit(1);
**/
    }
#endif
    return(num0+num1);
}

static void deleterows(int *status)

{
    int i,j;

    j=i=1;
    while(1){
        while(i<mar && status[i]){
            rind[i]->stat = FIX_LB;
            rind[i]->lp   = -1;
            i++;
        }
        if(i==mar)break;
        ( rind[j]=rind[i] )->lp = j;
        j++;
        i++;
    }
    mar = j;
    JJdelsetrows(lp,status);
}

void deletelastrow()
{
    int *status;
#ifdef STAMP
    std::cout << " ... deleting last rows\n";
#endif
    status = (int *)calloc( mar , sizeof(int) );
    mar--;
    status[mar] = 1;
    JJdelsetrows(lp,status);
    free(status);
    status = NULL; /*PWOF*/
    rind[ mar ]->lp = -1;
    rind[ mar ] = NULL;
}


double     solve_child_col(VARIABLE   *col,double val,int sol)

{
    double lb;
    int    index[2];
    double value[2];
    char   lu[2];

    if(col->lp<0){        
        std::cout << "ERROR: variable not in the LP: stat=" << col->stat << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    index[0] = col->lp;
    lu[0]    = 'B';
    value[0] = val;
    if( JJchgbds(lp,1,index,lu,value) ){        
        std::cout << "ERROR: it is not possible to fix LP-var " << col->index << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }

    pricing_util = 0;
    lb = solve_lp(0);
    if( sol ){
        get_solution();
        get_base();
    }
    integer_solution();

    index[1] = col->lp;
    lu[0]    = 'L';
    lu[1]    = 'U';
    value[0] = 0.0;
    value[1] = 1.0;
    if( JJchgbds(lp,2,index,lu,value) ){        
        std::cout << "ERROR: it is not possible to fix LP-var " << col->index << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    return(lb);
}

double     solve_child_row(CONSTRAINT *row,int sol)

{
    int    i,rmatbeg,cont;
    double lb,rhs;
    int    *rmatind;
    double *rmatval;
    char   sense;
    
    cont    = row->card;
    rhs     = row->rhs;
    sense   = row->sense;
    rmatbeg = 0;
    rmatind = (int *)malloc( cont * sizeof(int) );
    if( rmatind==NULL ){        
        std::cout << " ERROR: not enough memory for branch-constraint" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(i=0;i<cont;i++) rmatind[i] = row->stack[i]->lp;
    rmatval = (double *)malloc( cont * sizeof(double) );
    if( rmatval==NULL ){        
        std::cout << " ERROR: not enough memory for branch-constraint" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(i=0;i<cont;i++) rmatval[i] = 1.0;

    if (JJaddrows(lp,0,1,cont,&rhs,&sense,
        &rmatbeg,rmatind,rmatval,NULL,NULL) ) {        
        std::cout << " ERROR: it was not possible to add branch constraint" << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    free( (void *)rmatval );
    rmatval = NULL; /*PWOF*/
    free( (void *)rmatind );
    rmatind = NULL; /*PWOF*/

    pricing_util = 0;
    lb = solve_lp(0);
    if( sol ) {
        get_solution();
        get_base();
    }
    integer_solution();

    if( JJdelrows(lp,mar,mar) ){        
        std::cout << " ERROR: it was not possible to delete branch constraint" << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    return(lb);
}



int integer_solution(void)
{
    int      i,ub,nsup;
    double   *xval;
    double   obj;
    VARIABLE **sup;

    JJoptimize(lp);
    if( JJgetstat(lp)!=JJ_OPTIMAL && JJgetstat(lp) !=JJ_OPTIMAL_INFEAS ){
#ifdef STAMP        
        std::cout << "WARNING: LP status = " << JJgetstat(lp)  << std::endl;
#endif
        return(1);
    }
    xval = (double *)malloc( mac*sizeof(double) );
    JJgetx(lp, xval, 0, mac-1);
    i = mac;
    while(--i)
        if( xval[i]>ZERO && xval[i]<1-ZERO ) {
            free(xval);
            xval = NULL; /*PWOF*/
            return(0);
        }
    sup = (VARIABLE **)malloc( ncols*sizeof(VARIABLE *) );
    for(nsup=0;nsup<nfixed1;nsup++) sup[ nsup ] = fixed1[ nsup ];
    for(i=1;i<mac;i++) if( xval[i] > ZERO ) sup[ nsup++ ] = cind[i];
    if( protected_flow(nsup,sup) ) {
        JJgetobjval( lp , &obj );
        ub = lowerb1 + (int)ceil(obj-ZERO);
        if( upperb > ub ){
#ifdef STAMP
                std::cout << " Integer solution of value " << ub << std::endl;
                std::cout << "    BETTER THAN THE CURRENT ONE" << std::endl;
#endif
                ubtype  = 'L';
                topti   = seconds()-t0;
                upperb  = ub;
                nbetter = nsup;
                for(i=0;i<nsup;i++)
                    better[i] = sup[i];
#ifdef STAMP                
                write_heu(fheuristi);
#endif
                CSPnewsolution();
        }
        free( xval );
        xval = NULL; /*PWOF*/
        free( sup );
        sup = NULL; /*PWOF*/
        return(1);
    }
    free( xval );
    xval = NULL; /*PWOF*/
    free( sup );
    sup = NULL; /*PWOF*/
    return(0);
}





void setup_lp()
{
#ifdef STAMP
    double obj;
#endif
    JJoptimize(lp);
#ifdef STAMP
    if( JJgetitc(lp) ){
        JJgetobjval( lp , &obj );
        std::cout << "WARNING: " << JJgetitc(lp) << " LP iterations in setup_lp! LB=" << lowerb1+obj << std::endl;
    }
#endif
}



void control()
{
   int i,j;

   if( mar != JJgetmar(lp) || mac != JJgetmac(lp) ){       
       std::cout << "ERROR: mar=" << mar << " mac=" << mac << " getmar=" << JJgetmar(lp) << " getmac=" << JJgetmac(lp) << std::endl;
       CSPexit(EXIT_ERROR); //exit(1);
   }
   for(i=0;i<ncols;i++){
       j=columns[i].lp;
       if(j>0 && cind[j]->index!=i){           
           std::cout << "\n ERROR in control: cind[" << j << "]=" << cind[j]->index << " stat[" << i << "]=" << columns[i].lp << std::endl;
           CSPexit(EXIT_ERROR);
       }
   }
}


void activa_pricing()
{
    pricing_util = 0;
}

