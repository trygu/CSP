/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPHEUR.c                                        */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/

/*
#define CHECKLP
*/

#include <stdio.h>
#include <stdlib.h>
#include "Cspdefns.h"
#include "CSPGLOB2.H"
#include "Jjsolver.h"
#ifdef CHECKLP
#include "check.c"
#endif
#include "cspheur.h"
#include "cspsep.h"
#include "cspdebug.h"
#include "cspback.h"



static JJLPptr Nlp; // included by Salome 12/2/2012

#define   MAXLPITER 1000
#define   EPSILON   1.0
#define   BIGVALUE  1.0e+13  
                      //2140000000.0           //1.0e+11

/* PROTOTYPES OF FUNCTIONS */

static void   load_2network(char*);
static void   unload_2network(void);
static double protection_2level(VARIABLE *,int,double,int);
static void   update_heuristic(int,int,VARIABLE **);
static void   remove_cell(int,double*);
static void   restore_cell(int,double*);
static void   potential_innecesary(int*,int*);


/*================================================================*/

int heuristic(int improv)

{
    int        k,l,nstack,nsup,weight;
    PROT_LEVEL *pro;
    char       *status;
    int        *stack;
    double     bd[4];
    VARIABLE   **sup;
    float      t1,t2;
    double     *oldval;
    double     obj;

#ifdef STAMP
    std::cout << "    .. heuristic .." << std::endl;
#endif
    t1 = seconds();



    status = (char *)calloc( Rncells, sizeof(char) );
    if(status==NULL){        
        std::cout << "There is not enough memory for STATUS" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    oldval = (double *)malloc( Rncells * sizeof(double) );
    if( oldval==NULL ){        
        std::cout << "There is not enough memory for OLDVAL" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(l=0;l<Rncells;l++)
        oldval[l] = columns[l].val;

    /*** maximum set of variables for the heuristic ***/
    for(l=0;l<Rncells;l++)
        if( oldval[l]>ZERO )
            status[l]=1;
    
    for(l=0;l<nbetter;l++)
	status[ better[l]->index ]=1;
    load_2network(status);
    free( status );
    status = NULL; /*PWOF*/


    /*** initial set of variables for the heuristic  ***/
    for(l=0;l<Rncells;l++){
        if( columns[l].val > EPSILON-ZERO ) columns[l].val = 1.0;
        else                                columns[l].val = 0.0;
    }


/*****       constructive part    *****/


    if( Rncells>500 && Rnsums>500 && nprot_level>500 )
        JJlpiterlimit(MAXLPITER);
    

#ifdef STAMP        
        std::cout << " " << nprot_level << ":";
#endif
   
    for(k=0;k<nprot_level;k++){
        pro = prot_level + k;
#ifdef STAMP        
        if(k%1000==0) std::cout << " " << k;
#endif        
        obj = protection_2level(pro->sen->var,pro->sense,pro->level,1);  // k==2 in problem subtable      

        if( obj+ZERO>=BIGVALUE ){
#ifdef STAMP        
            std::cout << "WARNING:   not enough cells in the initial maximal set, or infeasible" << std::endl;
#endif
            weight = INF;
	    nsup   = -1;
            goto SALIR;
        }
#ifdef STAMP        
//        if( obj>1-ZERO ){
//            printf("WARNING: not enough cells in the initial potential set\n");
//        }
#endif
    }
#ifdef STAMP        
    std::cout << "\n:";
#endif
    if( Rncells>500 && Rnsums>500 && nprot_level>500 )
        JJlpiterlimit(1000000);

    
/*****       saving    *****/

    sup = (VARIABLE **)malloc( Rncells * sizeof(VARIABLE *) );
    if(sup==NULL){            
            std::cout << "There is not enough memory for SUP in heuristic" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
    }
    weight = 0;
    nsup   = 0;
    for(k=0;k<Rncells;k++)
        if( columns[k].val>1-ZERO ){
            sup[ nsup++ ] = columns+k;
            weight += columns[k].weight;
        }
    update_heuristic(weight,nsup,sup);
    free( (void*)sup );
    //sup = NULL; /*PWOF*/

    if( !improv ) goto SALIR;

#ifdef STAMP
    t2 = seconds();
    std::cout << "HEUR: weight=" << weight << " sup=" << nsup << " time = " << t2-t1 << std::endl;
#endif


/*****       clean-up part   *****/

    stack = (int *)malloc( Rncells * sizeof(int) );
    if( stack==NULL ){        
        std::cout << "There is not enough memory for STACK" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }

	
	
	/* this condition must be strengthened */
    potential_innecesary(&nstack,stack);

//#ifdef STAMP
//	printf("%d variables to be studied\n",nstack);
//#endif


    for(l=0;l<nstack;l++){
//#ifdef STAMP
//    	printf(" studying variable = %d\n",stack[l]);
//#endif
        remove_cell(stack[l],bd);
        for(k=0;k<nrows;k++)
            if( violated(rows[k]) > ZERO )
                break;
        if( k<nrows ){

//#ifdef STAMP
//            printf("CLEAN-UP: variable %d needed because POOL\n",stack[l]);
//#endif

            restore_cell(stack[l],bd);
        }else{
            for(k=0;k<nprot_level;k++){
                pro = prot_level + k;                
                if( protection_2level(pro->sen->var,pro->sense,pro->level,0)>ZERO )
                    break;
            }
            if( k<nprot_level ){
//#ifdef STAMP
//                printf("CLEAN-UP: variable %d needed because Network-LP\n",stack[l]);
//#endif
                restore_cell(stack[l],bd);
            }
        }
    }
    free(stack);
    stack = NULL; /*PWOF*/
    

/*****       saving    *****/

    sup = (VARIABLE **)malloc( Rncells * sizeof(VARIABLE *) );
    if(sup==NULL){            
            std::cout << "There is not enough memory for SUP in heuristic" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
    }
    weight = 0;
    nsup   = 0;
    for(k=0;k<Rncells;k++)
        if( columns[k].val>1-ZERO ){
            sup[ nsup++ ] = columns+k;
            weight += columns[k].weight;
        }
    update_heuristic(weight,nsup,sup);
    free( (void*)sup );
    //sup = NULL; /*PWOF*/

SALIR:

    unload_2network();
    for(k=0;k<Rncells;k++)
        columns[k].val = oldval[k];
    free( (void*)oldval );
    oldval = NULL; /*PWOF*/

    t2 = seconds();
    theur += t2-t1;
#ifdef STAMP
    std::cout << "HEUR: weight=" << weight << " sup=" << nsup << " time = " << t2-t1 << "  (nrows=" << nrows << ")" << std::endl;
#endif
    bad_heuristic = 0;
    return(weight);
}


static void   update_heuristic(int weight,int nsup,VARIABLE      **sup)

{
#ifdef STAMP
    int k,l;
#endif

    if( weight >= upperb ) return;

    ubtype      = 'H';
    upperb      = weight;
    topti       = seconds()-t0;
    nbetter     = nsup;
    while(nsup--) better[nsup] = sup[nsup];
#ifdef STAMP    
    write_heu(fheuristi);

    for(k=0;k<Rncells;k++)
        if(columns[k].sensitive){
            for(l=0;l<nbetter;l++)
                if( better[l]==columns+k) break;
            if(l==nbetter){
                std::cout << "ERROR: sensitive cell " << k << " not in heuristic" << std::endl;
                CSPexit(EXIT_ERROR); //exit(1);
            }
        }
#endif
    CSPnewsolution();
}









/*====================================================================*/



/* PRIVATE GLOBAL DATA TO INTERPRETATE THE DUAL SOLUTION */

static int    *net2cell;   /* cell number in the network of a node         */
static int    *cell2net;   /* node associated with a cell of the network   */
static int    *net2sum;    /* sum number in the network of a node          */
static int    *sum2net;    /* node associated with a sum of the network    */

/* PRIVATE GLOBAL DATA FOR LP */

static char     Nprobname[]="HNETWORK";
static double   *Nobjx;
static double   *Nrhsx;
static char     *Nsenx;
static int      *Nmatbeg;
static int      *Nmatind;
static int      *Nmatcnt;
static double   *Nmatval;
static double   *Nbdl;
static double   *Nbdu;


/**
***   Set on the network.
**/


static void   load_2network(char *status)

{
        int      i,l;
        int      mar,mac,macsz,marsz,matsz;
        VARIABLE *var;
        struct   CELDA *c;

/*
        printf("\nNetwork with NR=%d and NC=%d cols ;",NR,NC);
*/
/* allocation memory */

        macsz = 2*Rncells;
        marsz = Rnsums;
        matsz = 8*macsz;

        net2cell = (int *)malloc( Rncells * sizeof(int) );
        cell2net = (int *)malloc( Rncells * sizeof(int) );
        sum2net = (int *)malloc( Rnsums * sizeof(int) );
        net2sum = (int *)malloc( Rnsums * sizeof(int) );

        Nobjx=(double *)malloc(sizeof(double)*macsz);
        if(Nobjx==NULL){                
                std::cout << "There is not enough memory for OBJX" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nrhsx=(double *)malloc(sizeof(double)*marsz))==NULL){                
                std::cout << "There is not enough memory for RHSX" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nsenx=(char *)malloc(sizeof(char)*marsz))==NULL){                
                std::cout << "There is not enough memory for SENX" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nmatbeg=(int *)malloc(sizeof(int)*macsz))==NULL){                
                std::cout << "There is not enough memory for MATBEG" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nmatcnt=(int *)malloc(sizeof(int)*macsz))==NULL){                
                std::cout << "There is not enough memory for MATCNT" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nmatind=(int *)malloc(sizeof(int)*matsz))==NULL){                
                std::cout << "There is not enough memory for MATIND" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nmatval=(double *)malloc(sizeof(double)*matsz))==NULL){                
                std::cout << "There is not enough memory for MATVAL" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nbdl=(double *)malloc(sizeof(double)*macsz))==NULL){                
                std::cout << "There is not enough memory for BDL" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }
        if((Nbdu=(double *)malloc(sizeof(double)*macsz))==NULL){                
                std::cout << "There is not enough memory for BDU" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
        }

/* LP row construction */

        for(i=0;i<Rnsums;i++) sum2net[i]=0;
        for(i=0;i<Rncells;i++)
            if( status[i] )
                for( c=columns[i].next ; c ; c=c->next )
                    sum2net[ c->index ]++;

        mar=0;
        for(i=0;i<Rnsums;i++)
            if( sum2net[i] ){
                sum2net[i]   = mar;
                net2sum[mar] = i;
                Nrhsx[mar]   = 0;
                Nsenx[mar]   = 'E';
                mar++;
            } else
                net2sum[i]   = -1;

/* LP column construction */

        mac = 0;
        l   = 0;
        for(i=0;i<Rncells;i++){
            var = columns+i;
            if( status[i] ){
                net2cell[mac/2]  = i;
                cell2net[i]      = mac/2;
                Nobjx[mac]       = (var->val>EPSILON-ZERO ? 0.0 : (1-var->val) * var->weight +1);
//                Nobjx[mac]       = (var->val>EPSILON-ZERO ? 0.0 : var->weight +1);
                Nmatbeg[mac]     = l;
                Nmatcnt[mac]     = 0;
                for( c=var->next ; c ; c=c->next ){
                    Nmatcnt[mac]++;
                    Nmatval[l]     = c->coef;
                    Nmatind[l++]   = sum2net[ c->index ];
                    if( l==matsz ){                        
                        std::cout << "ERROR 999: not enough space " << std::endl;
                        CSPexit(EXIT_ERROR); //exit(1);
                    }
                }
                Nbdl[mac]      = 0;
                Nbdu[mac]      = var->uvalue;
                mac++;
                Nobjx[mac]     = (var->val>EPSILON-ZERO ? 0.0 : (1-var->val) * var->weight +1);
//                Nobjx[mac]     = (var->val>EPSILON-ZERO ? 0.0 : var->weight +1);
                Nmatbeg[mac]   = l;
                Nmatcnt[mac]   = 0;
                for( c=var->next ; c ; c=c->next ){
                    Nmatcnt[mac]++;
                    Nmatval[l]     = -c->coef;
                    Nmatind[l++]   = sum2net[ c->index ];
                    if( l==matsz ){                        
                        std::cout << "ERROR 888: not enough space " << std::endl;
                        CSPexit(EXIT_ERROR); //exit(1);
                    }
                }
                Nbdl[mac]      = 0;
                Nbdu[mac]      = var->lvalue;
                mac++;
            } else
                cell2net[i]    = -1;
        }


#ifdef STAMP
        std::cout << "HEUR: mac=" << mac << "  mar=" << mar << "  (nrows=" << nrows << ")" << std::endl;
#endif


#ifdef CHECKLP
        l = JJcheckprob (Nprobname, mac, mar, 0, 1, Nobjx, Nrhsx,
                       Nsenx, Nmatbeg, Nmatcnt, Nmatind, Nmatval,
                       Nbdl , Nbdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       macsz, marsz, matsz, 0, 0, 0, 0, 0, NULL);
        if(l)CSPexit(EXIT_ERROR); //exit(1);
#endif

        Nlp = JJloadprob (Nprobname, mac, mar, 0, 1, Nobjx, Nrhsx,
                       Nsenx, Nmatbeg, Nmatcnt, Nmatind, Nmatval,
                       Nbdl , Nbdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       macsz, marsz, matsz, 0, 0, 0, 0, 0);

        if ( Nlp == NULL ){
                std::cout << "ERROR: There is not enought memory for the first LP" << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
        }
        /* The protection-level network LP runs at the legacy CPLEX
           feasibility tolerance, not the 1e-4 that JJloadprob sets
           for the data-valued models. Its rows are zero-rhs balance
           equations, so rounding can never make them inconsistent,
           and the loose tolerance is actively harmful: a forced
           protection flow below 1e-4 is accepted as a plain
           violation of those rows, the LP reports optimal with
           objective 0, and no flow is rerouted through the other
           cells, so no secondary suppression is left to mark (the
           MOD D=0 regression). */
        if (JJsetoptd(Nlp, "primal_feasibility_tolerance", 1.0e-6) != 0){
                std::cout << "ERROR: could not set the network LP feasibility tolerance" << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
        }
#ifdef STAMP    
        JJsetscr_ind(Nlp, 0);
#endif

#ifdef CHECKLP
        JJmpswrite(Nlp,fmpsnet);
#endif
////        JJlpwrite(Nlp,"sdcnet.lp");
}

/**
***   Set off the network.
**/

static void unload_2network()
{
#ifdef STAMP
        JJsetscr_ind(Nlp, 0);
#endif
        JJfreeprob(/*(void*)*/&Nlp);
/*****
		JJfreedata(NULL, Nobjx, Nrhsx,
                       Nsenx, Nmatbeg, Nmatcnt, Nmatind, Nmatval,
                       Nbdl , Nbdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL);
*****/
        free((void *)Nobjx);
        Nobjx = NULL; /*PWOF*/
        free((void *)Nrhsx);
        Nrhsx = NULL; /*PWOF*/
        free((void *)Nsenx);
        Nsenx = NULL; /*PWOF*/
        free((void *)Nmatbeg);
        Nmatbeg = NULL; /*PWOF*/
        free((void *)Nmatind);
        Nmatind = NULL; /*PWOF*/
        free((void *)Nmatcnt);
        Nmatcnt = NULL; /*PWOF*/
        free((void *)Nmatval);
        Nmatval = NULL; /*PWOF*/
        free((void *)Nbdl);
        Nbdl = NULL; /*PWOF*/
        free((void *)Nbdu);
        Nbdu = NULL; /*PWOF*/

        free((void *)net2cell);
        net2cell = NULL; /*PWOF*/
        free((void *)cell2net);
        cell2net = NULL; /*PWOF*/
        free((void *)net2sum);
        net2sum = NULL; /*PWOF*/
        free((void *)sum2net);
        sum2net = NULL; /*PWOF*/

}


/**
***  Compute the minimum ('type'=1) or maximum ('type'=-1)
***  protection level 'opt' of a suppressed cell.
***  The primal solution is adapted if 'mode'==1.
**/

static double    protection_2level(VARIABLE  *var,int type,double goal,int mode)

{
    int       k,l,mac;
    //int       netstatus,netnodes,netarcs,netiter,lpiter;
    int       netstatus,lpiter;
    double    opt;
    double    *x;
    int       index1,index2;
    char      lu;
    double    bd;

    index1 = 2* cell2net[ var->index ];
    index2 = index1+1;
    if( type < 0 ){
        index1 ++;
        index2 --;
    }
    lu = 'L';         // Modified on september 2001, to allow protection leves > external bounds
    bd = goal;
    JJchgbds(Nlp,1,&index1,&lu,&bd);
//    lu = 'U';
//    bd = goal;
//    JJchgbds(Nlp,1,&index1,&lu,&bd);
//    lu = 'L';
//    bd = 0.0;
//    JJchgbds(Nlp,1,&index2,&lu,&bd);
    lu = 'U';
    bd = 0.0;
    JJchgbds(Nlp,1,&index2,&lu,&bd);

/*
    if ( JJnetopt(Nlp,&netstatus,&netnodes,&netarcs,&netiter) ){
		if( netstatus==JJ_INForUNB ){
			opt =INF;
			goto FUERA;
		}			

        puts(" it was not possible to solve with NETOPT ");
        JJlpwrite(Nlp,fsdcnetlp);
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }


    if( netstatus!=CPX_NETOPTIMAL || netnodes!=JJgetmar(Nlp) || netarcs!=JJgetmac(Nlp) ){
        printf(" WARNING in sdcnet.c - netopt\n");
        printf(" netstatus=%d netnodes=%d netarcs=%d netiter=%d mar=%d mac=%d\n",
             netstatus,netnodes,netarcs,netiter,JJgetmar(Nlp),JJgetmac(Nlp));
    }

    printf(" netstatus=%d netnodes=%d netarcs=%d netiter=%d mar=%d mac=%d\n",
             netstatus,netnodes,netarcs,netiter,JJgetmar(Nlp),JJgetmac(Nlp));
    printf(" ... solving LP:stat=%d opt=%f tit=%d Iit=%d mac=%d mar=%d nz=%d\n",
JJgetstat(Nlp),(float)opt,JJgetitc(Nlp),JJgetitci(Nlp),JJgetmac(Nlp),JJgetmar(Nlp),JJgetmat(Nlp));
*/
    /* The legacy CPLEX backend reported definitive infeasible/
       unbounded models through the return code; the HiGHS wrapper
       returns 0 for any definitive status and reports it through
       JJgetstat. Both paths must take the same branch: drop the
       level (opt = INF, below BIGVALUE) and keep the current
       suppression set. */
    if ( JJoptimize(Nlp) ){
        netstatus = JJgetstat(Nlp);
        //puts(" CELL=%s  type=%d  prot=%f (lpstat=%d)",var->name,type,goal,netstatus);
        std::cout << " CELL=" << var->name << "  type=" << type << "  prot= " << goal << " (lpstat=" << netstatus;
		if( netstatus==JJ_INForUNB || netstatus==JJ_INFEASIBLE || netstatus==JJ_UNBOUNDED ){
            opt =INF;
			goto FUERA;
		}

        //puts(" it was not possible to solve with -optimize- (lpstat=%d)",netstatus);
        std::cout << " it was not possible to solve with -optimize- (lpstat=" << netstatus << ")";
        JJlpwrite(Nlp,fsdcnetlp);
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }

//    printf(" ... solving LP:stat=%d opt=%f tit=%d Iit=%d mac=%d mar=%d nz=%d\n",
//JJgetstat(Nlp),(float)opt,JJgetitc(Nlp),JJgetitci(Nlp),JJgetmac(Nlp),JJgetmar(Nlp),JJgetmat(Nlp));

//    printf(" ... solving LP:stat=%d opt=%f tit=%d Iit=%d mac=%d mar=%d nz=%d\n",
//JJgetstat(Nlp),(float)opt,JJgetitc(Nlp),JJgetitci(Nlp),JJgetmac(Nlp),JJgetmar(Nlp),JJgetmat(Nlp));
    mac = JJgetmac(Nlp);
    netstatus = JJgetstat(Nlp);
    lpiter = JJgetitc(Nlp);

    if( netstatus==JJ_INForUNB || netstatus==JJ_INFEASIBLE || netstatus==JJ_UNBOUNDED ){
        opt = INF;
        goto FUERA;
    }

    if( netstatus != JJ_OPTIMAL && netstatus != 5/*11*/ && lpiter != MAXLPITER ){  //changed by Salome 08/02/12

#ifdef STAMP
        std::cout << "WARNING: network stat=" << JJgetstat(Nlp) << std::endl;
#endif
        opt = BIGVALUE;
    } else {
        JJgetobjval(Nlp,&opt);

        if( mode ){
            x = (double *)malloc(sizeof(double)*mac);
            if(x==NULL){                
                std::cout << "There is not enough memory for RED COST" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
            }
            if ( JJgetx(Nlp, x, 0, mac-1) ){
                std::cout << "ERROR: no solution is avalaible " << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
            }
            for(k=0;k<mac;k++)
                if( x[k]>ZERO ){
                    l = k/2;
                    columns[ net2cell[l] ].val = 1.0;
                    JJchgcoef(Nlp,-1,2*l,0.0);
                    JJchgcoef(Nlp,-1,2*l+1,0.0);
                }
            free( x );
            x = NULL; /*PWOF*/
        }
    }
FUERA:
    lu = 'L';
    bd = 0.0;
    JJchgbds(Nlp,1,&index1,&lu,&bd);
    lu = 'U';
    if( type < 0 )
        bd = var->uvalue;
    else 
        bd = var->lvalue;
    JJchgbds(Nlp,1,&index2,&lu,&bd);

//    lu = 'L';
//    bd = 0.0;
//    JJchgbds(Nlp,1,&index2,&lu,&bd);
//    lu = 'U';
//    if( type < 0 )
//        bd = var->lvalue;
//    else 
//        bd = var->uvalue;
//    JJchgbds(Nlp,1,&index1,&lu,&bd);

    return(opt);
}



static void remove_cell(int index,double *bd)

{
    int j;
    int    ind[2];
    char   lu[2] = {'B','B'};
    double value[2] = {0.0,0.0};
#ifdef STAMP
    if(cell2net[index] == -1){
        std::cout << "ERROR: variable " << index << " not in the NLP" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
#endif
    columns[index].val = 0.0;
    j = 2* cell2net[index];
    JJgetbdl(Nlp,bd,j,j+1);
    JJgetbdu(Nlp,bd+2,j,j+1);
    ind[0] = j;
    ind[1] = j+1;
    JJchgbds(Nlp,2,ind,lu,value);
}

static void restore_cell(int index,double *bd)

{
    int j;
    int ind[4];
    char lu[4] = {'L','L','U','U'};
#ifdef STAMP
    if(cell2net[index] == -1){
        std::cout << "ERROR: variable " << index << " not in the NLP" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
#endif
    columns[index].val = 1.0;
    j = 2* cell2net[index];
    ind[0] = ind[2] = j;
    ind[1] = ind[3] = j+1;
    JJchgbds(Nlp,4,ind,lu,bd);
}

int sort_cells(const void *p, const void *q/*int *p,int *q*/)

{
    int cp,cq;

    cp = columns[*(int *)p].weight;
    cq = columns[*(int *)q].weight;
    if( cp>cq ) return(-1);
    if( cp<cq ) return(1);
    return(0);
}
static void potential_innecesary(int    *nlist,int    *list)

{
    int  l,card;
    int  *rnumber;
    VARIABLE **rvar,**stack;
    char *status;
    //int  sort_cells(const void*,const void*);
    
    rnumber = (int *)calloc( (nrows+1) , sizeof(int) );
    if(rnumber==NULL){
        std::cout << "There is not enough memory for RNUMBER" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    rvar = (VARIABLE **)malloc( (nrows+1) * sizeof(VARIABLE*) );
    if(rvar==NULL){
        std::cout << "There is not enough memory for RVAR" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(l=0;l<nrows;l++)
        if( rows[l]->type==CAPACITY ){
            card = rows[l]->card;
            stack = rows[l]->stack;
            while(card--)
                if( stack[card]->val>1-ZERO && stack[card]->sensitive==0 ){
                    rnumber[l]++;
                    rvar[l] = stack[card];
                }
        }

    status = (char *)calloc( Rncells, sizeof(char) );
    if(status==NULL){
        std::cout << "There is not enough memory for STATUS" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(l=0;l<nrows;l++)
        if( rnumber[l]==1 ){
            status[rvar[l]->index]=1;
/*
#ifdef STAMP
            printf("CLEAN-UP: variable %d not need to be considered\n",l);
#endif
*/
        }
    free(rvar);
    rvar = NULL; /*PWOF*/
    free(rnumber);
    rnumber = NULL; /*PWOF*/
        
    *nlist = 0;
    for(l=0;l<Rncells;l++)
        if( columns[l].val>1-ZERO && columns[l].sensitive==0 && status[l]==0 )
            list[(*nlist)++] = l;
    free((void*)status);
    status = NULL; /*PWOF*/


    //qsort( (char *)list , *nlist , sizeof(int *) , sort_cells );
    //qsort( (char *)list , *nlist , sizeof(int) , sort_cells ); // list consists of int, not int*
    qsort( (void *)list , *nlist , sizeof(int) , sort_cells );
}



