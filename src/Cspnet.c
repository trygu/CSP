/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPNET.c                                         */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/


/* #define CHECKLP */


#include <stdio.h>
#include <stdlib.h>
#include "Cspdefns.h"
#include "CSPGLOB2.H"
#include "Jjsolver.h"
#ifdef CHECKLP
#include "\cplex50\check.c"
#endif
#include "cspnet.h"
#include "cspback.h"



/* PROTOTYPES OF FUNCTIONS */



/* PRIVATE GLOBAL DATA TO INTERPRETATE THE DUAL SOLUTION */

static int    *net2cell;   /* cell number in the network of a node         */
static int    *cell2net;   /* node associated with a cell of the network   */
static int    *net2sum;    /* sum number in the network of a node          */
static int    *sum2net;    /* node associated with a sum of the network    */

/* PRIVATE GLOBAL DATA FOR LP */

static char     Nprobname[]="N-NETWORK";
static double   *Nobjx;
static double   *Nrhsx;
static char     *Nsenx;
static int      *Nmatbeg;
static int      *Nmatind;
static int      *Nmatcnt;
static double   *Nmatval;
static double   *Nbdl;
static double   *Nbdu;
static char     *Nxctype;
static JJLPptr  Nlp;

/*====================================================================*/


int      protected_flow(int nlist,VARIABLE **list)

{
    int    k;
    double *status;

    status = (double *)calloc(sizeof(double),Rncells);
    if(status==NULL){        
        std::cout << "There is not enough memory for STATUS" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(k=0;k<nlist;k++) status[ list[k]->index ]=1;
    k = protected1(status);
    free(status);
    status = NULL; /*PWOF*/
    return(k);
}


int    protected1(double *status)

{
    int        k,l;
    double     *primal;
    char       *sleep;
    PROT_LEVEL *pro,*pro0;

    load_network(status,'C');
    sleep = (char *)calloc( nprot_level , sizeof(char) );
    primal = (double *)malloc( Rncells * sizeof(double) );
    l = nprot_level;
    while(l--){
        pro = prot_level+l;
        if( !sleep[l] ){
            if( protection_level(pro->sen->var,pro->sense,NULL,NULL,NULL,primal,'C') < pro->sense * pro->sen->var->nominal + pro->level-ZERO )
                break;
            k = l;
            while(k--){
                pro0 = prot_level+k;
                if( pro0->sense * primal[ pro0->sen->var->index ] > pro0->sense * pro0->sen->var->nominal + pro0->level - ZERO )
                    sleep[k] = 1;
            }
        }
    }
    free( primal );
    primal = NULL; /*PWOF*/
    free( sleep );
    sleep = NULL; /*PWOF*/
    unload_network();
    if(l<0) return(1);
#ifdef STAMP    
    std::cout << " violated protection level " << l << std::endl;
#endif
    return(0);
}


/**
***   Set on the network.
**/


void   load_network(double *status,char type)
    //type   C = real numbers  ;   I = integer numbers
{
        int      i,l;
        int      mar,mac,macsz,marsz,matsz;
        VARIABLE *var;
        struct   CELDA *c;

/*
        printf("\nNetwork with NR=%d and NC=%d cols ;",NR,NC);
*/

#ifdef PARTIAL
        for(i=0;i<Rncells;i++)
            if(columns[i].sensitive && status[i]<ZERO) status[i]=2*ZERO;
#endif


/* allocation memory */

        macsz = Rncells;
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

		if( type=='I' ){
            if((Nxctype=(char *)malloc(sizeof(char)*macsz))==NULL){                
                std::cout << "There is not enough memory for XCTYPE" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
            }
		}

/* LP row construction */

        for(i=0;i<Rnsums;i++) sum2net[i]=0;
        for(i=0;i<Rncells;i++)
            if( status[i]>ZERO )
                for( c=columns[i].next ; c ; c=c->next )
                    sum2net[ c->index ]++;

        mar=0;
        for(i=0;i<Rnsums;i++)
            if( sum2net[i] ){
                sum2net[i] = mar;
                net2sum[mar] = i;
                Nrhsx[mar]   = Rrhs[i];
                Nsenx[mar]   = 'E';
                mar++;
            } else
                sum2net[i] = -1;

/* LP column construction */

        mac = 0;
        l   = 0;
        for(i=0;i<Rncells;i++){
            var = columns+i;
            if( status[i]>ZERO ){
                net2cell[mac]  = i;
                cell2net[i]    = mac;
                Nobjx[mac]     = 0;
                Nmatbeg[mac]   = l;
                Nmatcnt[mac]   = 0;
                for( c=var->next ; c ; c=c->next ){
                    Nmatcnt[mac]++;
                    Nmatval[l]     = c->coef;
                    Nmatind[l++]   = sum2net[ c->index ];
                    if( l==matsz ){                        
                        std::cout << "ERROR: not enough space " << std::endl;
                        CSPexit(EXIT_ERROR); //exit(1);
                    }
                }
                Nbdl[mac]      = var->nominal - var->lvalue * status[i];
                Nbdu[mac]      = var->nominal + var->uvalue * status[i];
                mac++;
            } else {
                cell2net[i]    = -1;
                for( c=var->next ; c ; c=c->next )
                    if( sum2net[c->index] != -1 )
                        Nrhsx[ sum2net[c->index] ] -= c->coef * var->nominal;
            }
        }

       //printf("  mac=%d   mar=%d  ",mac,mar);

       /* LP loading */


#ifdef CHECKLP
        l = JJcheckprob (Nprobname, mac, mar, 0, -1, Nobjx, Nrhsx,
                       Nsenx, Nmatbeg, Nmatcnt, Nmatind, Nmatval,
                       Nbdl , Nbdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       macsz, marsz, matsz, 0, 0, 0, 0, 0, NULL);
        if(l)CSPexit(EXIT_ERROR); //exit(1);
#endif

        Nlp = JJloadprob (Nprobname, mac, mar, 0, -1, Nobjx, Nrhsx,
                       Nsenx, Nmatbeg, Nmatcnt, Nmatind, Nmatval,
                       Nbdl , Nbdu , NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL,
                       NULL, NULL, NULL, NULL, NULL, NULL,
                       macsz, marsz, matsz, 0, 0, 0, 0, 0);

        if ( Nlp == NULL ){                
                std::cout << "ERROR: There is not enough memory for the first LP" << std::endl;
                CSPexit(EXIT_LPSOLVER); //exit(1);
        }
#ifdef STAMP    
        JJsetscr_ind(Nlp,0);
#endif
        if( type == 'I'){
			for(i=0;i<mac;i++)
				Nxctype[i] = 'I';
			l = JJcopyctype( Nlp , Nxctype );
			free(Nxctype);
                        Nxctype = NULL; /*PWOF*/ 
		}


#ifdef CHECKLP
        JJmpswrite(Nlp,fmpsnet);
#endif
}

/**
***   Set off the network.
**/

void unload_network()
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
*** Bounding the lower/upper cell variations
**/

void   bounding_bds(int bound)

{
    int    k,l,mac;
    VARIABLE *var;
    int    *index;
    char   *type;
    double *value;

    mac = JJgetmac(Nlp);
    index=(int *)malloc(sizeof(int)*2*mac);
    if(index==NULL){        
        std::cout << "There is not enough memory for INDEX" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    type=(char *)malloc(sizeof(char)*2*mac);
    if(type==NULL){        
        std::cout << "There is not enough memory for TYPE" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    value=(double *)malloc(sizeof(double)*2*mac);
    if(value==NULL){        
        std::cout << "There is not enough memory for VALUE" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(l=k=0;k<mac;k++){
        var = columns + net2cell[k];
        index[l] = k;
        type[l]  = 'L';
        value[l] = var->nominal - MIN( bound , var->lvalue ) * var->val;
        l++;
        index[l] = k;
        type[l]  = 'U';
        value[l] = var->nominal + MIN( bound , var->uvalue ) * var->val;
        l++;
    }
    JJchgbds(Nlp,l,index,type,value);
    free(index);
    index = NULL; /*PWOF*/
    free(type);
    type = NULL; /*PWOF*/
    free(value);
    value = NULL; /*PWOF*/
}
 
void   bounding_1()
{
    int    k,l,mac;
    VARIABLE *var;
    int    *index;
    char   *type;
    double *value;

    mac = JJgetmac(Nlp);
    index=(int *)malloc(sizeof(int)*2*mac);
    if(index==NULL){        
        std::cout << "There is not enough memory for INDEX" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    type=(char *)malloc(sizeof(char)*2*mac);
    if(type==NULL){        
        std::cout << "There is not enough memory for TYPE" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    value=(double *)malloc(sizeof(double)*2*mac);
    if(value==NULL){        
        std::cout << "There is not enough memory for VALUE" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(l=k=0;k<mac;k++){
        var = columns + net2cell[k];
        index[l] = k;
        type[l]  = 'L';
        value[l] = var->nominal - var->val;
        l++;
        index[l] = k;
        type[l]  = 'U';
        value[l] = var->nominal + var->val;
        l++;
    }
    JJchgbds(Nlp,l,index,type,value);
    free(index);
    index = NULL; /*PWOF*/
    free(type);
    type = NULL; /*PWOF*/
    free(value);
    value = NULL; /*PWOF*/
}
 




/**
***  Compute the minimum ('type'=1) or maximum ('type'=-1)
***  protection level 'opt' of a suppressed cell.
***  The optimal dual solution is returned if it is requered.
**/

double    protection_level(VARIABLE *var,int goal,int *nvar,VARIABLE **lvar,double *cvar,double *primal,char type)
 //type 'C' = real numbers  ;  'I' = integer numbers
{
    int       k,i,num,mar,mac,card;
    int       netstatus,netnodes,netarcs,netiter;
    double    opt,valor;
    double    *pi,*x,*cost;
    struct    CELDA *c;
    std::string netlpname=fsdcnetlp;

    num = cell2net[var->index];
    if( num==-1 ){        
        std::cout << "ERROR: variable " << var->name << " not in the network" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }


    if( goal < 0 )
        JJchgcoef(Nlp,-1,num,-1.0);
    else
        JJchgcoef(Nlp,-1,num,1.0);

    if( type == 'I' ){
	JJlpwrite(Nlp,netlpname);
        if ( JJmipopt(Nlp) ){            
            std::cout << " it was not possible to solve with mip-opt " << std::endl;
            JJlpwrite(Nlp,netlpname);
            CSPexit(EXIT_LPSOLVER); //exit(1);
        }
        k = JJgetstat(Nlp);
	    JJgetobjval(Nlp,&opt);
        JJchgcoef(Nlp,-1,num,0.0);
		return opt;
    }

    if ( JJnetopt(Nlp,&netstatus,&netnodes,&netarcs,&netiter) ){        
        std::cout << " it was not possible to solve with NETOPT " << std::endl;
        JJlpwrite(Nlp,netlpname);
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    
/*
    if( netstatus!=CPX_NETOPTIMAL || netnodes!=JJgetmar(Nlp) || netarcs!=JJgetmac(Nlp) ){
        printf(" WARNING in sdcnet.c - netopt\n");
        printf(" netstatus=%d netnodes=%d netarcs=%d netiter=%d mar=%d mac=%d\n",
             netstatus,netnodes,netarcs,netiter,JJgetmar(Nlp),JJgetmac(Nlp));
    }
*/
    if ( JJdualopt(Nlp) ){        
        std::cout << " it was not possible to solve with dualopt " << std::endl;
        JJlpwrite(Nlp,netlpname);
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }

#ifdef STAMP
    k = JJgetstat(Nlp);
    if( k!=JJ_OPTIMAL && k!=5/*11*/ ){  // changed by Salome 08/02/12
        std::cout << "ERROR: attacker problem with CPX status = " << k  << std::endl;
        JJlpwrite(Nlp,netlpname);
        CSPexit(EXIT_ERROR); //exit(1);
    }
#endif    
    mar = JJgetmar(Nlp);    // Number of rows (constraints)
    mac = JJgetmac(Nlp);    // Number de colums (variable)
    JJgetobjval(Nlp,&opt);  // Value of obj. solution (opt)
/*
    printf(" ... solving LP:stat=%d opt=%f tit=%d Iit=%d mac=%d mar=%d nz=%d\n",
JJgetstat(Nlp),(float)opt,JJgetitc(Nlp),JJgetitci(Nlp),JJgetmac(Nlp),JJgetmar(Nlp),JJgetmat(Nlp));
*/

    if(nvar){
        cost = (double *)malloc(sizeof(double)*Rncells);
        if(cost==NULL){            
            std::cout << "There is not enough memory for RED COST" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        pi = (double *)malloc(sizeof(double)*mar);
        if(pi==NULL){            
            std::cout << "There is not enough memory for RED COST" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        if ( JJgetpi(Nlp, pi, 0, mar-1) ){    // pi gets dual values
            std::cout << "ERROR: no solution is avalaible " << std::endl;
            CSPexit(EXIT_LPSOLVER); //exit(1);
        }

        card = 0;
        for(k=0;k<Rncells;k++){
            valor = 0;
            for( c=columns[k].next ; c ; c=c->next ){
                i = sum2net[c->index];
                if( i!=-1 )
                    valor -= pi[i] * c->coef;
            }
            cost[k] = valor;            
        }
        //JJgetdj(Nlp,cost,0,Rncells);
        cost[ var->index ] += goal;

        for(k=0;k<Rncells;k++){
            if (cost[k]>ZERO && columns[k].uvalue ){
                lvar[card] = &columns[k];
                cvar[card] = cost[k] * columns[k].uvalue;
                card++;
            } else if( -cost[k]>ZERO && columns[k].lvalue ) {
                lvar[card] = &columns[k];
                cvar[card] = - cost[k] * columns[k].lvalue;
                card++;
            }
            //if (cost[k]>ZERO) printf(" ind %d uval %lf ",k,columns[k].uvalue);
            //if (-cost[k]>ZERO) printf(" ind %d lval %lf ",k,columns[k].lvalue);
        }        
        *nvar = card;
        free( cost );
        cost = NULL; /*PWOF*/
        free( pi );
        pi = NULL; /*PWOF*/
    }
    if( primal ){
        for(k=0;k<Rncells;k++) primal[k] = columns[k].nominal;
        x = (double *)malloc(sizeof(double)*mac);
        if(x==NULL){            
            std::cout << "There is not enough memory for RED COST" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        if ( JJgetx(Nlp, x, 0, mac-1) ){    // x contains the solution values ​​of the variables listed            
            std::cout << "ERROR: no solution is avalaible " << std::endl;
            CSPexit(EXIT_LPSOLVER); //exit(1);
        }
        for(k=0;k<mac;k++) primal[ net2cell[k] ] = x[k];
        free( x );
        x = NULL; /*PWOF*/
    }

    JJchgcoef(Nlp,-1,num,0.0);

    return(opt);
}

void free_col(int index,double *bd)

{
    int j;
    int    ind[2];
    char   lu[2] = {'L','U'};
    double value[2];
#ifdef STAMP
    if(cell2net[index] == -1){
        std::cout << "ERROR: variable " << index << " not in the NLP" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
#endif
    j = cell2net[index];
    JJgetbdl(Nlp,bd,j,j);
    JJgetbdu(Nlp,bd+1,j,j);
    ind[0] = ind[1] = j;
#ifdef VSCIP
    value[0] = -1e20;
    value[1] =  1e20;
#endif
#ifdef CPLEX5 // Also CPLEX7
    value[0] = -CPX_INFBOUND; //-INFBOUND;
    value[1] = CPX_INFBOUND;  //INFBOUND;
#endif
#ifdef XPRESS
    value[0] = -XPRS_PLUSINFINITY;
    value[1] = XPRS_PLUSINFINITY;
#endif
#ifdef XPRESS_13
    value[0] = -XPRS_PLUSINFINITY;
    value[1] = XPRS_PLUSINFINITY;
#endif
    JJchgbds(Nlp,2,ind,lu,value);
}

void unfree_col(int index,double *bd)

{
    int ind[2];
    char lu[2] = {'L','U'};
#ifdef STAMP
    if(cell2net[index] == -1){
        std::cout << "ERROR: variable " << index << " not in the NLP" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
#endif
    ind[0] = ind[1] = cell2net[index];
    JJchgbds(Nlp,2,ind,lu,bd);
}

