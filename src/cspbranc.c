/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPBRANCH.c                                      */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "Cspdefns.h"
#include "CSPGLOB2.H"
#include "cspsolve.h"
#include "cspbranc.h"
#include "cspdebug.h"
#include "cspback.h"

/*  PROTOTYPES OF FUNCTIONS  */


static void       write_prob(void);
static void       insert_branch(long);      /* function for sort the tree   */
static VARIABLE   *variable_branching(int); /* choosing the branch-variable */
static void       constraint_branching(int,VARIABLE**,CONSTRAINT**,CONSTRAINT**);



#define NEAR05    1          /* with fractional value nearest to 0.5 */
#define BETTERLP  2          /* bigger PROP * LP_0 + (1-PROP) * LP_1 */

#define FRAC    0.5          /* maximum distance with respect to 0.5 */
#define PROP    0.5          /* proportion for weighting LP values   */
#define NUMB     10          /* number of frac. variables to study   */


/*  PRIVATE GLOBAL VARIABLES */

struct BRANCH *tree; /* array to save the branch-plant          */


/*  FUNCTIONS                */

int load_branch_tree()
{
    FILE     *pfile;

    mar = 1;
    mac = ncols;
    tree = NULL;

    pfile = fopen(fbranch,"w");
    if(pfile==NULL){          
          std::cout << "ERROR: not possible to write on " << fbranch << std::endl;
          CSPexit(EXIT_ERROR); //exit(1);
    }
    fclose(pfile);

    write_prob();
    return(0);
}

int unload_branch_tree()
{
    struct BRANCH  *ptr;

#ifdef STAMP
    if( tree ) std::cout << "ERROR: the branch-decision tree is not empty " << std::endl;
#endif

    while(tree){
        ptr    = tree->next;
        free(tree);
        tree   = ptr;
    }
    return(0);
}


static void write_prob()
{
    int      k,cont;
    VARIABLE *col;
    FILE     *pfile;
//    float    t1;
    
//    t1 = seconds();

    pfile = fopen(fbranch,"a");
    if(pfile==NULL){          
          std::cout << "ERROR: not possible to write on " << fbranch << std::endl;
          CSPexit(EXIT_ERROR); //exit(1);
    }
    fseek(pfile,0L,SEEK_END);
    insert_branch( ftell(pfile) );

    for(cont=k=0;k<ncols;k++)
         if( columns[k].stat != FIX_LB ) cont++;
    fprintf(pfile," %d",cont);
    col = columns;
    for(k=0;k<ncols;k++){
         if( col->stat != FIX_LB )
             fprintf(pfile," %p %d",col,col->stat);
         col++;
    }
    for( cont=0 , k=1 ; k<mar ; k++ )
         if( rind[k]->stat ==0 ) cont++;
    fprintf(pfile," %d",cont);
    for(k=1;k<mar;k++)
         if( rind[k]->stat ==0 )
             fprintf(pfile," %p",rind[k]);
    fclose(pfile);
}

static void   insert_branch(long  offset )

{
    struct BRANCH *branch,*ptr1,*ptr2;

    branch = (struct BRANCH *)malloc( sizeof(struct BRANCH) );
    if( branch==NULL ){        
        std::cout << "Not enought memory for branching "<< std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    branch->file   = offset;
    branch->val    = lowerb;

    if(tree==NULL || tree->val>lowerb){
        branch->next = tree;
        tree = branch;
        return;
    }
    ptr1 = tree;
    ptr2 = ptr1->next;
    while( ptr2!=NULL && ptr2->val<lowerb ){
        ptr1 = ptr2;
        ptr2 = ptr1->next;
    }
    branch->next = ptr2;
    ptr1->next = branch;
}


int read_prob()
{
    int            k,l,lcuts,cont;
    CONSTRAINT     **stack = NULL;
    long           offset;
    struct BRANCH  *ptr;
    VARIABLE       *col;
    FILE           *pfile;
//    float          t1;

    do{
        if( tree==NULL ) return(0);
        offset = tree->file;
        lowerb = tree->val;
        ptr    = tree->next;
        free(tree);
        tree   = ptr;
    }while( ceil(lowerb-ZERO)+ZERO > upperb );
#ifdef STAMP
    std::cout << "  >>>>>>>>>>>>> reading problem " << (float)lowerb << std::endl;
#endif
 //   t1 = seconds();
    pfile = fopen(fbranch,"r");
    if(pfile==NULL){          
          std::cout << "ERROR: not possible to write on " << fbranch << std::endl;
          CSPexit(EXIT_ERROR); //exit(1);
    }
    fseek(pfile,offset,SEEK_SET);

    for(k=0;k<ncols;k++)
         columns[k].stat = FIX_LB;
    if (fscanf(pfile,"%d",&cont) != 1) goto invalid_branch;
    for(k=0;k<cont;k++){
         if (fscanf(pfile,"%p %d",&col,&l) != 2) goto invalid_branch;
         col->stat = l;
    }
    if (fscanf(pfile,"%d",&lcuts) != 1) goto invalid_branch;

    if(lcuts){
         stack = (CONSTRAINT **)malloc( lcuts * sizeof( CONSTRAINT *) );
         if( stack==NULL ){              
              std::cout << "ERROR: not enought memory for STACK" << std::endl;
              CSPexit(EXIT_MEMO); //exit(1);
         }
         for(k=0;k<lcuts;k++){
              if (fscanf(pfile,"%p",stack+k) != 1) goto invalid_branch;
              stack[k]->stat = 0;
         }
    }
    fclose(pfile);

    load_lp();
    if(lcuts){
        add_rows(lcuts,stack);
        free( (void *)stack );
        stack = NULL; /*PWOF*/
    }
    put_base();
    setup_lp();
    get_solution();
    return(1);

invalid_branch:
    free(stack);
    fclose(pfile);
    CSPexit(EXIT_ERROR);
    return 0;
}

int var_branching()
{
    VARIABLE *col;
    double   old_lowerb = lowerb;

#ifdef STAMP
    std::cout << "  >>>>>>>>>>>>> var branching " << std::endl;
    if(branchs==0) write_sol(fsolution);
#endif

    if( upperb-ZERO < ceil(lowerb-ZERO) )
        return(1);
    col = variable_branching(BETTERLP);
    if( upperb-ZERO < ceil(lowerb-ZERO) )
        return(1);
    if(col==NULL)
        return(0);
#ifdef STAMP
    std::cout << branchs << " on variable " << col->index << "=" << col->val << std::endl;
#endif
    branchs++;
/* making the left child */

    if( ceil( solve_child_col( col , (double)0 , 1)-ZERO) < upperb-ZERO ) {
        col->stat = FIX_LB;
        write_prob();
    }

/* making the left child */

    if( ceil( solve_child_col( col , (double)1 , 1)-ZERO) < upperb-ZERO ) {
        col->stat = FIX_UB;
        write_prob();
    }
    lowerb = old_lowerb;
    return(1);
}


static VARIABLE *variable_branching(int type)

{
    int      k,card;
    double   val,val0,val1,lb0,lb1,minim;
    VARIABLE *col;
    VARIABLE **stack;
    int      sort_var(const void*,const void*);

    col = NULL;

    card = 0;
    for(k=0;k<nsupport;k++)
        if( fabs( support[k]->val - 0.5 ) < FRAC-ZERO ) card++;
    if( card==0 ){        
        std::cout << "ERROR: no fractional variables (FRAC=" << FRAC << std::endl;
        write_sol(fsolution);
        CSPexit(EXIT_ERROR); //exit(1);
    }

    switch(type){
    case NEAR05:
        minim = FRAC;
        for(k=0;k<nsupport;k++){
            val = fabs( support[k]->val - 0.5 );
            if( val < minim ){
                minim = val;
                col   = support[k];
            }
        }
        return( col );
    case BETTERLP:
        stack = (VARIABLE **)malloc( card * sizeof( VARIABLE * ) );
        if( stack==NULL ){            
            std::cout << "ERROR: not memory for 'frac' in branching" << std::endl;
            CSPexit(EXIT_ERROR); //exit(1);
        }
        card = 0;
        for(k=0;k<nsupport;k++)
            if( fabs( support[k]->val - 0.5 ) < FRAC-ZERO )
                stack[card++] = support[k];
        //qsort( (char *)stack , card , sizeof(VARIABLE *) , sort_var );
        qsort( (void *)stack , card , sizeof(VARIABLE *) , sort_var );

        if(card > NUMB) card = NUMB;

        minim = lowerb - 1;
        for(k=0;k<card;k++){
            activa_pricing();
            val0 = solve_child_col( stack[k] , (double)0 , 0);
            lb0  = ceil( val0 -ZERO ) + ZERO;
            activa_pricing();
            val1 = solve_child_col( stack[k] , (double)1 , 0);
            lb1  = ceil( val1 -ZERO ) + ZERO;
            if( lb0 > upperb && lb1 > upperb ) {
                lowerb = upperb;
                free( stack );
                stack = NULL; /*PWOF*/
                return(NULL);
            }
            if( lb0 > upperb ) {
                del_col( stack[k] , FIX_UB );
                free( stack );
                stack = NULL; /*PWOF*/
                return(NULL);
            }
            if( lb1 > upperb ) {
                del_col( stack[k] , FIX_LB );
                free( stack );
                stack = NULL; /*PWOF*/
                return(NULL);
            }
            val  = PROP * val0 + (1-PROP) * val1;
            if( val > minim ){
                minim = val;
                col   = stack[k];
            }
        }
        free( stack );
        stack = NULL; /*PWOF*/
        return( col );
    default:        
        std::cout << "ERROR: unknown type of branching" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    return( NULL );
}

int sort_var(const void *p, const void *q/*VARIABLE **p,VARIABLE **q*/)

{
    double vp,vq;
    int cp,cq;

    vp = fabs( (*(VARIABLE **)p)->val -0.5);
    vq = fabs( (*(VARIABLE **)q)->val -0.5);
    if( vp>vq ) return(1);
    if( vp<vq ) return(-1);
    cp = (*(VARIABLE **)p)->weight;
    cq = (*(VARIABLE **)q)->weight;
    if( cp<cq ) return(1);
    if( cp>cq ) return(-1);
    return(0);
}

int con_branching()
{
    CONSTRAINT *row0,*row1;
    double     old_lowerb = lowerb;


#ifdef STAMP
    std::cout << "  >>>>>>>>>>>>> constraint branching " << std::endl;
    if(branchs==0) write_sol(fsolution);
#endif

    if( upperb < ceil(lowerb-ZERO)+ZERO ) return(1);
////    constraint_branching(nbetter,better,&row0,&row1);
    constraint_branching(nsupport,support,&row0,&row1);
    if( upperb < ceil(lowerb-ZERO)+ZERO ) return(1);
    if(row0==NULL && row1==NULL)return(-1);
    if(row0==NULL || row1==NULL)return(0);


/***
 to do not use the branch-cuts
***/
////    return(-1);



    branchs++;
#ifdef STAMP
    std::cout << "branching " << branchs << " on constraint" << std::endl;
#endif


/* making the left child */

    row0->stat = LP_BA;
    add_rows(1,&row0);
    lowerb = solve_lp(0);
    get_solution();
    get_base();
    write_prob();
    deletelastrow();
    row0->stat = FIX_LB;

/* making the left child */


    row1->stat = LP_BA;
    add_rows(1,&row1);
    lowerb = solve_lp(0);
    get_solution();
    get_base();
    write_prob();
    deletelastrow();
    row1->stat = FIX_LB;

    lowerb = old_lowerb;
    return(1);
}

static void constraint_branching(int nvar,VARIABLE **xvar,CONSTRAINT **con0,CONSTRAINT **con1)

{
    int        k,card;
    double     rhs,rhs0,rhs1,lb0,lb1;
    VARIABLE   **x;
    double     *c;
    CONSTRAINT *row0,*row1;

    *con0 = *con1 = NULL;

    rhs  = 0;
    card = 0;
    for(k=0;k<nvar;k++)
        if( fabs( xvar[k]->val - 0.5 ) < FRAC-ZERO ){
            rhs += xvar[k]->val;
            card++;
        }

    if( card==0 ){
#ifdef STAMP
        std::cout << "Warning: no fractional variables for branch-constraint" << std::endl;
#endif
        return;
    }
    rhs0 = ceil(rhs-ZERO);
    rhs1 = floor(rhs+ZERO);
#ifdef STAMP
    std::cout << rhs1 << "<" << rhs << "<" << rhs0 << std::endl;
#endif
    if( rhs0 - rhs1 < ZERO ){
#ifdef STAMP
        std::cout << "Warning: no fractional branch-constraint (sum=" << rhs << std::endl;
#endif
        return;
    }
    x = (VARIABLE **)malloc( card*sizeof(VARIABLE *) );
    if( x==NULL ){        
        std::cout << " ERROR: not enough memory for the SORT" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    c = (double *)malloc( card*sizeof(double) );
    if( c==NULL ){        
        std::cout << " ERROR: not enough memory for the SORT" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    card = 0;
    for(k=0;k<nvar;k++)
        if( fabs( xvar[k]->val - 0.5 ) < FRAC-ZERO ){
            x[card] = xvar[k];
            c[card] = 1.0;
            card++;
        }
    row0 = (CONSTRAINT *)malloc( sizeof(CONSTRAINT) );
    if(row0==NULL){        
        std::cout << "There is not enough memory for ROW" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    row0->rhs    = rhs0;
    row0->card   = card;
    row0->stack  = x;
    row0->coef   = c;
    row0->sense  = 'G';
    row0->type   = BRANCHCUT;
    row0->stat   = FIX_LB;
    row0->lp     = -1;
    row0->con    = NULL;


    row1 = (CONSTRAINT *)malloc( sizeof(CONSTRAINT) );
    if(row1==NULL){        
        std::cout << "There is not enough memory for ROW" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    row1->rhs    = rhs1;
    row1->card   = card;
    row1->stack  = x;
    row1->coef   = c;
    row1->sense  = 'L';
    row1->type   = BRANCHCUT;
    row1->stat   = FIX_LB;
    row1->lp     = -1;
    row1->con    = NULL;

    add_rows(1,&row0);
    activa_pricing();
    lb0 = ceil( solve_lp(0) - ZERO ) + ZERO;
    deletelastrow();
    add_rows(1,&row1);
    activa_pricing();
    lb1 = ceil( solve_lp(0) - ZERO ) + ZERO;
    deletelastrow();


    if( lb0 > upperb && lb1 > upperb ) {
        lowerb = upperb;
        free( row0 );
        row0 = NULL; /*PWOF*/
        free( row1 );
        row1 = NULL; /*PWOF*/
        return;
    }
    if( lb0 > upperb ) {
        row1->stat = LP_BA;
        add_rows(1,&row1);
        free( row0 );
        row0 = NULL; /*PWOF*/
        rows[nrows++] = *con1 = row1;
        return;
    }
    if( lb1 > upperb ) {
        row0->stat = LP_BA;
        add_rows(1,&row0);
        rows[nrows++] = *con0 = row0;
        free( row1 );
        row1 = NULL; /*PWOF*/
        return;
    }
    rows[nrows++] = *con0 = row0;
    rows[nrows++] = *con1 = row1;
}
    
double     violation_branch(CONSTRAINT *con)

{
    int       k;
    double    viola;
    VARIABLE  **xvar;
    double    *cvar;

    if( con->type != BRANCHCUT ){        
        std::cout << " ERROR: not branching constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    xvar  = con->stack;
    cvar  = con->coef;
    k     = con->card;
    viola = con->rhs;
    while(k--)
        viola -= cvar[k] * xvar[k]->val;
    if( con->sense=='L' ) viola = -viola;
    return( viola );
}


double     get_coeficient_branching(VARIABLE   *col,CONSTRAINT *con)

{
    int       k;
    VARIABLE  **xvar;
    double    *cvar;

    xvar = con->stack;
    cvar = con->coef;
    k = con->card;
    while(k--)
        if( xvar[k] == col )
            return( cvar[k] );
    return 0;
}


