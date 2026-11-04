/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPDEBUG.c                                       */
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
#include "cspdebug.h"
#include "cspsep.h"
#include "cspprice.h"
#include "cspnet.h"
#include "cspcover.h"
#include "cspback.h"



static int    in_better(VARIABLE *);
static void   testing(CONSTRAINT *);


void write_sol(char *text)

{
    int      k;
    FILE     *fich;
    VARIABLE *col;

    fich=fopen(text,"w");
    fprintf(fich,"2 %d %d\n%d\n",Rncells,Rnsums,nsensitive);
    for(k=0;k<nsensitive;k++){
        col = sensitive[k].var;
        if( col->name )
            fprintf(fich,"%12s",col->name);
        else
            fprintf(fich,"%d",col->index);
        fprintf(fich," %f %f %f\n",sensitive[k].lpl,sensitive[k].upl,col->nominal);
    }
    fprintf(fich,"%d\n",nsupport-nsensitive);
    for(k=0;k<nsupport;k++){
        col = support[k];
        if( !col->sensitive ){
            if( col->name )
                fprintf(fich,"%12s",col->name);
            else
                fprintf(fich,"%d",col->index);
            fprintf(fich," %f\n",(float)col->val);
        }
    }
    fprintf(fich,"zopt = %f\n",(float)lowerb);
    fprintf(fich,"time = %f\n",(float)seconds()-t0);
    fclose(fich);
}

void write_heu(const char *text)

{
    int      k,weight;
    FILE     *fich;
    VARIABLE *col;

    fich=fopen(text,"w");
    fprintf(fich,"%d %d\n%d\n",Rncells,Rnsums,nsensitive);
    for(k=0;k<nsensitive;k++){
        col = sensitive[k].var;
        if( col->name )
            fprintf(fich,"%12s",col->name);
        else
            fprintf(fich,"%d",col->index);
        fprintf(fich," %f %f %f\n",sensitive[k].lpl,sensitive[k].upl,col->nominal);
    }
    fprintf(fich,"%d\n",nbetter-nsensitive);
    weight = 0;
    for(k=0;k<nbetter;k++){
        col = better[k];
        if( !col->sensitive ){
            if( col->name )
                fprintf(fich,"%12s",col->name);
            else
                fprintf(fich,"%d",col->index);
            fprintf(fich," 1.0 %f\n",col->nominal);
        }
        weight += col->weight;
    }
    fprintf(fich,"%f\n",(float)upperb);
    fclose(fich);
    if( weight != upperb ){        
        std::cout << "ERROR:  weight=" << weight << "   upperb=" << upperb << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
}

void print_sol()
{
    int      k;
    VARIABLE *col;
    for(k=0;k<nsupport;k++){
        col = support[k];
        if( !col->sensitive ){
            if( col->name )
                std::cout << col->name;
            else
                std::cout << col->index;
            std::cout << (float)col->val << std::endl;
        }
    }
    std::cout << "zopt = " << (float)lowerb << std::endl;
}

int control_ind()
{
    int i,lpt,stat;
    double val;
/*
    printf(" :::: control mac=%d mar=%d\n",mac,mar);
*/
    for(i=0;i<ncols;i++){
        stat = columns[i].stat;
        lpt  = columns[i].lp;
        val  = columns[i].val;
        switch( stat ){
            case WAITING:
/**
                 if( lpt != -1 ){
                     printf("ERROR 1c: var=%d  with stat=%d lp=%d\n",i,stat,lpt);
                     return(1);
                 }
**/
                 if( val>ZERO ){
                     std::cout << "ERROR 1c: var=" << i << " with stat=" << stat << " val=" << val << std::endl;
                     return(1);
                 }
                 break;
            case FIX_LB:
/**
                 if( lpt != -1 ){
                     printf("ERROR 1a: var=%d  with stat=%d lp=%d\n",i,stat,lpt);
                     return(1);
                 }
**/
                 if( val>ZERO ){
                     std::cout << "ERROR 1b: var=" << i << " with stat=" << stat << " val=" << val << std::endl;
                     return(1);
                 }
                 break;
            case FIX_UB:
/**
                 if( lpt != -1 ){
                     printf("ERROR 1c: var=%d  with stat=%d lp=%d\n",i,stat,lpt);
                     return(1);
                 }
**/
                 if( val<1-ZERO ){
                     std::cout << "ERROR 1d: var=" << i << " with stat=" << stat << " val=" << val << std::endl;
                     return(1);
                 }
                 break;
            case LP_LB:
            case LP_UB:
            case LP_BA:
                 if( lpt<0 || lpt>=mac ){
                     std::cout << "ERROR 2: var=" << i << " with stat=" << stat << " lp=" << lpt << std::endl;
                     return(1);
                 }
                 if( cind[lpt] != columns+i ){
                     std::cout << "ERROR 3: var=" << i << " with lp=" << lpt << std::endl;
                     return(1);
                 }
                 break;
            default:
                 std::cout << "ERROR 4: " << i << " with unknown stat=" << stat << std::endl;
                 return(1);
        }
    }
    for(i=0;i<nrows;i++){
        stat = rows[i]->stat;
        lpt  = rows[i]->lp;
        switch( stat ){
            case FIX_LB:
            case FIX_UB:
                 if( lpt != -1 ){
                     std::cout << "ERROR 5: row=" << i << " with stat=" << stat << " lp=" << lpt << std::endl;
                     return(1);
                 }
                 break;
            case LP_LB:
            case LP_UB:
            case LP_BA:
                 if( lpt<0 || lpt>=mar ){
                     std::cout << "ERROR 6: row=" << i << " with stat=" << stat << " lp=" << lpt << std::endl;
                     return(1);
                 }
                 if( rind[lpt] != rows[i] ){
                     std::cout << "ERROR 7: row=" << i << " with lp=" << lpt << std::endl;
                     return(1);
                 }
                 break;
            default:
                 std::cout << "ERROR 8: " << i << " with unknown stat=" << stat << std::endl;
                 return(1);
        }
    }

    for(i=1;i<mac;i++)
        if(i!=cind[i]->lp){
            std::cout << "ERROR 5: col i=" << i << "  lp=" << cind[i]->lp << std::endl;
            return(1);
        }
    for(i=1;i<mar;i++)
        if(i!=rind[i]->lp){
            std::cout << "ERROR 6: row i=" << i << "  lp=" << rind[i]->lp << std::endl;
            return(1);
        }
    return(0);
}


void print_pool()
{
    int i;
    for(i=0;i<nrows;i++)
        print_row( rows[i] );
}

void print_lp()
{
    int i;
    for(i=0;i<mar;i++)
        print_row( rind[i] );
}


int        print_row(CONSTRAINT *con)

{
    int      card=0;
    double   *coef=NULL;
    VARIABLE **stack=NULL;

    if(con==NULL) return(1);

    std::cout << " Row type=" << con->type << "  lp=" << con->lp << "  vio=" << violated(con) << std::endl;

    switch( con->type ){

    case CAPACITY:
    case BRIDGE:
        stack = con->stack;
        coef  = con->coef;
        card  = con->card;
        break;
    default:        
        std::cout << "ERROR print: unknown type " << con->type << " of constraint in paint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    while(card--){
        std::cout << coef[card];
        print_col(stack[card]);
    }

    switch( con->type ){
    case CAPACITY:
    case BRIDGE:
        break;
    default:        
        std::cout << "ERROR print: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    switch( con->sense ){
    case 'G': std::cout << " > ";
              break;
    case 'E': std::cout << " = ";
              break;
    case 'L': std::cout << " < ";
              break;
    default:  std::cout << " extrange" << std::endl;
              CSPexit(EXIT_ERROR);
              //exit(1);
    }

    std::cout << (float)con->rhs << std::endl;
    return(0);
}

int print_col(VARIABLE *var)

{
    if(var==NULL) return(1);
//    printf(" {%d,%lf}",var->index,var->val);
    std::cout << var->name;
    return(0);
}

void print_card_cover()
{
    int i;

    std::cout << "  card of covers in pool: ";
    for(i=0;i<nrows;i++)
        if( rows[i]->type==COVER )
            std::cout << rows[i]->card;
    std::cout << "\n";
}


static int in_better(VARIABLE * col )

{
    int k;
    for(k=0;k<nbetter;k++) if( col==better[k] ) return(1);
    return(0);
}


void read_heuristic(char *name)

{
    int    i,col;
    FILE   *file;

    file = fopen(name,"r");
    if(file==NULL){        
        std::cout << "ERROR: not possible to open " << name << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    if (fscanf(file,"%d\n",&nbetter) != 1) {
        fclose(file);
        CSPexit(EXIT_ERROR);
        return;
    }
    upperb = 0;
    for(i=0;i<nbetter;i++){
        if (fscanf(file,"%d\n",&col) != 1 || col < 0 || col >= ncols) {
            fclose(file);
            CSPexit(EXIT_ERROR);
            return;
        }
        std::cout << col << " ;";
        better[i] = columns + col;
        upperb += columns[ col ].weight;
    }
    std::cout << " ... feasible heuristic of value " << upperb << std::endl;
    fclose(file);

    if( !protected_flow(nbetter,better) ){        
        std::cout << "ERROR: not feasible solution" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
}

void writing_lp()
{
    int    k;

    for(k=0;k<mac;k++)
        print_col( cind[k] );
} 


void control_pool()
{
    int    k;

    std::cout << "    Controlling pool with heuristic:" << std::endl;
    std::cout << " variables in heuristic:" << std::endl;
    for(k=0;k<nbetter;k++)
        std::cout << " " << better[k]->index;
    std::cout << "\n";
    for(k=0;k<nrows;k++)
        if( violated_by_heur(rows[k]) > 0 ) {
             std::cout << "ERROR: " << k << " invalid constraint"  << std::endl;
             print_row( rows[k] );
        }
}

double violated_by_heur(CONSTRAINT *con)

{
    int      card=0;
    double   *coef=NULL;
    VARIABLE **stack=NULL;
    double   val;


    switch( con->type ){

    case CAPACITY:
        stack = con->stack;
        coef  = con->coef;
        card  = con->card;
        break;
    default:        
        std::cout << "ERROR d: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    val   = -(con->rhs);
    while(card--) if( in_better( stack[card] ) ) val += coef[card];

    switch( con->type ){
    case CAPACITY:
        break;
    default:        
        std::cout << "ERROR e: unknown type "<< con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    if(con->sense=='G') val = -val;
    return(val);
}


int control_constraint(CONSTRAINT *con)

{
    VARIABLE **stack=NULL;
    double   *coef=NULL,*value=NULL;
    int      card=0,i;
    double   c;

    value = (double *)malloc( ncols*sizeof(double) );
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
        //stack = (VARIABLE **)malloc( ncols*sizeof(int) );
        stack = (VARIABLE **)malloc( ncols*sizeof(VARIABLE *) ); /*PWOF*/
        coef = (double *)malloc( ncols*sizeof(double) );
        card = extend_cover(con,stack,coef);
        break;
    default:        
        std::cout << "ERROR f: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    for(i=0;i<ncols;i++)
        value[i] =0;
    for(i=0;i<card;i++)
        value[stack[i]->index]=coef[i];
        
    for(i=0;i<ncols;i++){
        c = get_coeficient(columns+i,con);
        if( fabs(c-value[i])>ZERO ){
            std::cout << "ERROR: val=" << i << "  old_coef=" << value[i] << "  new_coef=" << c << std::endl;
            CSPexit(EXIT_ERROR); //exit(1);
        }
    }
            
    switch( con->type ){
    case COVER:
        free( coef );
        coef = NULL; /*PWOF*/
        free( stack );
        stack = NULL; /*PWOF*/
        break;
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        break;
    default:        
        std::cout << "ERROR 8: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    free(value);
    value = NULL; /*PWOF*/
    return(1);
}


void   control_cut(int nnod,int nedg,int source,int sink,int *fr,int *to,double *weight,double minvalue,int ncut,char   *cutlst)

{
    int  k;
    char s,t;
    double val;

    while(ncut--){
        val = 0.0;
        s = cutlst[source-1];
        t = cutlst[sink-1];
        for(k=0;k<nedg;k++)
            if( cutlst[ fr[k]-1 ] == s && cutlst[ to[k]-1 ] == t )
                val += weight[k];
        if ( val<minvalue-ZERO || val>minvalue+ZERO ){            
            std::cout << " ERROR: in control cut: read=" << val << " virtual=" << minvalue << std::endl;
            CSPexit(EXIT_ERROR); //exit(1);
        }
        cutlst += nnod;
    }
}


static void testing(CONSTRAINT *con)

{
    int      i,j,card;
    VARIABLE **stack;
    double   *coef;
    double   val;
    VARIABLE *list[9];

    list[0] = columns+6;
    list[1] = columns+7;
    list[2] = columns+12;
    list[3] = columns+13;
    list[4] = columns+14;
    list[5] = columns+19;
    list[6] = columns+18;
    list[7] = columns+21;
    list[8] = columns+24;

    stack = con->stack;
    coef  = con->coef;
    card  = con->card;
    val   = -(con->rhs);
    
    for(i=0;i<card;i++)
        for(j=0;j<9;j++)
            if( list[j]==stack[i] ){
                val += coef[i];
                break;
            }
    if( (val<-ZERO && con->sense=='G')||(val>ZERO && con->sense=='L') ){
        std::cout << "\n\nWARNING: violated constraint<<<<<<<<<<<<<<<<<<<<" << std::endl;
    }
    for(j=0;j<9;j++)
        if( list[j]->stat == FIX_LB )
            std::cout << "\n\nWARNING: " << j << " fixed <<<<<<<<<<<<<<<<<<<<" << std::endl;

}

void testing_pool()
{
    int i;
    for(i=0;i<nrows;i++)
        if( rows[i]->type!=BRANCHCUT ) testing( rows[i] );
}
