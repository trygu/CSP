/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPBACK.c                                       */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/


#include <stdlib.h>
#include "cspback.h"
#include "cspmain.h"

static int (*external_stop)(void) = {NULL};
static int (*external_heur)(void) = {NULL};
static int (*external_exit)(int) = {NULL};
static int (*external_time)(void) = {NULL};

/* First error of the current run. Legacy code called exit(1) at fatal
   points; hosts never registered an exit handler, so the flag is the
   real error channel: the optimizer checks it at each node of the
   search and the host maps a nonzero CSPoptimize result to its
   failure codes. First error wins, so the root cause survives. */
static int csp_error = 0;

int CSPerror(void)
{
    return csp_error;
}

void CSPclearerror(void)
{
    csp_error = 0;
}

int CSPstoptime()
{
    if ( external_time==NULL )
        return 0;
    else
        return (*external_time)();
}

int CSPstopcondition()
{
    if( external_stop==NULL )
        return 0;
    else
        return (*external_stop)();
}

int CSPnewsolution()
{
    if( external_heur==NULL )
        return 0;
    else
        return (*external_heur)();
}

int CSPexit(int flag)

{
    if( csp_error==0 )
        csp_error = flag;
    if( external_exit==NULL )
        return 0;
    else
        return (*external_exit)(flag);
}

int CSPdefinestop(int (*stop)(void))

{
    external_stop = stop;
    return 0;
}

int CSPdefineheur(int (*heur)(void))

{
    external_heur = heur;
    return 0;
}

int CSPdefineexit(int (*fexit)(int))

{
    external_exit = fexit;
    return 0;
}

int CSPdefinestoptime(int (*fstoptime)(void))
{
    external_time = fstoptime;
    return 0;
}
