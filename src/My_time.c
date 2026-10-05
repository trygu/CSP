/***********************************************************
 *                                                         *
 *   cputime.c                          J.J.Salazar        *
 *                                                         *
 *   Different functions for getting CPU time since        *
 *   the program start for different computers with        *
 *   different C compilers:                                *
 *                                                         *
 *      BORLAND    for PC with DOS and Borland C           *
 *      WATCOM     for PC with DOS and Watcom C            *
 *      HPUX       for HP with HP-UX                       *
 *      SUN        for SUN with UNIX                       *
 *      ULTRIX     for DEC Station with ULTRIX             *
 *      VAX        for DEC Station with VMS                *
 *                                                         *
 *   float seconds(void)     return the CPU time           *
 *   float elapsed(void)     return the elapsed time       *
 *                                                         *
 ***********************************************************/

#include <time.h>


#ifdef MICROSOFT2
    float seconds(void)
    {
         return( (float)clock()/CLOCKS_PER_SEC );
    }
    float elapsed()
    {
         return( (float)clock()/CLOCKS_PER_SEC );
    }
#endif

#ifdef WATCOM
    float seconds()
    {
         return( (float)clock()/CLOCKS_PER_SEC );
    }
    float elapsed()
    {
         return( (float)clock()/CLOCKS_PER_SEC );
    }
#endif

#ifdef BORLAND
    float seconds()
    {
         return( (float)clock()/CLK_TCK );
    }
    float elapsed()
    {
         return( (float)clock()/CLK_TCK );
    }
#endif

#ifdef HPUX
#include <sys/times.h>
    float   seconds()
    {
        clock_t times();
        struct tms tt;
        times(&tt);
        return( (float)tt.tms_utime/CLK_TCK );
    }
    float   elapsed()
    {
        clock_t times();
        struct tms tt;
        return( (float)times(&tt)/CLK_TCK );
    }
#endif

#ifdef SUN
#include <sys/times.h>
    float   seconds()
    {
        clock_t times();
        struct tms tt;
        times(&tt);
        return( (float)tt.tms_utime/CLK_TCK );
    }
    float   elapsed()
    {
        clock_t times();
        struct tms tt;
        return( (float)times(&tt)/CLK_TCK );
    }
#endif

#ifdef ULTRIX
#include <sys/times.h>
#define    HZ  60.0
    float   seconds()
    {
        clock_t times();
        struct tms tt;
        times(&tt);
        return( (float)tt.tms_utime/HZ );
    }
    float   elapsed()
    {
        clock_t times();
        struct tms tt;
        return( (float)times(&tt)/HZ );
    }
    float   _seconds()
    {
        clock_t times();
        struct tms tt;
        times(&tt);
        return( (float)tt.tms_utime/HZ );
    }
    float   _elapsed()
    {
        clock_t times();
        struct tms tt;
        return( (float)times(&tt)/HZ );
    }
#endif

#ifdef VAX
    float seconds(void)
    {
         return( (float)clock()/CLK_TCK );
    }
    float elapsed()
    {
         return( (float)clock()/CLK_TCK );
    }
#endif

/* Portable fallback for modern compilers (MinGW, macOS, Linux, etc.) where
 * none of the legacy platform macros above is defined. */
#if !defined(MICROSOFT2) && !defined(WATCOM) && !defined(BORLAND) && \
    !defined(HPUX) && !defined(SUN) && !defined(ULTRIX) && !defined(VAX)
    float seconds(void)
    {
         return( (float)clock()/CLOCKS_PER_SEC );
    }
    float elapsed()
    {
         return( (float)clock()/CLOCKS_PER_SEC );
    }
#endif

