#ifdef VSCIP
#include <highs/interfaces/highs_c_api.h>
#endif

#if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #define EXPORTFUNC __declspec(dllexport)
#else
    #define EXPORTFUNC __attribute__ ((visibility("default")))
#endif

#ifdef CPLEX7
#include <cplex.h>
#endif

#ifdef XPRESS_13
#include "xprs.h"
#endif

#include <iostream>
#include <fstream>

// Header files from tauhitas for callback facility with Java
#include "IProgressListener.h"
#include "ICallback.h"

// Different implementation for different solvers
// Assigned to different namespaces to be able to choose solver at runtime
#ifdef CPLEX7
namespace CPLEXv
{
extern CPXENVptr Env;
EXPORTFUNC void CSPSetFileNames(const char*);    
EXPORTFUNC void CSPFreeFileNames();
EXPORTFUNC void CSPSetDoubleConstant(const int, double);
EXPORTFUNC double CSPGetDoubleConstant(const int);
EXPORTFUNC void CSPSetIntegerConstant(const int, int);
EXPORTFUNC int CSPGetIntegerConstant(const int);
EXPORTFUNC int CSPloadprob(int,double*,int,double*,int*,char*,double*,double*,double*,double*,char**,int*,int*,signed char*);
EXPORTFUNC int CSPoptimize(IProgressListener*);
EXPORTFUNC int CSPfreeprob();
EXPORTFUNC int CSPsolution(int*,int*,char*);
EXPORTFUNC int CSPrelbounds(int,int*,double*,double*,char);
}  //namespace CPLEXv end
#endif

#ifdef VSCIP
namespace SCIPv
{
extern void *_scip;
EXPORTFUNC void CSPSetFileNames(const char*);    
EXPORTFUNC void CSPFreeFileNames();
EXPORTFUNC void CSPSetDoubleConstant(const int, double);
EXPORTFUNC double CSPGetDoubleConstant(const int);
EXPORTFUNC void CSPSetIntegerConstant(const int, int);
EXPORTFUNC int CSPGetIntegerConstant(const int);
EXPORTFUNC int CSPloadprob(int,double*,int,double*,int*,char*,double*,double*,double*,double*,char**,int*,int*,signed char*);
EXPORTFUNC int CSPoptimize(IProgressListener*);
EXPORTFUNC int CSPfreeprob();
EXPORTFUNC int CSPsolution(int*,int*,char*);
EXPORTFUNC int CSPrelbounds(int,int*,double*,double*,char);
}  //namespace SCIPv end
#endif

#ifdef XPRESS_13
namespace XPRESSv
{
EXPORTFUNC void CSPSetFileNames(const char*);    
EXPORTFUNC void CSPFreeFileNames();
EXPORTFUNC void CSPSetDoubleConstant(const int, double);
EXPORTFUNC double CSPGetDoubleConstant(const int);
EXPORTFUNC void CSPSetIntegerConstant(const int, int);
EXPORTFUNC int CSPGetIntegerConstant(const int);
EXPORTFUNC int  CSPloadprob(int,double*,int,double*,int*,char*,double*,double*,double*,double*,char**,int*,int*,signed char*);
EXPORTFUNC int  CSPoptimize(IProgressListener*);
EXPORTFUNC int  CSPfreeprob();
EXPORTFUNC int  CSPsolution(int*,int*,char*);
EXPORTFUNC int  CSPrelbounds(int,int*,double*,double*,char);
}  //namespace SCIPv end
#endif

// Functions to be used inside own code
int CSPtestprob(int,double*,int,double*,int*,char*,double*,double*,double*,double*,char**,int*,int*,signed char*);
int CSPabsbounds(int,int*,double*,double*,char);
int CSPpartialbounds();
int CSPwrite(char*);