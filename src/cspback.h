int CSPstopcondition(void);
int CSPstoptime(void);
int CSPnewsolution(void);
int CSPexit(int);
int CSPdefinestop(int (*)(void));
int CSPdefineheur(int (*)(void));

#if defined(TAUARGUS_STATIC)
    #define EXPORTFUNC
#elif defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #define EXPORTFUNC __declspec(dllexport)
#else
    #define EXPORTFUNC __attribute__ ((visibility("default")))
#endif


EXPORTFUNC int CSPdefineexit(int (*)(int));
EXPORTFUNC int CSPdefinestoptime(int (*)(void));

#define EXIT_MEMO       1
#define EXIT_ERROR      2
#define EXIT_LPSOLVER   3



