#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(asm18_r)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dupmean_r)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(eha_rv6)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mwin_r)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mwincenter_r)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(peak_r)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(peakfilter_r)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(trough_r)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(tune_r)(void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"asm18_r",      (DL_FUNC) &F77_NAME(asm18_r),      20},
    {"dupmean_r",    (DL_FUNC) &F77_NAME(dupmean_r),     6},
    {"eha_rv6",      (DL_FUNC) &F77_NAME(eha_rv6),      22},
    {"mwin_r",       (DL_FUNC) &F77_NAME(mwin_r),        9},
    {"mwincenter_r", (DL_FUNC) &F77_NAME(mwincenter_r),  9},
    {"peak_r",       (DL_FUNC) &F77_NAME(peak_r),        6},
    {"peakfilter_r", (DL_FUNC) &F77_NAME(peakfilter_r), 12},
    {"trough_r",     (DL_FUNC) &F77_NAME(trough_r),      6},
    {"tune_r",       (DL_FUNC) &F77_NAME(tune_r),        6},
    {NULL, NULL, 0}
};

void R_init_astrochron(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
