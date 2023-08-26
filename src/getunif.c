/* see https://cran.r-project.org/doc/manuals/R-exts.html#Calling-C-from-Fortran-and-vice-versa */

#include <R.h>
#include <Rmath.h>

void F77_SUB(getseed)(void) { GetRNGstate(); }
void F77_SUB(putseed)(void) { PutRNGstate(); }
void F77_SUB(getunif)(double* e) { *e = unif_rand(); }

/* another approach: 
double F77_SUB(getunif)(void) { return unif_rand(); }
*/
