/**transProd - takes in matrices x and y and multiplies trans(x) and y */

#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

SEXP transProd(SEXP x, SEXP y){

int nrx, ncx, nry, ncy;
int *xdims, *ydims;
double *ansptr, *xptr, *yptr;

SEXP ans;

xdims = getDims(x); ydims = getDims(y);
nrx = xdims[0]; ncx = xdims[1];
nry = ydims[0]; ncy = ydims[1];

PROTECT(x = coerceVector(x, REALSXP));
PROTECT(y = coerceVector(y, REALSXP));
xptr = REAL(x); yptr = REAL(y);

PROTECT(ans = allocMatrix(REALSXP, ncx,ncy));
ansptr= REAL(ans);

char *transa = "T", *transb = "N";
double one = 1.0; double zero = 0.0;

F77_CALL(dgemm) (transa, transb, &ncx, &ncy, &nrx, &one, xptr, &nrx, yptr, &nry, &zero, ansptr, &ncx);

UNPROTECT(3);
return(ans);

}//end matProd


