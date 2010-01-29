/**Matrix Product of Matrices 2*x and y*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

SEXP matProd2(SEXP x, SEXP y){

int nrx, ncx, nry, ncy;
int *xdims, *ydims;
double *ansptr, *xptr, *yptr;

SEXP ans;

xdims = getDims(x); ydims = getDims(y);
nrx = xdims[0]; ncx = xdims[1];
nry = ydims[0]; ncy = ydims[1];

if(ncx !=nry){Rprintf("Matrix Product Error: matrix Dimensions unmatched! ncx = %d, nry = %d\n",ncx,nry);}

PROTECT(x = coerceVector(x, REALSXP));
PROTECT(y = coerceVector(y, REALSXP));
xptr = REAL(x); yptr = REAL(y);

PROTECT(ans = allocMatrix(REALSXP, nrx,ncy));
ansptr= REAL(ans);

char *transa = "N", *transb = "N";
double two = 2.0; double zero = 0.0;

F77_CALL(dgemm) (transa, transb, &nrx, &ncy, &ncx, &two, xptr, &nrx, yptr, &nry, &zero, ansptr, &nrx);

UNPROTECT(3);
return(ans);

}//end matProd
