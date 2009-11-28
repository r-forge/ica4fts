/**LogCosh Function*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif


SEXP logCosh(SEXP X, double *c){

int nrx, ncx;
int *xdims;
double *xptr;


xdims = getDims(X);
nrx = xdims[0]; ncx = xdims[1];


PROTECT(X= coerceVector(X, REALSXP));
xptr = REAL(X);
	
	for(int i =0; i< nrx*ncx;i++){		
		xptr[i] = 1/(*c) * log10(cosh(*c * xptr[i]));		
		}//end for 
		
		UNPROTECT(1);
		return(X);
	
}//end logCosh

