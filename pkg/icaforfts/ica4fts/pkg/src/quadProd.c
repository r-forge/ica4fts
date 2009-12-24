/**Quadratic Form multiplication. Inputs vector x and matrix A */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

/**
#ifndef tra
#define tra
#include "transProd.c"
#endif

#ifndef mat
#define mat
#include "matProd.c"
#endif

*/


#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif



/** quadProduct. Takes in two arguments. Argument 1 is nx1 vector x and square postive definite matrix nxn y  
There is no check on pos def of matrix A.*/

SEXP quadProduct(SEXP x, SEXP y){	
	
	SEXP mid;
	SEXP ans;	
			
	PROTECT(mid = transProd(x,y));	
	PROTECT(ans = matProd(mid,x));
	
	UNPROTECT(2);
	return(ans);
}//end quadProd