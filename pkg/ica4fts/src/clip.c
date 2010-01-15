/**Clip function - Derivative of the Huber*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

#define absolute(a) (a<0? (-a):(a)) 

SEXP clip(SEXP X, double *c){
	int nrx, ncx;
	int *xdims;
	double *xptr;

	xdims = getDims(X);
	nrx = xdims[0]; ncx = xdims[1];

	PROTECT(X= coerceVector(X, REALSXP));
	xptr = REAL(X);
	double cdoub = 2* *c;
	
	for(int i=0; i< nrx*ncx;i++){
		if(xptr[i] < -*c){
			xptr[i] = -cdoub;			
			}else if(absolute(xptr[i])<=*c){
				xptr[i] = 2*xptr[i];
				}else{
					xptr[i] = cdoub;					
					}			
		}//end for 
		UNPROTECT(1);
		return(X);
}//end clip
