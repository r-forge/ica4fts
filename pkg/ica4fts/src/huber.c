#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

#ifndef absolute
#define absolute(a) (a<0? (-a):(a)) 
#endif

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

/**Huber substitution function*/
SEXP Huber(SEXP X, double *c){
	int nrx, ncx;
	int *xdims;
	double *xptr;
	xdims = getDims(X);
	nrx = xdims[0]; ncx = xdims[1];

	PROTECT(X= coerceVector(X, REALSXP));
	xptr = REAL(X);

	double csqr = *c * *c; //calculates constants first
	double cdoub = 2* *c;
	
	for(int i =0; i< nrx*ncx;i++){
			double abx = absolute(xptr[i]);
			if( abx <= *c){
				xptr[i] = xptr[i]* xptr[i];
			}else{				
				xptr[i] = cdoub * abx - csqr;  				
		}//end if-else
	}//end for 
		
		UNPROTECT(1);
		return(X);
}//end Huber