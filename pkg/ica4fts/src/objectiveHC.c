#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

#ifndef absolute
#define absolute(a) (a<0? (-a):(a)) 
#endif


/**T: px1  matrix with p angles. E: nxd data matrix. N: column vector containing arbitrary lags. C: constant to 
be used in Huber (double scalar)*/
SEXP objectiveHC(SEXP W, SEXP E, SEXP N, SEXP C){	

	//Extract information from C
	double *Cptr;	
	PROTECT(C = coerceVector(C,REALSXP)); //PROTECT 1
	Cptr = REAL(C);	
	double c = Cptr[0]; 		
	
	SEXP S;
	PROTECT(S = matProd(E,W)); //multiplication of dxd matrix and dxn matrix (transposed from nxd) PROTECT 2
		
	SEXP H; 
	PROTECT(H= Huber(S,&c)); //Huber transform 		PROTECT 3
	
	SEXP CV; 
	PROTECT(CV= crossCovariance(H,N)); //PROTECT 4
		
	UNPROTECT(4);
	return CV;	
	
}//end objectiveHC