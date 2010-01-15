/**Objective functionLS. - Similar to objective_C but uses logcosh function instead. 
Using arbitrary weight matrix P and N, matrix containing arbitrary lags formed according
to rules as stated in crossCovariance_B.c. Phi must be of suitable dimension q which must match N.*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm
#include "matProd.h"
#include "logcosh.h"
#include "crossCovariance.h"

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif


/**T: px1  matrix with p angles. E: nxd data matrix. N: column vector containing arbitrary lags. C: constant to 
be used in logCosh (double scalar)*/
SEXP objectiveLS(SEXP W, SEXP E, SEXP N, SEXP C){
	
	//Extract information from C
	double *Cptr;	
	PROTECT(C = coerceVector(C,REALSXP)); //PROTECT 1
	Cptr = REAL(C);	
	double c = Cptr[0]; 	
	
	SEXP S;
	PROTECT(S = matProd(E,W)); //multiplication of dxd matrix and dxn matrix (transposed from nxd) PROTECT 2
		
	SEXP H; 
	PROTECT(H= logCosh(S,&c)); //Logcosh transform 		PROTECT 3
	
	SEXP CV; 
	PROTECT(CV= crossCovariance(H,N)); //PROTECT 4
		
	UNPROTECT(4);
	return CV;	
	
}//end objectiveLS
