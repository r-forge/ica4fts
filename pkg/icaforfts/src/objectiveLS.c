/**Objective functionLS. - Similar to objective_C but uses logcosh function instead. 
Using arbitrary weight matrix P and N, matrix containing arbitrary lags formed according
to rules as stated in crossCovariance_B.c. Phi must be of suitable dimension q which must match N.*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

#ifndef logC
#define logC
/**logCosh substitution function*/
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
#endif


/**T is px1  matrix with p angles. E is nxd data matrix. N is column vector containing arbitrary lags. C is constant to 
be used in Huber (double scalar)*/
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