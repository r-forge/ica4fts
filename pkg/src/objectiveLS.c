/**Objective functionLS. - Similar to objective_C but uses logcosh function instead. 
Using arbitrary weight matrix P and N, matrix containing arbitrary lags formed according
to rules as stated in crossCovariance_B.c. Phi must be of suitable dimension q which must match N.*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif


/**
#ifndef croC
#include "crossCovariance.c"
#define croC
#endif

#ifndef logC
#include "logCosh.c"
#define logC
#endif

#ifndef qua
#include "quadProd.c"
#define qua
#endif
*/


/**T is px1  matrix with p angles. E is nxd data matrix. N is column vector containing arbitrary lags. C is constant to 
be used in Huber (double scalar), P is arbitrary square Phi matrix which must be of dim q = d(d-1)/2 + L*d(d-1).*/
SEXP objectiveLS(SEXP W, SEXP E, SEXP N, SEXP C, SEXP PH){
	
	//Extract information from E	
	int* Edims;	
	Edims = getDims(E);
	int d = Edims[1];	
	
	//Extract information from C
	double *Cptr;	
	PROTECT(C = coerceVector(C,REALSXP)); //PROTECT 1
	Cptr = REAL(C);	
	double c = Cptr[0]; 		
	
	//Extract information from PH
	int* Pdims;	
	Pdims = getDims(PH);
	int Pc = Pdims[1];	

	//Extract information from N	
	int* Ndims;	
	Ndims = getDims(N);
	int Nr = Ndims[0]; 	
	
	int *nptr;
	PROTECT(N = coerceVector(N,INTSXP)); //PROTECT 2
	nptr = INTEGER(N);	
	int slag = nptr[0]; //smallest lag			

	//Checks
	if(slag ==0){		
		int q = d*(d-1)/2 + (Nr-1)*d*(d-1);
		if(Pc!=q) {Rprintf("Slag = 0; Weight matrix not of right dimension.q = %d, Pc = %d \n",q,Pc);UNPROTECT(2); return R_NilValue; }		
		}else{		
		int q = Nr * d*(d-1);
		if(Pc!=q) {Rprintf("Weight matrix not of right dimension.q = %d, Pc = %d \n",q,Pc); UNPROTECT(2); return R_NilValue; }
		}	
	

	
	SEXP S;
	PROTECT(S = matProd(E,W)); //multiplication of dxd matrix and dxn matrix (transposed from nxd) PROTECT 3
		
	SEXP H; 
	PROTECT(H= logCosh(S,&c)); //Huber transform 		PROTECT 4
	
	SEXP CV; 
	PROTECT(CV= crossCovariance(H,N)); //PROTECT 5
	
 	SEXP ans;
 	PROTECT(ans = quadProduct(CV,PH)); //PROTECT 6	
		
	UNPROTECT(6);
	return ans;	
	
}//end objectiveLS