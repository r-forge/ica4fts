/**Function theta2w. takes in vector of angles and translates to matrix of dimension n */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif




/**

#ifndef sub
#define sub
#include "substitute.c"
#endif

#ifndef mat
#define mat
#include "matProd.c"
#endif

*/

/**Assumed that vector W is of a valid length such that # of elements = # of elements in Up-Trig portion of dxd matrix*/
SEXP theta2w(SEXP W){	
	
	int p;
	int *Wdims;
	double *Wptr;
	
	Wdims = getDims(W); //extract dimensions of W
	p = Wdims[0]; //number of entries
	PROTECT(W = coerceVector(W,REALSXP));
	Wptr = REAL(W);	
	
	int d = (sqrt(8*p+1) + 1)/2;	
		
	SEXP ans = identity(&d); //initialize
	SEXP interm; //intermediate matrix 
	int k = 0; //index for W
	
	for(int j=2;j<=d;j++){
		for(int i=j-1;i>=1;i--){			
			
			interm = subs(&d,&i,&j,&Wptr[k]);
			ans = matProd(interm,ans);			
	
			k++;
			
			}//end inner for		
		}//end outer for 		
	
		
		
		UNPROTECT(1);		
		return ans;
		//return ans;	
}//end function theta2w