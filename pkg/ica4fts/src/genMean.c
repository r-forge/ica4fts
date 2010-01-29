/**Generates the column means of a nxd input matrix. Outputs real vector with dimension d */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

/**Functions takes in nxd matrix x and generates a 1xd matrix with means of each column in x*/
SEXP genMean(SEXP x){
	
	int n, d;
	int *xdims;
	double sum;
	double *ansptr, *xptr;
	SEXP ans; 
	
	xdims = getDims(x);
	n = xdims[0];
	d = xdims[1];
	
	PROTECT(x = coerceVector(x, REALSXP)); //protect 1
	xptr = REAL(x);
	
	PROTECT(ans = allocMatrix(REALSXP,1,d)); //protect 2
	ansptr = REAL(ans);
	
	int k =0; //index for ans
	
	for(int i=0;i<=(d-1)*n;i=i+n){
		sum = 0;
				
		for(int j=i;j<=i+n-1;j++){			
			sum = sum+xptr[j];			
			}//end inner for		
					
			ansptr[k] = sum/n;	//calculate mean and store value
			k++;	
		}//end outer for
	
		UNPROTECT(2);
		return(ans);
		
}//end genMean
