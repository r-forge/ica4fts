#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

/**Functions takes in nxd matrix x and generates a dxd matrix with means of each column in x*/
SEXP crossMean(SEXP x){
	
	int n, d;
	int *xdims;
	double sum;
	double *meanptr, *xptr;
	SEXP mean; 
	
	xdims = getDims(x);
	n = xdims[0];
	d = xdims[1];
	
	PROTECT(x = coerceVector(x, REALSXP)); //PROTECT 1
	xptr = REAL(x);
	
	PROTECT(mean = allocMatrix(REALSXP,1,d)); //PROTECT 2
	meanptr = REAL(mean);
	
	int k =0; //index for mean
	
	for(int i=0;i<=(d-1)*n;i=i+n){
		sum = 0;				
		
		for(int j=i;j<=i+n-1;j++){			
			sum = sum+xptr[j];			
			}//end inner for					
		
		meanptr[k] = sum/n;	//calculate mean and store value
		k++;	
	}//end outer for
	
		
	/**Mean vector completed. Now generate crossMean vector.*/
	
	SEXP crossM;
	double *crossMptr;
	PROTECT(crossM = allocMatrix(REALSXP,d,d)); //PROTECT 3
	crossMptr = REAL(crossM);
	double prod;
	
	
	int c = 0; //index for crossMean
	
	/**Calculates the mean. First fills the diagonal, then fills the upper and lower triangular portions using symmetry.*/
	for(int i=0; i<d;i++){
		
		c = i*(d+1);
		crossMptr[c] = meanptr[i] * meanptr[i];
		
		for(int j=i+1;j<d;j++){			
			prod = meanptr[i] * meanptr[j];				
			
			c = i + j*d;			
			crossMptr[c] = prod;
			
			c = j + i*d;
			crossMptr[c] = prod; 
			
		}//end inner for
	}//end outer for
	
	
	
	UNPROTECT(3);
	return(crossM);
		
}//end crossMean