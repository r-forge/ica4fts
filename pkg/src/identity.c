/**Create an identity matrix*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/**Creates an identity matrix indexed by the SEXP format. Problem is that we require passing address.*/
SEXP identity(int *n){
	
SEXP mat; 

//Rprintf("Function 'identity' read in dimension is %d \n",*n);

PROTECT(mat = allocMatrix(REALSXP,*n,*n));

for(int i=0;i< *n * *n;i++){
	
	if((i)%(*n+1)==0){
		REAL(mat)[i] = 1;
	}else{
		REAL(mat)[i] = 0;
	}//end if-else		
}//end for

UNPROTECT(1);
return(mat);
}//end identity

