#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

/**Creates an identity matrix indexed by the SEXP format.*/
SEXP identity(int *n){		
	SEXP mat; 
	
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
