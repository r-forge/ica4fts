/**Substitutes (i,j) entry of identity matrix to make a rotation matrix with specified input angle*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/**
#ifndef ide
#define ide
#include "identity.c"
#endif

*/


#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

/**Inputs: m (dimension) , i (index 1), j (index 2), ag (angle in radians) */
SEXP subs(int *m, int *i, int *j, double *ag){
	
	int nr=*m;
	double *Mptr;
	
	//Check i and j
	if(*i >= *j){ Rprintf("Input indices are not in right order"); return R_NilValue;}
	
	SEXP M;
	M = identity(&nr);
	PROTECT(M = coerceVector(M,REALSXP));
	Mptr = REAL(M);
		
	//Substitutions
	int ii = (*i-1) * (nr+1); //index of (i,i)
	int jj = (*j-1) * (nr+1); //index of (j,j)
	
	Mptr[ii] = cos(*ag); 
	Mptr[jj] = Mptr[ii]; 
	
	Mptr[ii+(*j-*i)*nr]= -sin(*ag); 
	Mptr[jj-(*j-*i)*nr]= sin(*ag); 
	

	
	UNPROTECT(1);
	return(M);
	
}//end subs

	