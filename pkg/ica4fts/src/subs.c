#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

/**Substitutes identity matrix at i,j coordinates
Inputs: m (dimension) , i (index 1), j (index 2), ag (angle in radians) */
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