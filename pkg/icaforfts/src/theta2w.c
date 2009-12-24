/**Function theta2w. takes in vector of angles and translates to matrix of dimension n */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

#ifndef iden
#define iden
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
#endif

#ifndef sub
#define sub
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
#endif



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
