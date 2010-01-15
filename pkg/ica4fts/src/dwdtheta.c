/**dwdtheta function: takes the derivative of identified rotation matrix, specified by inputs a and b.*/
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif


/**Creates an identity matrix indexed by the SEXP format. Problem is that we require passing address.*/
SEXP zeros(int *n){
	
	SEXP mat; 
	PROTECT(mat = allocMatrix(REALSXP,*n,*n));

	for(int i=0;i< *n * *n;i++){		
			REAL(mat)[i] = 0;
	}//end for

	UNPROTECT(1);
	return(mat);
}//end zeros


/**Inputs: m (dimension) , i (index 1), j (index 2), ag (angle in radians) */
SEXP derivative(int *m, int *i, int *j, double *ag){
	
	int nr=*m;
	double *Mptr;
	
	//Check i and j
	if(*i >= *j){ Rprintf("Input indices are not in right order"); return R_NilValue;}
	
	SEXP M;
	M = zeros(&nr);
	PROTECT(M = coerceVector(M,REALSXP));
	Mptr = REAL(M);
		
	//Substitutions
	int ii = (*i-1) * (nr+1); //index of (i,i)
	int jj = (*j-1) * (nr+1); //index of (j,j)
	
	Mptr[ii] = -sin(*ag); 
	Mptr[jj] = Mptr[ii];	
	Mptr[ii+(*j-*i)*nr]= -cos(*ag); 
	Mptr[jj-(*j-*i)*nr]= cos(*ag);  
	
	UNPROTECT(1);
	return(M);	
}//end derivative




/**W is the vector of angles and a,b specifies which matrix to take the derivative of. Result is dxd matrix - product
	of rotational matrices and one derivative of a rotational matrix.*/
SEXP dwdtheta(SEXP W,int *a,int *b){

	int p;
	int *Wdims;
	double *Wptr;
	
	Wdims = getDims(W); //extract dimensions of W
	p = Wdims[0]; //number of entries
	PROTECT(W = coerceVector(W,REALSXP));
	Wptr = REAL(W);	
	
	int d = (sqrt(8*p+1) + 1)/2;	
	
	int A = *a;
	int B = *b;
	
	SEXP ans = identity(&d); //initialize
	SEXP interm; //intermediate matrix 
	int k = 0; //index for W

	for(int j=2;j<=d;j++){
		for(int i=j-1;i>=1;i--){			
			
			if(i==A && j ==B){
				interm = derivative(&d,&i,&j,&Wptr[k]);	
			}else{
			interm = subs(&d,&i,&j,&Wptr[k]);			
			}
			ans = matProd(interm,ans);
			
			k++;
			}//end inner for		
		}//end outer for 
	UNPROTECT(1);
	return ans;

}//end dwdtheta
