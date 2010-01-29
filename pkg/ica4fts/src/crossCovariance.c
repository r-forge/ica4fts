#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> //for dgemm

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

/**Cross Mean. Function takes in nxd matrix x and generates a dxd matrix with means of each column in x*/
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



/**Cross Covariance Funciton. W is the data matrix and N is the vector of lags to take cross covariance over*/
SEXP crossCovariance(SEXP W, SEXP N){
	
	//Extract information from N
	int Nrows;
	int *Ndims;
	int *nptr;
	Ndims = getDims(N);
	Nrows = Ndims[0];	
	PROTECT(N = coerceVector(N,INTSXP)); //PROTECT 1
	nptr = INTEGER(N);	
	int slag = nptr[0]; //smallest lag	

	//Extract dimensions of W
	int n, d;
	int *Wdims;
	double *Wptr;	
	Wdims = getDims(W); 
	n = Wdims[0]; 
	d = Wdims[1]; 		
	
	PROTECT(W = coerceVector(W,REALSXP)); //PROTECT 2
	Wptr = REAL(W);
	
	//Generate cross-mean of matrix 
	SEXP mean;
	PROTECT(mean = crossMean(W));//PROTECT 3
	double *mptr = REAL(mean);
	
	int b;	
	//Allocate space for output	
	if(slag ==0){
		b = (d*(d-1))/2 + (Nrows-1)*d*(d-1); //number of entries in output
	}else{
		b = Nrows * d*(d-1);
		}
	double *outptr;
	SEXP out;
	PROTECT(out = allocMatrix(REALSXP,b,1));//PROTECT 4
	outptr = REAL(out);
	
	int z=0; //index for out
	int u = 0;//index for Nptr
	double ans= 0.0; 

	//Only do this if smallest lag is 0
	if(slag ==0){	
	//Lag 0. A total of d(d-1)/2 entries in out will be filled. i and j are column numbers. k is row number 
	for(int i=0;i<d;i=i++){
		for(int j=i+1;j<d;j++){			
			ans = 0.0; //next column so restart													
			for(int k=0;k<n;k++){				
				ans = ans + Wptr[i*n+k]*Wptr[j*n+k]; 								
			} //end k for 									
			outptr[z] = ans/n - mptr[i+j*d];//adjust by cross mean and store in out 			
			z++;							
		}//end j for
	}//end i for  -- END Lag 0
		u++;
	}//end if(slag ==0)
	
	//z contains the index of the next position for outptr
	int lag;
	//Lag 1 and above. 	
	while(u<Nrows){	
		
		lag = nptr[u];
		if(lag >=n){Rprintf("Invalid lag in vector N, lag is greater than # rows in available data. %d \n",lag); break;}
		
	for(int i=0;i<d;i++){		
			int q = i*n +lag; //starting coordinate in 1st matrix
		for(int j=0;j<d;j++){				
				ans=0.0;				
				if(i!=j){				
				 	int r = j*n; //starting coordinate in 2nd matrix				 
				 	for(int k=0;k<n-lag;k++){					 
					 	ans = ans + Wptr[q+k] *Wptr[r+k];			 
						 } //end k for -- Finished down a column					 	 
				 	outptr[z] = ans/n - mptr[i+j*d];				
					 z++;
				}else{ 
						//Do nothing				
					} //end if-else 
								
			}//end j for	-- Finished column in 2nd matrix. Move to next one.		
		}//end i for -- Finished with column in 1st matrix. Move to next one.  
		u++;
} //while loop for index u
	UNPROTECT(4);
	return (out);	
}//end crossCovariance
