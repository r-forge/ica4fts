#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#ifndef getDims
#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol),INTSXP))
#endif

/**Returns W[i,j]*/
double extractElement(SEXP W,int *i,int *j){	
	int *Wdims;
	double *Wptr;
	PROTECT(W= coerceVector(W,REALSXP));
	Wdims = getDims(W);
	int EEnr = Wdims[0]; //rows
	int EEnc = Wdims[1]; //cols
	
	if(*i>EEnr ||*j>EEnc){
		Rprintf("extractElement error: Invalid extraction coordinates.");
		Rprintf("i = %d, j = %d , nr = %d, nc = %d \n",*i,*j,EEnr,EEnc);
		}
	Wptr = REAL(W);		
	UNPROTECT(1);	
	return Wptr[(*i-1)+(*j-1)*EEnr];	
}//end extractElement

/**Ag: vector of angles. E: data matrix. S: ExW transformed data-calculated in R. N: vector of lags. PHI: t(P) + P where P is weight matrix*/
/**Uses Huber and Clip functions*/
SEXP gradientHC(SEXP Ag,SEXP E,SEXP S,SEXP H,SEXP N,SEXP PHI){
	
	//Extract information from E (nxd)
	int *Edims;	Edims = getDims(E);
	int n = Edims[0]; //number of rows of E - n dim
	int d = Edims[1]; //number of cols of E - d dim
	double *eptr; eptr = REAL(E);
	int dist = d*(d-1);	
	int half_dist = dist/2;
		
	//Extract information from H and N
	double h; int slag; double* hptr; int* Nptr; 
	PROTECT(H = coerceVector(H,REALSXP)); hptr = REAL(H);//PROTECT 1
	PROTECT(N = coerceVector(N,INTSXP));  Nptr = INTEGER(N);//PROTECT 2
	h = hptr[0]; slag = Nptr[0]; 
	
	int Nrows = Rf_nrows(N);	
	int p = dist/2;	
	SEXP S2;PROTECT(S2 = Rf_duplicate(S)); //PROTECT 3 - Make a copy

	//Apply huber
	SEXP Hub; PROTECT(Hub = Huber(S,&h)); //PROTECT 4
	double *hubptr; hubptr = REAL(Hub);		

	//Apply clip
	SEXP CP; PROTECT(CP = Clip(S2,&h)); //PROTECT 5
	double *CPptr; CPptr = REAL(CP);
		
	int q;
	if(slag == 0){		
		q = dist/2 + (Nrows-1)*dist;
	}else{
		q = Nrows*dist;			
	}	

	int ell=0; //lag indexing - updated according to k
	
	//Generate constant hjtbar
	SEXP hjtbar;
	PROTECT(hjtbar = genMean(Hub)); //PROTECT 6
	double* HBptr; HBptr = REAL(hjtbar); 
	
	SEXP fbar; PROTECT(fbar = crossCovariance(E,N)); //PROTECT 7	
	SEXP dJdf; PROTECT(dJdf = matProd2(PHI,fbar)); //PROTECT 8
	double* Jfptr; Jfptr = REAL(dJdf);
	
	//dJdT vector of results
	SEXP dJdT; PROTECT(dJdT = allocMatrix(REALSXP,p,1)); //PROTECT 9
	double* JTptr; JTptr = REAL(dJdT);
	
		int a = 1; //index for angles
		int b = 2;			
		
	//Start of outermost loop
	for(int r=1;r<=p;r++){
		
		JTptr[r-1] = 0; //initialize to zero first	
		SEXP dWdT;
		PROTECT(dWdT = dwdtheta(Ag,&a,&b)); //dWdT matrix with derivative of matrix(a,b)
			
		int I=1;
		int J=2;
		
	for(int k =1;k<=q;k++){
		
		//adjust ell
		if(slag==0){
				if(k<=dist/2){
					ell=0;			
				}else{				
					int u = k - dist/2 - 1;
					u = u/dist + 1; 
					ell = Nptr[u];									
				}//end if-else
			}else{ //slag !=0						
				int u = k-1;
				u = u/dist; 
				ell = Nptr[u];			
				}		
					
		double dfdT =0;

	for(int j=1;j<=d;j++){	
			for(int i=1;i<=d;i++){				
			
				//Storage for g(I,j) - Fixed for each I and j
				SEXP G; G = allocMatrix(REALSXP,n,1);  //g(I,j)
				double *gptr; gptr = REAL(G);
				double sum=0; double gbar = 0;
			
				for(int p=0;p<n;p++){				
					gptr[p] = CPptr[(j-1)*n+p]* eptr[(i-1)*n+p];
					sum = sum+gptr[p];				
					}//end for p 								
				gbar = sum/(double)n;				

				double dfdW = 0;				
				double temp = 0;			
				
				if(j==I){	
												
					for(int t=ell+1;t<=n;t++){					
						temp = temp + hubptr[(J-1)*n + (t-ell-1)] * gptr[t-1];											
					}//end for t					
					
					temp = temp/(double)n;					
					double sec = HBptr[J-1]*gbar;						
					dfdW = dfdW + temp - sec;									
				}//end if(j==I)
						
				
				
				if(j==J){					
										
					for(int t=ell+1;t<=n;t++){	
						temp = temp + hubptr[(I-1)*n + t-1]* gptr[t-ell-1];													
					}//end for t		
					
					temp = temp/(double) n;
					double fourth = HBptr[I-1]*gbar;
				    dfdW = dfdW + temp - fourth;
				} //end if(j==J)	
		
		dfdT = dfdT + dfdW*extractElement(dWdT,&j,&i);	
			}//end for i
		}//end for j 
	
		JTptr[r-1] = Jfptr[k-1]*dfdT + JTptr[r-1]; 	
		
		//Adjustments for I and J
		if(k< half_dist){
				if(J<d){
					J++;
				}else{
					I++;
					J=I+1;
					}		
			}else if(k == half_dist|| (k-half_dist)%dist ==0){				
				I=1;J=2; //reset
					
			}else{				
				if(J<d){
					J++;
					}else{
					I++;J=1;	
						}				
				if(I==J) J++; //skip		
			}	
	}//end for k 
	
		//indexing for dWdT
		if(a==1){
			a=b;
			b=b+1;		
		}else{
			a=a-1;	
		}
	UNPROTECT(1); //Remove dWdT from memory stack. Create a new one in the next iteration.		
	}//end for r		
	
	UNPROTECT(9);
	return dJdT;	

}//end gradientHC