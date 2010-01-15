gradientICA <- function(T,E,N,C,PH,method = c("Huber","Cosh")){
	n = dim(E)[1]
	d = dim(E)[2]
	p = d*(d-1)/2
	
	#Check types
	if(!is.matrix(T)) stop("Error in argument 1: Vector of angles must be of format matrix.")
	if(!is.matrix(E)) stop("Error in argument 2: Time series data must be of format matrix.")
	if(!is.matrix(N)) stop("Error in argument 3: Vector of lags must be of format matrix.")
	if(!is.real(C)) stop("Error in argument 4: Input must be real nuber.")
	if(!is.matrix(PH)) stop("Error in argument 5: Weight matrix must be of format matrix.")	

	#Check data matrix
	if(d == 0) stop("Data matrix E must have at least 1 column.")

	#Check vector of angles T
	tlength = length(T)
	if(tlength != p) stop("Length of angle input non-comforming to data marix E")
		
	#Check PHI matrix	
	if(N[1]==0){ 
		L = d*(d-1)/2 + (length(N)-1)*d*(d-1)
	}else{ 
		L = length(N)*d*(d-1)
	}
	q = dim(PH)[1]
	if(q!=L) stop("Phi matrix of invalid dimensions")
	

	#Computations
	W <-.Call("theta2w",T)
	W= t(W)
	S <- E%*%W
	P <- t(PH)+PH #transform Phi matrix
	
	if(method=="Huber"){
		ans <-.Call("gradientHC",T,E,S,C,N,P)
		return(ans)
	}else if(method=="Cosh"){
		ans <-.Call("gradientLS",T,E,S,C,N,P)
		return(ans)
	}else{
		stop("Method specification is invalid. Use either Huber or Cosh")	
	}
}#end function