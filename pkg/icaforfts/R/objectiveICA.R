objectiveICA <- function(T,E,N,C,PH,method = c("Huber","Cosh")){
	
	n = dim(E)[1]
	d = dim(E)[2]
	p = d*(d-1)/2
	
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
	
	W <-.Call("theta2w",T)
	W= t(W)
	
	if(method=="Huber"){
		CV <-.Call("objectiveHC",W,E,N,C)
		ans <- t(CV)%*%PH%*%CV
		return(ans)
	}else if(method=="Cosh"){
		CV <-.Call("objectiveLS",W,E,N,C)
		ans <- t(CV)%*%PH%*%CV
		return(ans)
	}else{
		stop("Method specification is invalid. Use either Huber or Cosh")	
	}
}#end function