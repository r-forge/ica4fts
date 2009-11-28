objectiveICA <- function(T,E,N,C,PH,method = c("H","C")){
	n = dim(E)[1]
	d = dim(E)[2]
	p = d*(d-1)/2

	if(d == 0) stop("Data matrix E must have at least 1 column.")
	
	tlength = length(T)
	if(tlength != p) stop("Length of angle input non-comforming to data marix E")
	
	W <-.Call("theta2w",T)
		W= t(W)
	if(method=="H"){
		.Call("objectiveHC",W,E,N,C,PH)
	}else if(method=="C"){
		.Call("objectiveLS",W,E,N,C,PH)
	}else{
		stop("Method specification is invalid. Use either H or C")	
	}
}#end function


