\name{gradientICA}
\alias{gradientICA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculates the gradient for the objectiveICA function }
\description{
  Uses a series of C functions to calculate the derivative of the objectiveICA function. Offers a choice between the Huber substitution function or the logCosh function.
}
\usage{
gradientICA(T, E, N, C, PH, method = c("Huber", "Cosh"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{ Vector of angles of length p }
  \item{E}{ Financial time series data of dimension nxd }
  \item{N}{ Vector of lags that MUST be arranged in ascending order }
  \item{C}{ Real number C which is used for Huber substitution or Logcosh substitution }
  \item{PH}{ Phi matrix of weights of dimension q*q }
  \item{method}{ Choice between use of Huber or logCosh substitution methods }
}
\details{
  More details to help above
}
\value{
  Gradient of objective value
}
\author{ Erjie Ang \email{ ea75@cornell.edu } }
\seealso{ See }
\examples{
##
d = 2
n = 10000
p = d*(d-1)/2
E = matrix(rnorm(d*n),n,d)
T = as.matrix(rep(0,p))
N = as.matrix(c(1,2))
q = 2*2*p
C = 1
PH = matrix(rnorm(q*q),q,q)
HCans = gradientICA(T,E,N,C,PH,method = "Huber")
LSans = gradientICA(T,E,N,C,PH,method = "Cosh")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ finance }
