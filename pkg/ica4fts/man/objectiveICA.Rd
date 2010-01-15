\name{objectiveICA}
\alias{objectiveICA}
\title{ Calculates the objective function using independent component analysis}
\description{ 
  \code{objectiveICA} performs independent component analysis on financial time series data.
 It first calls C function theta2w to which multiplies a series of independent Given's rotational matrices. Depending on user input for the method, 
either Huber or Cosh, the objective function is calculated using C functions either objectiveHC or objectiveLS 
respectively.
}

\usage{
objectiveICA(T, E, N, C, PH, method = c("Huber", "Cosh"))
}

\arguments{
  \item{T}{ Vector of angles of length p }
  \item{E}{ Financial time series data of dimension nxd }
  \item{N}{ Vector of lags that MUST be arranged in ascending order  }
  \item{C}{ Real number C which is used for Huber substitution or Logcosh substitution }
  \item{PH}{ Phi matrix of weights of dimension q*q }
  \item{method}{ If method is Huber, use Huber substitution method. If method is Cosh, use logcosh substitution method }
}
\details{
  More details to help above
}

\value{
  Objective function value. 
}

\author{ Erjie Ang\email{ ea75@cornell.edu }}

\note{ further notes
}
\seealso{ see also 
}

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
HCans = objectiveICA(T,E,N,C,PH,method = "Huber")
LSans = objectiveICA(T,E,N,C,PH,method = "Cosh")
}
\keyword{ component }
\keyword{ financial }
