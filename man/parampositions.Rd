\name{parampositions}
\alias{parampositions}
\title{Position of the parameters}
\description{
  The function returns the internal positions of the model parameters 
}
\usage{
  parampositions(model, param, print=TRUE)
}
\arguments{
  \item{model}{see \command{\link{CovarianceFct}}}
  \item{param}{see \command{\link{CovarianceFct}}}
  \item{print}{if \code{FALSE} only an invisible list is returned}
}
\value{
  The model is printed and returned where the values of the parameters
  are the positions in the internal representation.

  An error appears if the model can be substantially simplified.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\seealso{
  \command{\link{CovarianceFct}}, \command{\link{fitvario}}
}
\examples{
## compare the output of the following commands
parampositions(model="exp", param=c(0,1,0,NA))
parampositions(model="exp", param=c(0,1,NA,NA)) 
parampositions(model="exp", param=c(0,0,1,NA)) 
## that is, the nugget in the standard model is removed if naught!
## the values of the other parameters do not matter. (First value
## of the returned vector refers to the mean position.)


parampositions(model="whi", param=cbind(c(1, 1, 1), c(2, 2, 2)))
parampositions(model="whi", param=cbind(c(1, 1, 1), c(2, 0, 2)))
## second lines, second column defines a nugget effect since scale is 0!

try(parampositions(model="whi", param=cbind(c(1, 1, 1), c(0, 0, 0))))
## leads to an error

try(parampositions(model="whi",
                   param=cbind(c(1, 1, 1), c(2, 0, 0), c(1, 0, 0)))) 
## leads to an error

try(parampositions(model="whi",
                   param=cbind(c(1, 1, 1), c(NA, 0, 0), c(1, 0, 0))))
## leads to an error

parampositions(model=list(list(model="exp", var=3, scale=6), "+",
                           list(model="whittle", var=2, scale=7,
                                kappa=NA))) ## again the values do not matter
}
\keyword{ spatial }%-- one or more ...
