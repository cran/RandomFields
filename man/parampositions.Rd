\name{parampositions}
\alias{parampositions}
\title{Position of the parameters}
\description{
  The function returns the internal positions of the model parameters 
}
\usage{
  parampositions(model, param, trend=NULL, dim, print=1)
}
\arguments{
  \item{model}{see \command{\link{CovarianceFct}}}
  \item{param}{see \command{\link{CovarianceFct}}}
  \item{trend}{trend}
  \item{dim}{dimension of the field}
  \item{print}{if \code{0} only an invisible list is returned}
}
\value{
  The model is printed and returned where the values of the parameters
  are the positions in the internal representation.

  An error appears if the model can be substantially simplified.
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\seealso{
  \command{\link{CovarianceFct}},
  \command{\link{fitvario}}
  \link{sophisticated}
}
\examples{
## basic models
parampositions(model="exp", param=c(1,2,3,NA), dim=1, print=3) 
parampositions(model="exp", param=c(0,1,NA,NA), dim=1, print=3) 

## nested 
parampositions(model="whi", param=rbind(c(1, NA, 3), c(4, 5, NA)),
               dim=1, print=3)

## complicated models
model <- list("+",
              list("$", aniso=matrix(c(5:7, NA), nc=2), var=NA, list("exp")), 
              list("$", var=2, scale=7, list("whittle", nu=NA)))
parampositions(model=model, dim=2, print=3) 
}
\keyword{ spatial }%-- one or more ...
