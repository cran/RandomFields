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
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\seealso{
  \command{\link{CovarianceFct}},
  \command{\link{fitvario}}
  \link{sophisticated}
}
\examples{
#####################
## output examples ##
#####################

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



########################################################
## internal estimation of log(scale) instead of scale ##
########################################################

## simulate
x <- c(0, 10, 0.1)
model <- list("$", var=1, scale=1, list("exponential"))
z <- GaussRF(x, model=model, gridtriple=TRUE, n=10)

## usual estimation
est.model <- list("+",
                 list("$", var=NA, scale=NA, list("exponential")),
                 list("$", var=NA, list("nugget")))
fit <- fitvario(x, gridtriple=TRUE, data=z, model=est.model)
str(fit$ml)

## estimating internally log(scale) instead of scale --
## the return is scale itself! However, in this setup
## the transformation leads to worse results
##
parampositions(est.model, dim=1) # indicates where scale is stored;
##                                 namely in the second position:
trafo <- function(param) {param[2] <- exp(param[2]); param}
lower <- list("+",
              list("$", var=NA, scale=-3, list("exponential")),
              list("$", var=NA, list("nugget")))
upper <- list("+",
              list("$", var=NA, scale=3, list("exponential")),
              list("$", var=NA, list("nugget")))
guess <- list("+",
              list("$", var=1, scale=0, list("exponential")),
              list("$", var=1, list("nugget")))
fitlog <- fitvario(x, gridtriple=TRUE, data=z, model=est.model,
                   transform=trafo, lower=lower, upper=upper,
                   users.guess=guess)
str(fitlog$ml) ## note that scale is returned, not log(scale) !



}
\keyword{ spatial }%-- one or more ...
