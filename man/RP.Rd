\name{RPprocess}
\alias{RP}
\alias{RPmodel}
\alias{RPmodels}
\alias{RPprocess}
\alias{RPprocesses}
\title{Models for classes of random fields (RP commands)}
\description{
 Here, all classes of random fields are described that can be
 simulated.
}
\section{Implemented processes}{
 \tabular{ll}{
    Gaussian Random Fields \tab see \link{Gaussian}\cr
    Max-stable Random Fields \tab see \link{Maxstable}\cr
    Other Random Fields 
    \tab \link[=RPbernoulli]{Binary field} \cr
    \tab \link[=RPchi2]{chi2 field}\cr
    \tab \link[=RPpoisson]{composed Poisson} (shot noise, random coin) \cr
   \tab \link[=RPt]{t field}\cr
 }
}
\seealso{
  \link{RC}, \link{RR}, \link{RM}, \link{RF}, \link{R.}
}

\me

\keyword{spatial}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
x <- seq(0, 10, 0.1)
model <- RMexp()

## a Gaussian field with exponential covariance function
z <- RFsimulate(model, x)
plot(z)

## a binary field obtained as a thresholded Gaussian field
b <- RFsimulate(RPbernoulli(model), x)
plot(b)

sum( abs((z@data$variabl1 >=0 ) - b@data$variable1)) == 0 ## TRUE,
## i.e. RPbernoulli is indeed a thresholded Gaussian process

\dontshow{FinalizeExample()}}

