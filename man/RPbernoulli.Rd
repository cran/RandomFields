\name{RPbernoulli}
\alias{RPbernoulli}

\title{Simulation of Binary Random Fields}
\description{
 Indicator or binary field which
 has the value 1, if an underfield field
 exceeds a given threshold, 0
 otherwise. 
}

\usage{

RPbernoulli(phi, stationary_only, threshold)

}

\arguments{
  \item{phi}{the \command{\link{RMmodel}}.
    Either a model for a process or a covariance model must be
    specified.
    In the latter case, a Gaussian process \command{\link{RPgauss}} is
    tacitely assumed.}
  \item{stationary_only}{optional arguments; same meaning as for 
 \command{\link{RPgauss}}. It is ignored if the submodel
 is a process definition.}

 \item{threshold}{real valued.
   \command{\link{RPbernoulli}} returns \eqn{1}
   if value of the random field given by \code{phi} is equal
   to or larger than the value of \code{threshold}, and \eqn{0}
   otherwise. In the multivariate case, a vector might be given.
   If the threshold is not finite, then the original field is returned.
 
   \code{threshold} defaults value is 0.
 }
}



\value{
 The function returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}

\details{
 \command{\link{RPbernoulli}} can applied to any field. If only
 a covariance model is given, a Gaussian field is simulated as
 underlying field. 
}

 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\seealso{
  \command{\link{Auxiliary RMmodels}},
  \link{RP},
 \command{\link{RMbernoulli}}
}


\keyword{spatial}



\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
x <- seq(0, 10, 0.1)
model <- RPbernoulli(RMexp(), threshold=0)
z <- RFsimulate(model, x, x, n=4)
plot(z)

\dontshow{\dontrun{ % to do sobald br wieder laueft
model <- RPbernoulli(RPbrownresnick(RMexp(), xi=1), threshold=1)
z <- RFsimulate(model, x, x, n=4)
plot(z)
}}
\dontshow{FinalizeExample()}
}
