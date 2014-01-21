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

RPbernoulli(phi, loggauss, stationary_only, threshold)

}

\arguments{
  \item{phi}{the \command{\link{RMmodel}}.
    Either a model for a process or a covariance model must be
    specified.
    In the latter case, a Gaussian process \command{\link{RPgauss}} is
    tacitely assumed.}
 \item{loggauss, stationary_only}{optional parameters; same meaning as for 
 \command{\link{RPgauss}}. They will be ignored if the submodel
 is a process definition.}
 \item{threshold}{real valued. \command{\link{RPbernoulli}} returns \eqn{1}
 if value of the random field given by \code{phi} is equal
 to or larger than the value of \code{threshold}, and \eqn{0}
 otherwise.
 
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
 \url{http://ms.math.uni-mannheim.de}
}


\examples{
RFoptions(seed=0)

model <- RPbernoulli(RMexp(), threshold=0)
z <- RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4)
plot(z)

\dontshow{\dontrun{
model <- RPbernoulli(RPbrownresnick(RMexp(), xi=1), threshold=1)
z <- RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4)
plot(z)
}}

\dontshow{RFoptions(seed=NA)}
}


\seealso{
  \command{\link{Auxiliary RMmodels}},
  \link{RP},
 \command{\link{RMbernoulli}}
}




\keyword{spatial}


