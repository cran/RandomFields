\name{RPmodel}
\alias{RPmodel}
\alias{RP}
\title{Models for classes of random fields}
\description{
 Here, all the classes of random fields are described that can be
 simulated
}
\section{Implemented models}{
 \tabular{ll}{
    Gaussian Random fields \tab see \link{Gaussian}\cr
    Max-stable Random Fields \tab see \link{Maxstable}\cr
    Other Random Fields 
    \tab \link[=RPbernoulli]{Binary field} \cr
    \tab \link[=RPchi2]{chi2 field}\cr
    \tab \link[=RPpoisson]{composed Poisson} (shot noise, random coin) \cr
  }
}
\seealso{
  \link{RR}, \link{RM}, \link{RF}
}
\examples{
RFoptions(seed=0)
  x <- seq(0, 10, 0.1)
  model <- RMexp()
  RFoptions(seed=0)
  z <- RFsimulate(model, x)
  plot(z)
  RFoptions(seed=0)
  b <- RFsimulate(RPbernoulli(model), x)
  plot(b)

str(z)

  sum( abs((z@data$variabl1 >=0 ) - b@data$variable1)) == 0 ## TRUE

\dontshow{RFoptions(seed=NA)}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}

\keyword{spatial}

