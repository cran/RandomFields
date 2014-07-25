\name{RPprocess}
\alias{RP}
\alias{RPmodel}
\alias{RPmodels}
\alias{RPprocess}
\alias{RPprocesses}
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
  \link{RC}, \link{RR}, \link{RM}, \link{RF}
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\keyword{spatial}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
x <- seq(0, 10, 0.1)
model <- RMexp()

z <- RFsimulate(model, x)
plot(z)

b <- RFsimulate(RPbernoulli(model), x)
plot(b)

sum( abs((z@data$variabl1 >=0 ) - b@data$variable1)) == 0 ## TRUE

\dontshow{FinalizeExample()}
}

