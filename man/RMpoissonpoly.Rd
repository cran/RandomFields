\name{RMpolygon}
\alias{RMpolygon}
\title{RMpolygon}
\description{
  \command{RMpolygon} referes to the indicator function of
  a typical Poisson polygon
}
\usage{
RMpolygon(lambda)
}
\arguments{
  \item{lambda}{
    intensity of the hyperplan process creating the random shape function

    The default value is \code{1}.
  }
}



    \seealso{
  \command{\link{RMball}}
  \command{\link{RMspheric}},
 \command{\link{RFsimulate}},
 \command{\link{RMmodel}}.
}

\author{
  Felix Ballani,
  \url{http://www.mathe.tu-freiberg.de/sto/mitarbeiter/felix-ballani}
  
  Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\keyword{spatial}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{\dontrun{
model <- RPcoins(shape=RMpolygon())
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(RFsimulate(model, x=x))
}}
\dontshow{FinalizeExample()}
}