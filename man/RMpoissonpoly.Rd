\name{RMpolygon}
\alias{RMpolygon}
\title{RMpolygon}
\description{
  \command{RMpolygon} referes to the indicator function of
  a typical Poisson polygon
}
\usage{
RMpolygon(lambda, safetyfactor)
}
\arguments{
  \item{lambda}{
    intensity of the hyperplan process creating the random shape function

    The default value is \code{1}.
  }
  \item{safetyfactor}{
    the shape function is cut after \code{safetyfactor} times the
    \eqn{d}th root of the area mean or mean volume of the polygon.
    Here, \eqn{d} is the dimension of the field.

    The default value is \code{15}. (The larger the more precise is the
    simulation at a higher prize of computing time.)
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
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}
\examples{
RFoptions(seed=0)
\dontshow{\dontrun{
model <- RPcoins(shape=RMpolygon())
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(RFsimulate(model, x=x))
}}
\dontshow{RFoptions(seed=NA)}
}