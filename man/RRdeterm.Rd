\name{RRdeterm}
\alias{RRdeterm}
\title{Random scaling used with balls}
\description{
 \command{RRdeterm} refers to the distribution of a deterministic variable.
 }
 \usage{
RRdeterm(mean) 
}

\arguments{
 \item{mean}{the deterministic value}
}

\value{
 \command{\link{RRdeterm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RRdistr}},
\command{\link{RRgauss}},
 }

 

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
x <- seq(-2, 2, 0.001)
p <- RFpdistr(RRdeterm(mean=1), q=x)
plot(x, p, type="l")
\dontshow{FinalizeExample()}
}

 
\keyword{spatial}
\keyword{models}




