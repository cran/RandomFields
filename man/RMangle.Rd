\name{RMangle}
\alias{RMangle}
\title{Anisotropy matrix given by angle}
\description{
  \code{A = diag(d_1, d_2) \%*\% matrix(ncol=2, c(cos(a), sin(a),
    -sin(a), cos(a)))}
}
 

\usage{
RMangle(angle, ratio, diag) 
}

\arguments{
 \item{angle}{angle \code{a}}
 \item{ratio}{equivalent to \code{diag=c(1, 1/ratio)}}
 \item{diag}{the diagonal components of the matrix}
}

\value{
 \command{\link{RMangle}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}} 
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMtrafo}}
  }

\examples{

RFoptions(seed=0)

x <- seq(0, 2, 0.05)
model <- RMexp(Aniso=RMangle(angle=pi/4, ratio=3))
z <- RFsimulate(x, x, grid=TRUE, model=model)
plot(z)


\dontshow{RFoptions(seed=NA)}

}


\keyword{spatial}
\keyword{models}




