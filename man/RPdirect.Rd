\name{Square root}
\alias{Direct}
\alias{RPdirect}
\title{Methods relying on square roots of the covariance matrix}



\description{
 Methods relying on square roots of the covariance matrix
}

\usage{
RPdirect(phi, boxcox) 
}

\arguments{
  \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
    specifies the covariance model to be simulated.}
  \item{boxcox}{the one or two parameters of the box cox transformation.
    If not given, the globally defined parameters are used.
    See \command{\link{RFboxcox}} for details.
  }
 }


\details{
 \command{RPdirect}
  is based on the well-known method for simulating
 any multivariate Gaussian distribution, using the square root of the
 covariance matrix. The method is pretty slow and limited to
 about 12000 points, i.e. a 20x20x20 grid in three dimensions. 
 This implementation can use the Cholesky decomposition and
 the singular value decomposition. 
 It allows for arbitrary points and arbitrary grids.
 }

 \value{
 \command{\link{RPdirect}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}

\references{
 \itemize{
 \item
 Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University. 

 }}
\me

\seealso{ \link{Gaussian},
  \link{RP}, \link{RPsequential}.
}


\keyword{methods}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMgauss(var=10, s=10) + RMnugget(var=0.01)
plot(model, xlim=c(-25, 25))

z <- RFsimulate(model=RPdirect(model), 0:10, 0:10, n=4)
plot(z)

\dontshow{FinalizeExample()}}
