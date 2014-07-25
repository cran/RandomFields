\name{RMave}
\alias{RMave}
\title{Space-time moving average model}
\description{
 \command{\link{RMave}} is a univariate stationary covariance model
 which depends on a normal scale mixture covariance model \eqn{phi}.
 
 The corresponding covariance function only depends on the difference
 \eqn{(h,u) \in {\bf R}^{d}}{} between two points in
 the \eqn{d}-dimensional space and is given by

 \deqn{C(h, u) = |E + 2Ahh^tA|^{-1/2}
 \phi(\sqrt(\|h\|^2/ 2 + (z^th + u)^2 (1 - 2h^tA (E+2Ahh^tA)^{-1} Ah)))}{
 C(h, u) = |E + 2 A h h^t A|^{-1/2}
 phi(sqrt[|h|^2/2 + (z^t h + u)^2 (1 - 2h^t A (E + 2 A h h^t A)^{-1} A h)])}
 where \eqn{E} is the identity matrix. The spatial dimension is
 \eqn{d-1} and \eqn{h} is real-valued.
}

\usage{
RMave(phi, A, z, spacetime, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a covariance model which is a normal mixture, that means an
   \command{\link{RMmodel}} whose \code{monotone} property equals
   \code{'normal mixture'}, see
   \code{\link{RFgetModelNames}(monotone="normal mixture")}
 }
 \item{A}{a symmetric \eqn{d-1 \times d-1}{d-1 x d-1}-matrix
   if the corresponding random field is in the \eqn{d}{d}-dimensional space}
 \item{z}{a \eqn{d-1} dimensional vector if the corresponding random
   field is on \eqn{d}-dimensional space}
 \item{spacetime}{logical. If FALSE then the model is interpreted as
   if \eqn{h=0}, i.e. the spatial dimension is \eqn{d}. Default is \code{TRUE}}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.}
}
\details{See Schlather, M. (2010), Example 13 with l=1)
}
\value{
 \command{\link{RMave}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Schlather, M. (2010)
 Some covariance models based on normal scale mixtures.
 \emph{Bernoulli}, \bold{16}, 780-797. 
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RFfit}}.
 \command{\link{RFsimulate}},
 \code{\link{RMmodel}},
 \code{\link{RMstp}}
}


\keyword{spatial}
\keyword{models}





\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## Example of an evaluation of the ave2-covariance function
## in three different ways
## ---------------------------------------------------------
## some parameters A and z
A <- matrix(c(2,1,1,2),ncol=2)
z <- c(1,2)
## h for evalutation
h <- c(1,2)
## some abbreviations
E <- matrix(c(1,0,0,1),ncol=2)
B <- A \%*\% h \%*\% t(h) \%*\% A
phi <- function(t){return(RFcov(RMwhittle(1), t))}
## ---------------------------------------------------------
## the following should yield the same value 3 times
## (also for other choices of A,z and h)
z1 <- RFcov( model=RMave(RMwhittle(1),A=A,z=z) , x=t(c(h,0)) )
z2 <- RFcov( model=RMave(RMwhittle(1),A=A,z=z,spacetime=FALSE) , x=t(h) )
z3 <- ( (det(E+2*B))^(-1/2) ) *
 phi( sqrt( sum(h*h)/2 + (t(z) \%*\% h)^2 *
 ( 1-2*t(h) \%*\% A \%*\% solve(E+2*B) \%*\% A \%*\% h) ) )
##
stopifnot(abs(z1-z2)<1e-12, abs(z2-z3)<1e-12)

\dontshow{FinalizeExample()}

}


