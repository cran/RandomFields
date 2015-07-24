\name{RMtbm}
\alias{RMtbm}
\title{Turning Bands Method}
\description{
 \command{\link{RMtbm}} is a univariate stationary isotropic covariance
 model in dimension \code{reduceddim} which depends on a univariate stationary
 isotropic covariance \eqn{\phi}{phi} in a bigger dimension \code{fulldim}.
 For formulas for the covariance function see details.
}
\usage{
RMtbm(phi, fulldim, reduceddim, layers, var, scale, Aniso, proj)
}
\arguments{ 
 \item{phi, fulldim, reduceddim, layers}{see \command{\link{RPtbm}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 \command{\link{RMtbm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\details{
  The turning bands method stems from the 1:1 correspondence between the
  isotropic covariance functions of different dimensions. See Gneiting
  (1999) and Strokorb and Schlather (2014).
  
  The standard case \code{reduceddim=1} and \code{fulldim=3}.
  If only one of the arguments are given, then the difference of the two
  arguments equals 2.
  
  For \code{d == n + 2}, where \code{n=reduceddim} and 
  \code{d==fulldim} the original dimension, we have
  
  \deqn{
    C(r) = \phi(r) + r \phi'(r) / n
  }{
    C(r) = phi(r) + r phi'(r) / n
  }
  which, for \code{n=1} reduced to the standard TBM operator
  \deqn{
    C(r) =\frac {d}{d r} r \phi(r)
  }{
    C(r) = d/dr [ r phi(r) ]
  }
  
  For \code{d == 2 && n == 1} we have
  
  \deqn{
    C(r) = \frac{d}{dr}\int_0^r \frac{u\phi(u)}{\sqrt{r^2 - u^2}} d u
  }{
    C(r) = d/dr int_0^r [ r phi(r) ] / [ sqrt{r^2 - u^2} ] d u
  }

  \sQuote{Turning layers} is a generalization of the turning bands
  method, see Schlather (2011).
 
}
\references{
  Turning bands
  \itemize{
    \item Gneiting, T. (1999)
    On the derivatives of radial positive definite function.
    \emph{J. Math. Anal. Appl}, \bold{236}, 86-99
    
    \item
    Matheron, G. (1973).
    The intrinsic random functions and their applications.
    \emph{Adv . Appl. Probab.}, \bold{5}, 439-468.
 
    \item
    Strokorb, K., Ballani, F., and  Schlather, M. (2014)
    Tail correlation functions of max-stable processes: Construction
    principles, recovery and diversity of some mixing max-stable
    processes with identical TCF. 
    \emph{Extremes}, \bold{} Submitted.    
  }

  Turning layers
  \itemize{
    \item Schlather, M. (2011) Construction of covariance functions and
    unconditional simulation of random fields. In Porcu, E., Montero, J.M.
    and Schlather, M., \emph{Space-Time Processes and Challenges Related
      to Environmental Problems.} New York: Springer.
  }

}
\seealso{
 \command{\link{RPtbm}},
 \command{\link{RFsimulate}}.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}
\keyword{models}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
x <- seq(0, 10, 0.02)
model <- RMspheric()
plot(model, model.on.the.line=RMtbm(RMspheric()), xlim=c(-1.5, 1.5))

z <- RFsimulate(RPtbm(model), x, x)
plot(z)
\dontshow{FinalizeExample()}
}
