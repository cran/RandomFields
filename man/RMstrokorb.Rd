\name{Strokorb's Functions}
\alias{RMstrokorb}
\alias{RMstrokorbMono}
\alias{RMstrokorbBall}
\alias{RMstrokorbPoly}
\alias{RMm2r}
\alias{RMm3b}
\alias{RMmps}
\title{Tail correlation function of the Brown-Resnick process}
\description{
  The models define various
  shape functions for max-stable processes for a given
  tail correlation function
}
\usage{
RMm2r(phi)
RMm3b(phi)
RMmps(phi)
}
\arguments{
  \item{phi}{a model for a tail correlation function belonging to the
    Gneiting class \eqn{H_d}}
}

\details{
  \command{RMm2r} used with \command{\link{RPsmith}} defines
  a monotone shape function that corresponds to a tail correlation
  function belonging to Gneiting's class \eqn{H_d}. Currently, the
  function is implemented for dimensions 1 and 3.
  Called as such it returns the corresponding monotone function.
  

  \command{RMm3b} used with \command{\link{RPsmith}} defines
  balls with random \emph{radius} that corresponds to a tail correlation
  function belonging to Gneiting's class \eqn{H_d}. Currently, the
  function is implemented for dimensions 1 and 3.
  (Note that in Strokorb et al. (2014) the density function for twice
  the radius is considered.)
  Called as such it returns the corresponding density function for the
  radius of the balls.

  \command{RMmps}  used with \command{\link{RPsmith}} defines
   random hyperplane polygons
   that corresponds a tail correlaton
   function belonging to Gneiting's class \eqn{H_d}.
   It currently only allows for
   \code{\link{RMbrownresnick}(\link{RMfbm}(alpha=1))} and dimension 2.
   Called as such it returns the tcf defined by the submodel -- this
   definition may change in future.
 
}
\value{
 object of class \code{\link[=RMmodel-class]{RMmodel}}
}

\references{
  \itemize{
    \item Strokorb, K. (2013) \emph{Properties of the Extremal Coefficient
    Functions.} Univ. Goettingen. PhD thesis. 

   \item
   Strokorb, K., Ballani, F. and Schlather, M. (2014)
   In Preparation.
 }
}
\seealso{
 \command{\link{RFsimulate}},
 \command{\link{RMmodel}}.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\keyword{spatial}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMbrownresnick(RMfbm(alpha=1.5, s=0.2))
plot(RMm2r(model))

x <- seq(0, 10, if (interactive()) 0.005 else 1)
z <- RFsimulate(RPsmith(RMm2r(model), xi=0), x)
plot(z, type="p", pch=20)

\dontshow{FinalizeExample()}
}