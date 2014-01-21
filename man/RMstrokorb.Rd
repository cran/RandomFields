\name{RMstrokorb}
\alias{RMstrokorb}
\alias{RMstrokorbMono}
\alias{RMstrokorbBall}
\alias{RMstrokorbPoly}
\title{Tail correlation function of the Brown-Resnick process}
\description{
  The \command{RMstrokorb} models define various
  shape functions for max-stable processes for a given
  tail correlation function
}
\usage{
RMstrokorbMono(phi);
RMstrokorbBall(phi);
%RMstrokorbPoly(phi);
}
\arguments{
  \item{phi}{a model for a tail correlation function belonging to the
    Gneiting class \eqn{H_d}}
}

\details{
  \command{RMstrokorbMoni} defines
  a monotone shape function that corresponds a tail correlaton
  function belonging to Gneiting's class \eqn{H_d}.

  \command{RMstrokorbBall} defines
  balls with random \emph{radius} that corresponds a tail correlaton
  function belonging to Gneiting's class \eqn{H_d}.

   \command{RMstrokorbPoly} ) defines
   that corresponds a tail correlaton
   function belonging to Gneiting's class \eqn{H_d}.
   It currently only allows for
   \code{\link{RMbrownresnick}(\link{RMfbm}(alpha=1))}.
 
}
\value{
 object of class \code{\link[=RMmodel-class]{RMmodel}}
}

\references{
  \itemize{
    \item Strokorb, K. (2013) \emph{Properties of the Extremal Coefficient
    Functions.} Univ.\ G\"ottingen. PhD thesis. 

   \item
   Strokorb, K., Ballani, F. and Schlather, M. (2014)
   In Preparation.
 }
}
\seealso{
 \command{\link{RFsimulate}},
 \command{\link{RMmodel}}.
}
\examples{
RFoptions(seed=0)
#plot shape function of an M3 process of type 'strokorb'
model <- RMbrownresnick(RMfbm(alpha=1.5, s=0.2))
plot(RMstrokorbMono(model))

x <- seq(0, 10, if (interactive()) 0.005 else 1)
z <- RFsimulate(RPsmith(RMstrokorbMono(model), xi=0), x)
plot(z, type="p", pch=20)

\dontshow{RFoptions(seed=NA)}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}
