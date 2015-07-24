\name{RMmodelsSpacetime}
\alias{RMmodelsSpacetime}
\alias{RMmodelsSpaceTime}
\alias{space-time}
\title{Space-time Covariance Models}
\description{
  Here, a collection of implemented space-time
  models is given.
}

\details{

\bold{Stationary space-time models}

Here, most of the models are composed models (operators).
Note that in space-time modelling the argument \code{proj} may
take also the values \code{"space"} for the projection on the space
and \code{"time"} for the projection onto the time axis.

\tabular{ll}{
  separable models \tab are easily constructed using
  \command{\link[=RMplus]{+}},   \command{\link[=RMmult]{*}},
  and \link[=RMmodels]{proj}, see also the example below
  \cr
  \command{\link{RMave}} \tab space-time moving average model
  \cr
  \command{\link{RMcoxisham}} \tab Cox-Isham model
  \cr
  \command{\link{RMcurlfree}} \tab curlfree (spatial) field (stationary
  and anisotropic)
  \cr
  \command{\link{RMdivfree}} \tab divergence free (spatial) vector valued
  field, (stationary and anisotropic)
  \cr
  \command{\link{RMgennsst}} \tab generalization of Gneiting's
  non-separable space-time model
  \cr 
  \command{\link{RMiaco}} \tab non-separabel space-time model
  \cr
  % obsolete -- included by Cauchy --- should be given by an example
  \command{\link{RMmastein}} \tab Ma-Stein model
  \cr
  \command{\link{RMnsst}} \tab Gneiting's non-separable space-time model
  \cr 
  \command{\link{RMstein}} \tab Stein's non-separabel space-time model
  \cr
  \command{\link{RMstp}} \tab Single temporal process
  \cr
  \command{\link{RMtbm}} \tab Turning bands operator
}
}


\references{
 \itemize{
   \item Schlather, M. (2011) Construction of covariance functions and
   unconditional simulation of random fields. In Porcu, E., Montero, J.M.
   and Schlather, M., \emph{Space-Time Processes and Challenges Related
     to Environmental Problems.} New York: Springer.
   }
}

\seealso{\command{\link{RFformula}},
  \link{RM},
  \command{\link{RMmodels}},
  \command{\link{RMmodelsAdvanced}}
}

\author{
 Alexander Malinowski, \email{malinowski@math.uni-mannheim.de}
 
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
## separable model with expontential model in space and gaussian in time
model <- RMexp(proj = "space") * RMgauss(proj = "time")
xT <- seq(0, 10, 0.1)
z <- RFsimulate(model, x=xT, T=xT)
plot(z)

\dontshow{FinalizeExample()}
}
