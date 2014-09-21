\name{Distribution Families}
\alias{RR}
\alias{RRmodel}
\alias{RRmodels}
\title{
   Distribution families (RR commands)
}
\description{
 Distribution families to specify random parameters in the model definition.
}


\note{
  The allowance of random parameters is a very recent, developing
  feature of \pkg{RandomField}.
  
  Future changings of the behaviour are not unlikely.
}

\details{
  When simulating Gaussian random fields, the random parameters
  are drawn only once at the very beginning.
  So, if the argument \code{n} in \command{\link{RFsimulate}}
  is greater than \code{1} then \code{n} simulations conditional
  on a single realisation of the random parameters are
  performed. See the examples below.

  There are (simple) multivariate version and additional version to the
  distributions families implemented. Further, \bold{any} distribution
  family defined in R can be used, see the examples below
  
  These function will allow for Baysian modelling. (Future project).
}

\section{Implemented models}{
\tabular{ll}{
  \command{\link{RRdeterm}} \tab no scattering
  \cr
  \command{\link{RRdistr}} \tab  families of distributions transferred
  from R
  \cr
  \command{\link{RRgauss}} \tab a (multivariate) Gaussian random variable
  \cr
  \command{\link{RRloc}} \tab modification of location and scale
  \cr
  \command{\link{RRspheric}} \tab random scale for
  the \command{\link{RMball}} to simulate \command{\link{RMspheric}}, etc.
  %on lower
%  dimensions by \command{\link{RPcoins}} and \command{\link{RPpoisson}}
  \cr
  \command{\link{RRunif}} \tab a (multivariate) uniform random variable
  \cr
}

}

\note{
  A further random element is \command{\link{RMsign}}, which is an
  operator on shape functions. As an exception its name starts with
  \code{RM} and not with \code{RR}.
}

\author{
  Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de/de/publications/software
}
}

\seealso{
  \link{RC},
  \link{RF},
  \link{RM},
  \link{RP},
  \command{\link{Other models}},
  \command{\link{RFdistr}},
 \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}},
 }

\keyword{spatial}
\keyword{distribution}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## here, the scale is given by an exponential variable:
model <- RMgauss(scale=exp())
for (i in 1:4) {
  RFoptions(seed = i)
  # each leads to a simulation with a different scale parameter
  plot(model) ## random
  plot(RFsimulate(model, x=seq(0,10,0.1)))
  readline("press return")
}

# but here, all 4 simulations have same (but random) scale:
plot(RFsimulate(model, x=seq(0,10,0.1), n=4)) 


## hierarchical models are also possible:
## here, the scale is given by an exponential variable whose
## rate is given by a uniform variable
model <- RMgauss(scale=exp(rate=unif()))
plot(model)
plot(RFsimulate(model, x=seq(0,10,0.1)))


## HOWEVER, the next model is deterministic with scale \eqn{e=2.718282}.
model <- RMgauss(scale=exp(1))
plot(model)
plot(RFsimulate(model, x=seq(0,10,0.1)))

\dontshow{FinalizeExample()}

}