\name{Random parameters}
\alias{RR}
\title{
   Random parameters
}
\description{
  Random parameters in the model definition
}


\note{
  The allowance of random parameters is a very recent, developing
  feature of \pkg{RandomField}.
  
  Future changings of the behaviour are not unlikely.
}

\details{
  When simulating Gaussian random fields, the random parameters
  are drawn only once at the very beginning.
  So, if the parameter \code{n} in \command{\link{RFsimulate}}
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
  \command{\link{RRdeterm}} \tab no scattering \cr
  \command{\link{RRdistr}} \tab  families of distributions transferred
  from R
  \cr
  \command{\link{RRgauss}} \tab a (multivariate) Gaussian random variable
  \cr
  \command{\link{RRspherical}} \tab random scale for the
  \command{\link{RMball}} to simulate \command{\link{RMspheric}}, etc.
  %on lower
%  dimensions by \command{\link{RPcoins}} and \command{\link{RPpoisson}}
  \cr
  \command{\link{RRunif}} \tab a (multivariate) uniform random variable
  \cr
}

}



\author{
Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de
}
}

\seealso{
  \command{\link{Other models}},
  \link{RF},
  \command{\link{RFdistr}},
    \link{RM},
 \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}},
  \link{RP},
 \command{\link{RRdistr}},
  \command{\link{RRgauss}},
  \command{\link{RRspherical}}
}

\examples{
RFoptions(seed=0)
## here, the scale is given by an exponential variable:
model <- RMgauss(scale=exp())
for (i in 1:4) {
  # each leads to a simulation with a different scale parameter
  plot(model) ## random !
  plot(RFsimulate(model, x=seq(0,10,0.1)))
  readline("press return")
}

# but here, all 4 simulation have same (but random) scale:
  plot(RFsimulate(model, x=seq(0,10,0.1), n=4)) 

\dontshow{\dontrun{
## hierarchical models are also possible:
## here, the scale is given by an exponential variable whose
## rate is given by a uniform variable
model <- RMgauss(scale=exp(rate=unif()))
plot(model)
plot(RFsimulate(model, x=seq(0,10,0.1)))
}}
\dontshow{RFoptions(seed=NA)}

}
\keyword{spatial}
