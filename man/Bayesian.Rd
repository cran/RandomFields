\name{Hierarchical Modelling}
\alias{bayesian}
\alias{Bayesian}
\alias{Bayesian Modelling}
\alias{Hierarchical}
\alias{Hierarchical Modelling}
\title{Bayesian Spatial Modelling}
\description{
  \pkg{RandomFields} provides Bayesian modelling to some extend:
  (i) simulation of hierarchical models at arbitrary depth;
  (ii) estimation of the parameteres of a hierarchical model of depth 1
  by means of maximizing the likelihood.
}

\details{
  A Bayesian approach can be taken for scalar, real valued model
  parameters, e.g. the shape parameter \code{nu} in the
  \link{RMmatern} model.
  A random parameter can be passed through a distribution
  of an existing family, e.g. (\code{dnorm}, \code{pnorm}, 
  \code{qnorm}, \code{rnorm}) or self-defined.
  It is passed without the leading letter
  \code{d}, \code{p}, \code{q}, \code{r}, but as a function call
  e.g \code{norm()}.
  This function call may contain arguments that must be
  named, e.g. \code{norm(mean=3, sd=5)}.
  
  Usage:
  \itemize{
    \item \code{exp()} denotes the exponential distribution family
    with rate 1,
    \item \code{exp(3)} is just the scalar \eqn{e^3} and
    \item \code{exp(rate=3)} is the exponential
    distribution family with rate  \eqn{3}.
  }
    
  The family can be passed in three ways:
  \itemize{
    \item implicitelty, e.g. \code{RMwhittle(nu=exp())} or    
    \item explicitely through \command{\link{RRdistr}}, e.g.
    \code{RMwhittle(nu=RRdistr(exp()))}.
    \item by use of \code{\link[=RR]{RRmodels}} of the package
  }
  The first is more convenient, the second more flexible and slightly safer.
  
}

\note{
  \itemize{
    \item
  While simulating any depth of hierarchical modelling is possible,
  estimation is currently restricted to one level of hierarchy.

  \item
  The effect of the distribution family varies between the different processes:

  \itemize{
    \item in Max-stable fields and
    \command{\link{RPpoisson}}, a new realisation of the prior
    distribution(s) is drawn for each shape function
    \item in all the other cases: a realisation of the prior(s)
    is only drawn once.
    This effects, in particular, Gaussian fields with argument
    \code{n>1}, where all the realisations are based on the same
    realisation out of the prior distribution(s).
  }
  
  Note that checking the validity of the
  arguments is rather limited for such complicated models, in general.
  
  }
  
}

%\references{Ribeiro}
\seealso{
  \link{RMmodelsAdvanced}
  For hierarchical modelling see \link{RR}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## See 'RRmodels'for hierarchical models

## the following model defines the argument nu of the Whittle-Matern
## model to be an expontential random variable with rate 5.
model <- ~ 1 + RMwhittle(scale=NA, var=NA, nu=exp(rate=5)) + RMnugget(var=NA)
data(soil)
fit <- RFfit(model, x=soil$x, y=soil$y, data=soil$moisture, modus="sloppy")
print(fit)


\dontshow{FinalizeExample()}

}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}
