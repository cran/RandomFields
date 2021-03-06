\name{Extremal t}
\alias{RPopitz}
\alias{extremal t}
\alias{extremal t process}

\title{Extremal t process}
\description{ 
  \command{RPopitz} defines an extremal t process.
}

\usage{
RPopitz(phi, xi, mu, s, alpha)
}

\arguments{
  \item{phi}{an \command{\link{RMmodel}};
    covariance model for a standardized
    Gaussian random field, or the field itself.
  }
  \item{xi,mu,s}{the extreme value index, the location parameter and the
    scale parameter, respectively, of the generalized extreme value
    distribution. See Details.
  }
  \item{alpha}{originally referred to the \eqn{\alpha}-Frechet marginal
    distribution, see the original literature for details.
  }
}

\details{\GEV
}
\me

\references{
  \itemize{
    \item
    Davison, A.C., Padoan, S., Ribatet, M. (2012).
    Statistical modelling of spatial extremes.
    \emph{Stat. Science} \bold{27}, 161-186.
   \item
   Opitz, T. (2012) A spectral construction of the extremal t process.
    \emph{arxiv} \bold{1207.2296}.
  }
}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RPgauss}},
 \command{\link{maxstable}},
 \command{\link{maxstableAdvanced}}.
}

\keyword{spatial}


\examples{\dontshow{StartExample()}

## sorry, does not work savely yet

\dontshow{\dontrun{
RFoptions(seed=0, xi=0)
## seed=0: *ANY* simulation will have the random xseed 0; set
##         RFoptions(seed=NA) to make them all random again
## xi=0: any simulated max-stable random field has extreme value index 0
x <- seq(0, 2, 0.01)
model <- RPopitz(RMgauss(), alpha=2)
z1 <- RFsimulate(model, x)
plot(z1, type="l")

}}

\dontshow{FinalizeExample()}}

