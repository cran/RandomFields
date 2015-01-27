\name{RFratiotest}
\alias{RFratiotest}
\alias{print.RFratiotest}
\title{Likelihood ratio test}
\description{
  The function performs an approximate $\eqn{chi^2} test or
  a Monte Carlo likelihood ratio test
  based on \command{\link{fitgauss}}. Currently it only works for
  Gaussian random fields.
}    
\usage{
RFratiotest(nullmodel, alternative, x, y = NULL, z = NULL, T = NULL,
grid=NULL, data,
alpha, n = 5 / alpha, seed = 0, 
lower = NULL, upper = NULL, bc_lambda, methods,
sub.methods, optim.control = NULL, users.guess = NULL,
distances = NULL, dim, transform = NULL, ...)
}
\arguments{  
  \item{nullmodel, alternative}{See Details.
     
    The set of parameters 
    to be estimated for \code{nullmodel}
    should be a subset of the parameters to be estimated for
    \code{alternative} if \code{alternative} is given.
  }
  \item{alpha}{value in [0,1] or missing. Significance level.
  }
 \item{n}{integer. The test is based on \code{n-1} simulations.}
 \item{seed}{integer. If not \code{NULL} and not \code{NA}, the
   \link[base]{.Random.seed} is set to \code{seed}. Otherwise,
   \command{\link[base]{set.seed}} is set to the value of 
   \code{RFoptions{}$general$seed}
   if the latter is not \code{NA}.
 }
 \item{x, y, z, T, grid, data, lower, upper, bc_lambda, methods,
   sub.methods, optim.control, users.guess, distances, dim, transform,
   ...}{see \command{\link{RFfit}}}
}

\note{
  An important \command{\link{RFoptions}} is \code{ratiotest_approx}.
}

\details{
  \code{nullmodel} (and the \code{alternative}) can be
  \itemize{
    \item
    a covariance model, see \command{\link{RMmodel}} or
    type \code{\link{RFgetModelNames}(type="variogram")} to get all options.

    Depending weather the \command{\link{RFoptions}}
    \code{ratiotest_approx} is TRUE the
    the chisq approximation is performed. Otherwise 
    a Monte Carlo ratio test is performed.
    \item
    \code{\link[=RFfit-class]{RFfit}} or
    \code{\link[=RMmodelFit-class]{RMmodelFit}}

    Here, a chisq approximative test is always performed on the
    already fitted models.
  }
  
  \command{RFratiotest} tries to detect whether
  \code{nullmodel} is a submodel of \code{alternative}.
  If it fails,
  \itemize{
    \item a message is printed that says that an
    \emph{automatic} detection has not been possible;
    \item it is not guaranteed anymore that the \code{alternative} model
    returns a (log) likelihood that is at least as large as
    that of the \code{nullmodel},
    even if \code{nullmodel} is a submodel of \code{alternative}.
    This is due to numerical optimisation which is never perfect.       
  }
  Otherwise it is guaranteed that the \code{alternative} model has
  a (log) likelihood that is at least as large as that of the \code{nullmodel}.
}

\note{
  Note that the likelihood ratio test may take a huge amount of time.
}

\section{Methods}{
 \describe{
   \item{print}{prints the summary}
   \item{summary}{gives a summary}
 }
}


\value{
  The test returns a message whether the null hypothesis, i.e. the
  smaller model is accepted.
  Invisibly, a list that also contains
  \itemize{
    \item \code{p}, the \eqn{p}-value
    \item \code{n}
    \item \code{data.ratio} the log ratio for the data
    \item \code{simu.ratio} the log ratio for the simulations
    \item \code{data.fit} the models fitted to the data
    \item \code{msg} the message that is also directly returned
  }
  It has S3 class \code{"RFratiotest"}.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\note{This function does not depend on the value of
 \command{\link{RFoptions}}\code{()$PracticalRange}. 
 The function \command{RFratiotest} always uses the standard specification
 of the covariance model as given in \command{\link{RMmodel}}.
}
\seealso{
  \command{\link{RFfit}}
 \command{\link{RMmodel}},
 \code{\link[=RandomFields-package]{RandomFields}},
 \command{\link{weather}}.
}
\examples{

\dontshow{\dontrun{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

data(soil)  ## see also ?soil
soil <- RFspatialPointsDataFrame(
 coords = soil[, c("x.coord", "y.coord")],
 data = soil[, c("moisture", "NO3.N", "Total.N",
 "NH4.N", "DOC", "N20N")],
 RFparams=list(vdim=6, n=1)
 )

model <- ~1 + RMplus(RMwhittle(scale=NA, var=NA, nu=NA), RMnugget(var=NA))
submodel <- ~1 + RMplus(RMwhittle(scale=NA, var=NA, nu=NA), RMnugget(var=0))

RFratiotest(submodel, model, data=soil["moisture"],
            modus_operandi="sloppy")
}}

\dontshow{FinalizeExample()}

}
\keyword{spatial}
\keyword{htest}

