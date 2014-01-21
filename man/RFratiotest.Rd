\name{RFratiotest}
\alias{RFratiotest}
%\alias{mleRF}% obsolete
\title{Likelihood ratio test}
\description{
  The function performs a Monte Carlo likelihood ratio test
  based on \command{\link{fitgauss}}. Currently it only works for
    Gaussian random fields.
}    
\usage{
RFratiotest(nullmodel, alternative, x, y = NULL, z = NULL, T = NULL, grid, data,
alpha=0.05, n = 5 / alpha, seed = 0, 
lower = NULL, upper = NULL, BC.lambda, methods,
sub.methods, optim.control = NULL, users.guess = NULL,
distances = NULL, dim, transform = NULL, ...)
}
\arguments{
  \item{nullmodel, alternative}{
    covariance model, see \command{\link{RMmodel}} or
 type \code{\link{RFgetModelNames}()} to get all options.

 
    The set of parameters to estimate for \code{nullmodel}
    should be a subset of the parameters  to estimate for \code{alternative}.
    }
  \item{alpha}{value in [0,1]. Significance level.
 }
 \item{n}{integer. The test is based on \code{n-1} simulations.}
 \item{seed}{integer. If not \code{NULL} and not \code{NA}, the
   \link[base]{.Random.seed} is set to \code{seed}. Otherwise,
   \command{\link[base]{set.seed}} is set to the value of 
   \code{RFoptions{}$general$seed}
   if the latter is not \code{NA}.
 }
 \item{x, y, z, T, grid, data, lower, upper, BC.lambda, methods,
   sub.methods, optim.control, users.guess, distances, dim, transform,
   ...}{see \command{\link{RFfit}}}
}
 
\details{
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
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}}

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
% options(error=recover)
% source("randomfields_2/tests/source.R")

\dontshow{{ %if (!.C("isAuthor", a=integer(1))$a) {}
RFoptions(modus_operandi="sloppy")
}}

\dontrun{
% library(RandomFields, lib="~/TMP");  source("randomfields_2/private/source.R")

RFoptions(seed=0)

## see also ?soil
data(soil)
str(soil)
soil <- RFspatialPointsDataFrame(
 coords = soil[, c("x.coord", "y.coord")],
 data = soil[, c("moisture", "NO3.N", "Total.N",
 "NH4.N", "DOC", "N20N")],
 RFparams=list(vdim=6, n=1)
 )

model <- ~1 + RMplus(RMwhittle(scale=NA, var=NA, nu=NA), RMnugget(var=NA))
submodel <- ~1 + RMplus(RMwhittle(scale=NA, var=NA, nu=NA), RMnugget(var=0))

test <- RFratiotest(submodel, model, data=soil["moisture"], n=20)
str(test)
}

\dontshow{{ %if (!.C("isAuthor", a=integer(1))$a) {}
RFoptions(modus_operandi="normal")
}}

\dontshow{RFoptions(seed=NA)}

}
\keyword{spatial}







% LocalWords: fitvario mleRF LSQ param NA naturalscaling PrintLevel optim pch
% LocalWords: RFoptions nphi ntheta ntime lowerbound minbounddistance var
% LocalWords: minboundreldist functioncalls eqn coord RMmodel emph cr
% LocalWords: RFgetModelNames cov variogram uote parampositions NaN kappas MLE
