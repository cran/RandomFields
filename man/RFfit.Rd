\name{RFfit}
\alias{RFfit}
%\alias{mleRF}% obsolete
\alias{RFfit.default}
\title{Fitting model parameters to spatial data (regionalised variables)
  and to linear (mixed) models}
\description{
  The function estimates arbitrary parameters of
  a random field specification with various methods.
  Currently, the models to be fitted can be
  \itemize{
    \item{\link[=RPgauss]{Gaussian random fields}}
    \item{\link[=RFformula]{linear models}}
  }
  The fitting of max-stable random fields and others
  has not been implemented yet.
}
\usage{
RFfit(model, x, y = NULL, z = NULL, T = NULL, grid=NULL, data,
      lower = NULL, upper = NULL, methods,
      sub.methods, optim.control = NULL, users.guess = NULL,
      distances = NULL, dim, transform = NULL, params=NULL, ...)
}

\arguments{
  \item{model,params}{\argModel
  
 All parameters that are set to \code{NA} will be estimated;
 see the examples below.

 Type \code{\link{RFgetModelNames}(type="variogram")} to get all options
 for \code{model}. 
 }
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{data}{\argData}
 \item{lower}{\argLower}
 \item{upper}{\argUpper}
 \item{methods}{\argFitmethods} 
 \item{sub.methods}{\argFitsubmethods. See Details.}
 \item{users.guess}{\argUsersguess}
 \item{distances,dim}{\argDistances}
 \item{optim.control}{\argOptimcontrol}
 \item{transform}{\argTransform}
 \item{...}{\argDots}
}
 
\details{
  For details on the simulation methods see
  \itemize{
    \item \link{fitgauss} for \link[=RPgauss]{Gaussian random fields}
    \item \link{fitgauss} for \link[=RFformula]{linear models}
  }

  If \code{x}-coordinates are not given, the function will check
  \code{data} for \code{NA}s and will perform imputing.

  The function has many more options to tune the optimizer,
  see \command{\link{RFoptions}} for details.

  If the model defines a Gaussian random field, the options
  for \code{methods} and \code{submethods} are currently
  \code{"ml"} and \code{c("self", "plain", "sqrt.nr", "sd.inv", "internal")},
  respectively.
}


\note{
  \itemize{
    \item
    An important optional argument is \code{boxcox} which indicates
    a Box-Cox transformation; see \code{boxcox} in \command{\link{RFoptions}}
    and \code{RFboxcox} for details.
    \item
    Instead of \command{\link[stats]{optim}}, other optimisers can be used,
    see \link{RFfitOptimiser}.
    \item
    Several advanced options can be found in sections \sQuote{General
      options} and \sQuote{fit} of \command{\link{RFoptions}}.
    \item
    In particular, \code{boxcox}, \code{boxcox_lb},  \code{boxcox_ub}
    allow Box-Cox transformation.
    \item
    This function does not depend on the value of
    \command{\link{RFoptions}}\code{()$PracticalRange}. 
    The function \code{RFfit} always uses the standard specification
    of the covariance model as given in \command{\link{RMmodel}}.
  }
  }

\value{
  The result depends on the logical value of
  \code{\link[=RFoptions]{spConform}}.
  If \code{TRUE}, an  S4 object is created. In case the model indicates
    a Gaussian random field, an
  \link[=RFfit-class]{RFfit} object is created.

  If \code{spConform=FALSE}, a list is returned.
  In case the model indicates
  a Gaussian random field, the details are given in \link{fitgauss}.
  

}
\references{
 \itemize{
 \item Burnham, K. P. and Anderson, D. R. (2002)
 \emph{Model selection and Multi-Model Inference: A Practical
   Information-Theoretic Approach.}
 2nd edition. New York: Springer.
 }
}

\me

\seealso{
 \command{\link{RFfitOptimiser}},
  \command{\link{RFlikelihood}},
 \command{\link{RFratiotest}},
 \command{\link{RMmodel}},
 \code{\link[=RandomFields-package]{RandomFields}},
 \command{\link{weather}}.
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
% options(error=recover)
% source("RandomFields/tests/source.R")
% to do bsp mit box cox

RFoptions(modus_operandi="sloppy")


#########################################################
## simulate some data first                            ## 
points <- 100
x <- runif(points, 0, 3)
y <- runif(points, 0, 3) ## random points in square [0, 3]^2
model <- RMgencauchy(alpha=1, beta=2)
d <- RFsimulate(model, x=x, y=y, grid=FALSE, n=100) #1000


#########################################################
## estimation; 'NA' means: "to be estimated"           ##
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=2) +
            RMtrend(mean=NA)
RFfit(estmodel, data=d)


#########################################################
## coupling alpha and beta                             ##
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA) + 
            RMtrend(NA)
RFfit(estmodel, data=d, transform = NA) ## just for information
trafo <- function(a) c(a[1], rep(a[2], 2))
fit <- RFfit(estmodel, data=d,
             transform = list(c(TRUE, TRUE, FALSE), trafo))
print(fit)
print(fit, full=TRUE)


\dontshow{\dontrun{
#########################################################
## Estimation with fixed sill (variance + nugget       ##
## equals a given constant)                            ##
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA) +
            RMnugget(var=NA) + RMtrend(mean=NA)
RFfit(estmodel, data=d, fit.sill=1, fit.optim_var_elimination="try")


#########################################################
## estimation in a anisotropic framework               ##
x <- y <- (1:3)/4
model <- RMexp(Aniso=matrix(nc=2, c(4,2,-2,1)), var=1.5)
d <- RFsimulate(model, x=x, y=y, n=n)
estmodel <- RMexp(Aniso=matrix(nc=2, c(NA,NA,-2,1)), var=NA) +
            RMtrend(mean=NA)
RFfit(estmodel, data=d, fit.nphi=20)


#########################################################
## AN EXAMPLE HOW TO USE OF PARAMETER 'transform'      ##
## estimation of coupled parameters (first column of   ##
##           the matrix 'Aniso' has identical entries) ##
# source("RandomFields/tests/source.R")
RFfit(estmodel, data=d, transform=list()) # shows positions of NAs
f <- function(param) param[c(1,2,2)]
RFfit(estmodel, data=d, transform=list(c(TRUE, TRUE, FALSE), f))


}}

\dontshow{RFoptions(modus_operandi="normal")}
\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{optimize}

