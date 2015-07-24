\name{RFfit}
\alias{RFfit}
%\alias{mleRF}% obsolete
\alias{RFfit.default}
\title{Fitting model parameters to spatial data (regionalised variables)
  and to linear (mixed) models}
\description{
  The function estimates arbitrary parameters of
  a random field specification with various methods.
  Currenty, the model to be fitted can be
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
 distances = NULL, dim, transform = NULL, ...)
}

\arguments{
 \item{model}{
 covariance model, see \command{\link{RMmodel}} or
 type \code{\link{RFgetModelNames}(type="variogram")} to get all options.
% See also t
 
 All parameters that are set to \code{NA} will be estimated;
 see the examples below.
 }
 \item{x}{vector of x coordinates, or object
   of class \code{\link[sp:GridTopology-class]{GridTopology}} or
   \code{\link[raster]{raster}};
 For more options see \link{RFsimulateAdvanced}.
 }
 \item{y}{vector of y coordinates}
 \item{z}{vector of z coordinates}
 \item{T}{vector of T coordinates; these coordinates are given in
 triple notation, see \command{\link{RFsimulate}}.}
 \item{data}{vector or matrix of values measured at \code{coord};
 If a matrix is given then the columns are interpreted as independent
 realisations.\cr
 If also a time component is given, then in the data the indices for
 the spatial components run the fastest.
 
 If an \code{m}-variate model is used, then each realisation is given as
 \code{m} consecutive columns of \code{data}.
 }
 \item{lower}{list or vector. Lower bounds for the parameters.
 If \code{param} is a vector, \code{lower} has to be a vector as well and
 its length must equal the number of parameters to be estimated. The order
 of \code{param} has to be maintained. A component being \code{NA} means
 that no manual lower bound for the corresponding parameter is set.
 
 If \code{param} is a list, \code{lower} has to be of (exactly) the same
 structure.
 }
 \item{upper}{list or vector. Upper bounds for the parameters.
 See also lower.
 }
 \item{grid}{boolean. Whether coordinates build a grid. In nearly all cases
   the value of \code{grid}
   is self-detected, so that \code{grid} need not be given,
   in general}
 

 \item{methods}{
  Main methods to be used for estimating. If several methods,
   estimation will be performed with each method and the results reported.
 }
 \item{sub.methods}{
 variants of the least squares fit of the variogram. See Details.
 variants of the maximum likelihood fit of the covariance function.
 See Details.
 }
 \item{users.guess}{
 User's guess of the parameters. All the parameters must be given
 using the same rules as for either \code{param} (except that no NA's should
 be contained) or \code{model}.
 }
 \item{distances,dim}{Instead of \code{x}-coordinates,  distances might
 be given. The the dimension of the space \code{dim} must be given explicitely.}
 \item{optim.control}{control list for \command{\link[stats]{optim}},
 which uses \sQuote{L-BFGS-B}. However \code{parscale} may not be given.}
 \item{transform}{
   this is an attempt to allow binding between
   parameters,
   e.g. one parameter value is supposed to equal another one,
   See examples below.

   \code{transform=list()} is not valid for estimating, but returns
   structural information to set up the correct function.
  See examples below.

 }
 \item{...}{further options and control arguments for the simulation
 that are passed to and processed by \command{\link{RFoptions}}.}
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
  \code{"ml"} and \code{c("self", "plain", "sqrt.nr", "sd.inv",
    "internal")},
  respectively.
}


\note{
  \itemize{
    \item
    An important optional argument is \code{boxcox} which indicates
    A Box-Cox transformation see \code{boxcox} in \command{\link{RFoptions}}
    and \code{RFboxcox} for details.
    \item
    Instead of \command{\link[stats]{optim}}, other optimiser can be used,
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
    of the covariance model as given in \command{\link{RMmodel}}
  }
  }

\value{
  The result depends on the logical value of
  \code{\link[=RFoptions]{spConform}}.
  If \code{TRUE}, an  S4 object is greated. In case the model indicates
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

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}


\seealso{
 \link{RFfitOptimiser}
  \command{\link{RFlikelihood}},
 \command{\link{RFratiotest}},
 \command{\link{RMmodel}},
 \code{\link[=RandomFields-package]{RandomFields}},
 \command{\link{weather}}.
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
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
\dontshow{FinalizeExample()}

}
\keyword{spatial}
\keyword{optimize}

