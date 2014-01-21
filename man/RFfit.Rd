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
RFfit(model, x, y = NULL, z = NULL, T = NULL, grid, data,
 lower = NULL, upper = NULL, BC.lambda, methods,
 sub.methods, optim.control = NULL, users.guess = NULL,
 distances = NULL, dim, transform = NULL, ...)
}

\arguments{
 \item{model}{
 covariance model, see \command{\link{RMmodel}} or
 type \code{\link{RFgetModelNames}()} to get all options.
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
 \item{grid}{boolean. Whether coordinates give a grid}
 \item{BC.lambda}{a vector of at most two numerical components (just
 one component 
 corresponds to two identical ones) which are the parameters of the 
 box-cox-transformation:
 \eqn{\frac{x^\lambda_1-1}\lambda_1+\lambda_2}{(x^\lambda_1-1)/\lambda_1+\lambda_2}
 If the model is univariate, the first parameter can be estimated by
 using \code{NA}.
 }
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
   e.g. a parameter is supposed to be twice another one.
   This option should not be used yet.
 }
 \item{...}{further options and control parameters for the simulation
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
 \url{http://ms.math.uni-mannheim.de}}

\note{This function does not depend on the value of
 \command{\link{RFoptions}}\code{()$PracticalRange}. 
 The function \code{RFfit} always uses the standard specification
 of the covariance model as given in \command{\link{RMmodel}}.
}
\seealso{
  \command{\link{RFratiotest}}
 \command{\link{RMmodel}},
 \code{\link[=RandomFields-package]{RandomFields}},
 \command{\link{weather}}.
}
\examples{
RFoptions(seed=0)
% options(error=recover)
% source("randomfields_2/tests/source.R")

\dontshow{{ %if (!.C("isAuthor", a=integer(1))$a) {}
RFoptions(modus_operandi="sloppy")
}}

n <- if (interactive()) 100 else 5

#########################################################
## simulate some data first
points <- if (interactive()) 100 else 40
x <- runif(points, 0, 3)
y <- runif(points, 0, 3) ## random points in square [0, 3]^2
model <- RMgencauchy(alpha=1, beta=2)
d <- RFsimulate(model, x=x, y=y, grid=FALSE, n=n) #1000


#########################################################
## estimation; 'NA' means: "to be estimated"
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=2) +
 RMtrend(mean=NA)
v <- RFfit(estmodel, data=d, spC=FALSE)
Print(v$ml, v$ml$ml.value)


\dontshow{\dontrun{

#########################################################
## Estimation with fixed sill (variance + nugget equals a
## given constant)
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA) +
 RMnugget(var=NA) + RMtrend(mean=NA)
v1 <- RFfit(estmodel, data=d, fit.sill=1, fit.optim_var_estimation="try",
            spC=FALSE)
Print(v1$ml, v1$ml$ml.value)



#########################################################
## estimation in a anisotropic framework
x <- y <- (1:3)/4
model <- RMexp(Aniso=matrix(nc=2, c(4,2,-2,1)), var=1.5)
z <- RFsimulate(model, x=x, y=y, grid=TRUE, n=n, spC=TRUE)
estmodel <- RMexp(Aniso=matrix(nc=2, c(NA,NA,-2,1)), var=NA) +
 RMtrend(mean=NA)
v2 <- RFfit(estmodel, data=z, fit.nphi=20, spC=FALSE)
Print(v2$ml, v2$ml$ml.value)
}}


\dontshow{ \dontrun{

#########################################################
## estimation of coupled parameters (alpha = beta, here)
# source("RandomFields/tests/source.R")
f <- function(param) param[c(1,2,3,3,4)]
RFfit(estmodel, data=d, transform=list(), spC=FALSE)# shows positions of the NAs
v3 <- RFfit(estmodel, data=d, spC=FALSE,
 transform=list(c(rep(TRUE, 3), FALSE, TRUE), f))
Print(v3$ml, v3$ml$ml.value)

########################################################
## estimation of anisotropy matrix ##
x <- seq(0, 5, 1.5)
model <- RMexp(var=1, scale=1)
z <- RFsimulate(model, x, x, x, grid=TRUE, n=n, spC=TRUE)
estmodel <- RMexp(var=NA, Aniso=diag(rep(NA, 3))) +
 RMnugget(var=NA) + RMtrend(mean=NA)
v4 <- RFfit(estmodel, data=z, spC=FALSE)
Print(v4$ml, v4$ml$ml.value)

########################################################
## estimation of anisotropy matrix where the two first##
## diagonal elements are identical ##
RFfit(estmodel, data=z, transform=list(), spC=FALSE)
 ## shows the positions of the NAs
trafo <- function(variab) { variab[c(1, 2, 2, 3, 4)] }
v5 <- RFfit(estmodel, data=z, critical=0, spC=FALSE,
 transform=list(c(TRUE, TRUE, FALSE, TRUE, TRUE),
 trafo))
Print(v5$ml, v5$ml$ml.value)

########################################################
## estimation of anisotropy matrix where all ##
## diagonal elements are identical ##
trafo <- function(variab) { variab[c(1, 2, 2, 2, 3)] }
v6 <- RFfit(estmodel, data=z, critical=0, spC=FALSE,
 transform=list(c(TRUE, TRUE, FALSE, FALSE, TRUE),
 trafo),
 split=FALSE)
Print(v6$ml, v6$ml$ml.value)
}}



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
