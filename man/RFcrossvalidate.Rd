\name{RFcrossvalidate}
\alias{RFcrossvalidate}
\alias{print.crossvalidate}
\alias{summary.crossvalidate}
\alias{print.summary.crossvalidate}
%\alias{mleRF}% obsolete
\alias{RFcrossvalidate.default}
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
\usage{D
RFcrossvalidate(model, x, y = NULL, z = NULL, T = NULL, grid, data,
 lower = NULL, upper = NULL, bc_lambda, method="ml",
 users.guess = NULL,
 distances = NULL, dim, optim.control = NULL,transform = NULL,
 full = FALSE, ...)
}

\arguments{
 \item{model, x, y, z, T, grid, data,lower, upper, bc_lambda, users.guess,
 distances, dim, optim.control, transform, ...}{see \command{\link{RFfit}}
 }
 \item{method}{
   Single method to be used for estimating, either
   one of the \code{methods} or one of the \code{sub.methods}
   see \command{\link{RFfit}}
 }
 \item{full}{logical.
   if \code{TRUE} then crossvalidatiaon is also performed
   for intermediate models used in
   \code{RFfit} (if any).
   }
}
 
\section{Methods}{
 \describe{
   \item{print}{prints the summary}
   \item{summary}{gives a summary}
 }
}

\note{
  An important option is \code{cross_refit} that determines
  whether the model is refitted for each location left out.
  Default is \code{FALSE}. See also \command{\link{RFoptions}}.
}


\value{
  An object of the \code{\link{class}} \code{"RFcrossvalidate"} which is
  a list with the following components, cf. \command{xvalid} in the
  package \pkg{geoR} :
  \item{data}{the original data.  }
  \item{predicted}{the values predicted by cross-validation.  }
  \item{krige.var}{the cross-validation prediction variance.  }
  \item{error}{the differences \code{data - predicted value}.   }
  \item{std.error}{the errors divided by the square root of the
    prediction variances.  }
  \item{p}{
    In contrast to \pkg{geoR} the p-value is returned,
    i.e. the probability
    that a difference with absolute value larger than the absolute
    value of the actual difference is observed.


    A method for \code{summary} returns summary statistics for the errors
    and standard errors similar to \pkg{geoR}.
  
    If \code{cross_refit = TRUE} and \code{detailed_output = TRUE}
    the returned object also constains a \code{fitted} which is
    a list of fitted models.
  }
}

\references{
  \itemize{
    \item Ribeiro, P.J., Jr. and Diggle, P.J (2014) R package \pkg{geoR}.
    
    \item Burnham, K. P. and Anderson, D. R. (2002)
    \emph{Model selection and Multi-Model Inference: A Practical
      Information-Theoretic Approach.}
    2nd edition. New York: Springer.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\note{This function does not depend on the value of
 \command{\link{RFoptions}}\code{()$PracticalRange}. 
 The function \code{RFcrossvalidate} always uses the standard specification
 of the covariance model as given in \command{\link{RMmodel}}.
}
\seealso{
  \command{\link{RFratiotest}}
  \command{\link{RFfit}}
  \command{\link{RMmodel}},
  \code{\link[=RandomFields-package]{RandomFields}},
  \command{\link{weather}}.
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
% options(error=recover)
% source("RandomFields/tests/source.R")

RFoptions(modus_operandi="sloppy")

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
RFcrossvalidate(estmodel, data=d, spC=FALSE)


\dontshow{\dontrun{

#########################################################
## Estimation with fixed sill (variance + nugget equals a
## given constant)
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA) +
 RMnugget(var=NA) + RMtrend(mean=NA)
RFcrossvalidate(estmodel, data=d, fit.sill=1, fit.optim_var_elimination="try",
            spC=FALSE)



#########################################################
## estimation in a anisotropic framework
x <- y <- (1:3)/4
model <- RMexp(Aniso=matrix(nc=2, c(4,2,-2,1)), var=1.5)
z <- RFsimulate(model, x=x, y=y, grid=TRUE, n=n, spC=TRUE)
estmodel <- RMexp(Aniso=matrix(nc=2, c(NA,NA,-2,1)), var=NA) +
 RMtrend(mean=NA)
RFcrossvalidate(estmodel, data=z, fit.nphi=20, spC=FALSE)
}}


\dontshow{ \dontrun{

#########################################################
## estimation of coupled parameters (alpha = beta, here)
# source("RandomFields/tests/source.R")
f <- function(param) param[c(1,2,3,3,4)]
RFcrossvalidate(estmodel, data=d, transform=list(), spC=FALSE)# shows positions of the NAs
RFcrossvalidate(estmodel, data=d, spC=FALSE,
 transform=list(c(rep(TRUE, 3), FALSE, TRUE), f))

########################################################
## estimation of anisotropy matrix ##
x <- seq(0, 5, 1.5)
model <- RMexp(var=1, scale=1)
z <- RFsimulate(model, x, x, x, grid=TRUE, n=n, spC=TRUE)
estmodel <- RMexp(var=NA, Aniso=diag(rep(NA, 3))) +
 RMnugget(var=NA) + RMtrend(mean=NA)
RFcrossvalidate(estmodel, data=z, spC=FALSE)


########################################################
## estimation of anisotropy matrix where the two first##
## diagonal elements are identical ##
RFcrossvalidate(estmodel, data=z, transform=list(), spC=FALSE)
 ## shows the positions of the NAs
trafo <- function(variab) { variab[c(1, 2, 2, 3, 4)] }
RFcrossvalidate(estmodel, data=z, critical=0, spC=FALSE,
 transform=list(c(TRUE, TRUE, FALSE, TRUE, TRUE),
 trafo))


########################################################
## estimation of anisotropy matrix where all ##
## diagonal elements are identical ##
trafo <- function(variab) { variab[c(1, 2, 2, 2, 3)] }
RFcrossvalidate(estmodel, data=z, critical=0, spC=FALSE,
 transform=list(c(TRUE, TRUE, FALSE, FALSE, TRUE),
 trafo),
 split=FALSE)
}}



\dontshow{RFoptions(modus_operandi="normal")}
\dontshow{FinalizeExample()}

}
\keyword{spatial}

