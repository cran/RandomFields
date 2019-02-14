\name{RFboxcox}
\alias{RFboxcox}
\title{Linear part of \command{\link{RMmodel}}}
\description{
  \command{\link{RFboxcox}} performs the Box-Cox transformation:
  \eqn{\frac{(x+\mu)^\lambda-1}{\lambda}}
}
\usage{
RFboxcox(data, boxcox, vdim = 1, inverse=FALSE, ignore.na=FALSE)
}
\arguments{
 \item{data}{matrix or list of matrices.
 }
 \item{boxcox}{the one or two parameters \eqn{(\lambda, \mu)}
   of the box cox transformation,
   in the univariate case; if \eqn{\mu} is not given, then \eqn{\mu} is
   set to \eqn{0}.
   If not given, the globally defined parameters are used, see Details.
   In the \eqn{m}-variate case \code{boxcox} should be a \eqn{2 \times
   m} matrix. If \eqn{\lambda =\infty} then no transformation is performed.
 }
 \item{vdim}{the multivariate dimensionality of the field;
 }
 \item{inverse}{logical. Whether the inverse transformation should be performed.
 }
 \item{ignore.na}{logical. If \code{FALSE} an error message is returned
   if any value of \code{boxcox} is \code{NA}. Otherwise the data are
   returned without being transformed.
 }
}

\details{
  The Box-Cox transfomation \code{boxcox} can be set
  globally through \command{\link{RFoptions}}. If it is set globally the
  transformation applies in the \bold{Gaussian} case to
  \command{\link{RFfit}},
  \command{\link{RFsimulate}},
  \command{\link{RFinterpolate}},
  \command{\link{RFvariogram}}.
  Always first, the Box-Cox transformation is applied to the data.
  Then the command is performed. The result is back-transformed before
  returned. 
  

  If the first value of the transformation is \code{Inf} no
  transformation is performed (and is identical to \code{boxcox = c(1,0)}).
  If \code{boxcox} has length 1, then the transformation parameter
  \eqn{\mu} is set to  \eqn{0}, which is the standard case.
}

\value{  
  \command{\link{RFboxcox}} returns a list
  of three components, \code{Y}, \code{X}, \code{vdim} returning
  the deterministic trend, the design matrix, and the multivariability,
  respectively.
  If \code{set} is positive, \code{Y} and \code{X} contain
  the values for the \code{set}-th set of coordinates.
  Else, \code{Y} and \code{X} are both lists containing
  the values for all the sets.
}

\me
\seealso{	
\link{Bayesian},
 \command{\link{RMmodel}},
  \command{\link{RFsimulate}},
 \command{\link{RFlikelihood}}.
}

\references{
  For the likelihood correction see 
  \itemize{
    \item Konishi, S.,  and Kitagawa, G. (2008)
    \emph{Information criteria and statistical modeling.}
    Springer Science & Business Media. Section 4.9.
  }
  
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

data(soil)
str(soil)
soil <- RFspatialPointsDataFrame(
 coords = soil[ , c("x.coord", "y.coord")],
 data = soil[ , c("moisture", "NO3.N", "Total.N", "NH4.N", "DOC", "N20N")],
 RFparams=list(vdim=6, n=1)
)
dta <- soil["moisture"]
\dontshow{if (RFoptions()$internal$examples_red) {
  warning("data have been reduced !")
  All <- 1:7
  rm(soil)
  data(soil)
  soil <- RFspatialPointsDataFrame(
     coords = soil[All, c("x.coord", "y.coord")],
     data = soil[All, c("moisture", "NO3.N", "Total.N",
      "NH4.N", "DOC", "N20N")],
      RFparams=list(vdim=6, n=1)
  )
  dta <- soil["moisture"]
}}

model <- ~1 + RMplus(RMwhittle(scale=NA, var=NA, nu=NA), RMnugget(var=NA))

\dontshow{\dontrun{
## Assuming log-Gaussian Data
print(fit <- RFfit(model, data=dta, loggaus=TRUE))
}}

## main Parameter in the Box Cox transformation to be estimated
print(fit <- RFfit(model, data=dta, boxcox=NA))

\dontshow{FinalizeExample()}}

\keyword{spatial}






