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
 \item{boxcox}{the one or two parameters of the box cox transformation,
   in the univariate case.
   If not given, the globally defined parameters are used, see Details.
   In the \eqn{m}-variate case \code{boxcox} should be a \eqn{2 \times
   m} matrix.
 }
 \item{vdim}{the multivariate dimensionality of the field;
 }
 \item{inverse}{logical. Wether the inverse transformation should be performed.
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
  \command{\link{RFempiricalvariogram}}.
  Always first, the Box-Cox transformation is applied to the data.
  Then the command is performed. The result is back-transformed before
  returned. 
  

  If the first value of the transformation is \code{Inf} no
  transformation is performed (and is identical to \code{boxcox = c(1,
    0)}).
  If \code{boxcox} has length 1, then the transformation parameters
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

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{	
\link{Bayesian},
 \command{\link{RMmodel}},
  \command{\link{RFsimulate}},
 \command{\link{RFlikelihood}}.
}

\references{
  For the likeilhood correction see 
  \itemize{
    \item Konishi, S.,  and Kitagawa, G. (2008)
    \emph{Information criteria and statistical modeling.}
    Springer Science & Business Media. Section 4.9.
  }
  
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}

data(soil)
str(soil)
soil <- RFspatialPointsDataFrame(
 coords = soil[ , c("x.coord", "y.coord")],
 data = soil[ , c("moisture", "NO3.N", "Total.N", "NH4.N", "DOC", "N20N")],
 RFparams=list(vdim=6, n=1)
)
data <- soil["moisture"]
\dontshow{if (RFoptions()$internal$examples_red) {
  warning("data have been reduced !")
  All <- 1:7 
  data(soil)
  soil <- RFspatialPointsDataFrame(
     coords = soil[All, c("x.coord", "y.coord")],
     data = soil[All, c("moisture", "NO3.N", "Total.N",
      "NH4.N", "DOC", "N20N")],
      RFparams=list(vdim=6, n=1)
  )
  data <- soil["moisture"]
}}

model <- ~1 + RMplus(RMwhittle(scale=NA, var=NA, nu=NA), RMnugget(var=NA))

## Assuming log-Gaussian Data
print(fit <- RFfit(model, data=data, loggaus=TRUE))

## main Parameter in the Box Cox transformation to be estimated
print(fit <- RFfit(model, data=data, boxcox=NA))


\dontshow{FinalizeExample()}

}
\keyword{spatial}






