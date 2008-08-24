\name{CheckXT}
\alias{CheckXT}
\alias{PrepareModel}
\alias{convert.to.readable}
\alias{plotWithCircles}
\alias{circle}
\alias{GetDistributionNames}
\alias{paramextract}
\title{Internal functions -- do not use them directly}
\description{
  \code{CheckXT} checks whether the coordinates of the data and related
  parameters are specified correctly and transforms the coordinates into
  a standard format

  \code{PrepareModel} checks whether the parameters of the covariance
  model and related parameters are specified correctly and transforms
  the parameters into a standard format

  \code{convert.to.readable} is the inverse function to
  \code{PrepareModel}; see Details

  \code{plotWithCircles} displays data values of marked point processes
  by circles

  \code{GetDistributionNames} returns the names of the currently
  available marginal distributions of the random fields

  \code{paramextract} extracts for some models some parameters
  from an internal parameter list
  

}
\usage{
CheckXT(x, y, z, T, grid, gridtriple)
PrepareModel(model, param, timespacedim, trend, method=NULL,
             named=FALSE)
convert.to.readable(l, allowed=c("standard", "nested", "list")) 
plotWithCircles(data, factor=1.0, xlim=range(data[,1])+c(-maxr,maxr),
                ylim=range(data[,2])+c(-maxr,maxr),col=1, fill=0, ...)
GetDistributionNames()
paramextract(p, model=c("cutoff"))
}
\arguments{
  \item{x}{\code{x} coordinates}
  \item{y}{\code{y} coordinates}
  \item{z}{\code{z} coordinates}
  \item{T}{time instances}
  \item{grid}{see \command{\link{GaussRF}}}
  \item{gridtriple}{see \command{\link{GaussRF}}}
  \cr 
  \item{model}{see \command{\link{GaussRF}}}
  \item{param}{see \command{\link{GaussRF}}}
  \item{timespacedim}{dimension of the random field including the time
    dimension, if there is any}
  \item{trend}{mean or trend of the random field}
  \item{method}{simulation method}
  \item{named}{logical. If \code{TRUE} \code{covnr} and \code{param}
    are returned with names}
 \cr
  \item{l}{list as returned by \code{PrepareModel}}
  \item{allowed}{allowed output formats, see
    \command{\link{CovarianceFct}}}
  \cr
  \item{data}{matrix of 3 columns; first two columns give the
    coordinates, the third the data}
  \item{factor}{enlargement factor for data}
  \item{xlim}{see \command{\link[graphics]{plot}}}
  \item{ylim}{see \command{\link[graphics]{plot}}}
  \item{col}{border colour of circles}
  \item{fill}{filling colour of circles}
  \item{...}{further graphical parameters}
  \item{p}{internal parameter list; e.g. the columns of
    \command{\link{CheckAndComplete}}\code{(...)$param}.
  }
  \item{model}{the name of a covariance model.}
}
\details{
  \code{convert.to.readable} is roughly speaking the inverse function to
  \code{PrepareModel}.  \code{convert.to.readable} also tries to 
  simplify the model definition, but cannot rediscover the given method for
  the simulation of the nugget effect in all cases.  Due to the
  simplification in \code{convert.to.readable} and the special
  definition of the nugget effect for nested models,
  \code{convert.to.readable} may return a correct model definition in case
  of incorrect input, namely if \code{scale} is set to \eqn{0} in a list
  definition, see Examples.
}
%\value{
%  lists of internal parameters
%}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
\seealso{\command{\link{CovarianceFct}}}
\keyword{ spatial }%-- one or more ...
\keyword{internal}
\examples{
% library(RandomFields)

x <- function(...) {
  str(PrepareModel(...))
  cat("--------------------------------\n")
  str(convert.to.readable(PrepareModel(...)))
}
model <- list(list(model="whi", kappa=5, var=2, s=4), "+",
    list(model="whi", kappa=1, var=3, s=0)) ## s=0 should not be used only in
##                             a model definition where the parameters are
##                             are given in a matrix, see the result
x(model=model, ti=1, me="ci")

## since convert.to.readable performs a one-step simplification,
## iterative calls may further simplify the model
xx <- convert.to.readable(PrepareModel(model=model, ti=1, me="ci"))
x(model=xx$mo, pa=xx$pa, ti=1, me=xx$me)

## back to the matrix definition of nested models
str(convert.to.readable(PrepareModel(xx, ti=1), allowed="nested"))

## back to the (correct) list definition
str(convert.to.readable(PrepareModel(xx, ti=1), allowed="list"))
}


