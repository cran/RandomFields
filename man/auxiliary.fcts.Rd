\name{CheckXT}
\alias{CheckXT}
\alias{PrepareModel}
\alias{plotWithCircles}
\alias{circle}
\alias{GetDistributionNames}
\alias{GetrfParameters}
\title{Internal functions -- do not use them directly}
\description{
  \code{CheckXT} checks whether the coordinates of the data and related
  parameters are specified correctly and transforms the coordinates into
  a standard format

  \code{PrepareModel} checks whether the parameters of the covariance
  model and related parameters are specified correctly and transforms
  the parameters into a standard format

  \code{plotWithCircles} displays data values of marked point processes
  by circles

  \code{GetDistributionNames} returns the names of the currently
  available marginal distributions of the random fields
  
  \code{GetrfParameters} returns some internal constants
  
 }
\usage{
CheckXT(x, y, z, T, grid, gridtriple)
PrepareModel(model, param, trend, method=NULL, nugget.remove=TRUE)
plotWithCircles(data, factor=1.0, xlim=range(data[,1])+c(-maxr,maxr),
                ylim=range(data[,2])+c(-maxr,maxr),col=1, fill=0, ...)
GetDistributionNames()
GetrfParameters(initcov)
}
\arguments{
  \item{x}{\code{x} coordinates}
  \item{y}{\code{y} coordinates}
  \item{z}{\code{z} coordinates}
  \item{T}{time instances}
  \item{grid}{logical; simulation on a grid; see \command{\link{GaussRF}}}
  \item{gridtriple}{logical; grid is given in abbreviated form;
    see \command{\link{GaussRF}}}
  \cr 
  \item{model}{name of the covariance model; see \command{\link{GaussRF}}}
  \item{param}{parameter for the covariance model; see \command{\link{GaussRF}}}
%  \item{timespacedim}{dimension of the random field including the time
%    dimension, if there is any}
  \item{trend}{mean or trend of the random field}
  \item{method}{simulation method}
  \item{nugget.remove}{if nugget is 0, the nugget part be removed in the
    returned model?}
  \cr
  \item{data}{matrix of 3 columns; first two columns give the
    coordinates, the third the data}
  \item{factor}{enlargement factor for data}
  \item{xlim}{see \command{\link[graphics]{plot}}}
  \item{ylim}{see \command{\link[graphics]{plot}}}
  \item{col}{border colour of circles}
  \item{fill}{filling colour of circles}
  \item{...}{further graphical parameters in \command{plotWithCircles}
    and any variables to be printed in \command{Print}}
  \item{initcov}{boolean, whether to call the C function InitModelList}
}
%\details{
%}
%\value{
%  lists of internal parameters
%}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\seealso{\command{\link{CovarianceFct}}}
\keyword{ spatial }%-- one or more ...
\keyword{internal}
\examples{
% library(RandomFields)

model <- list(model="whittle", param=c(1,2,3,4,5))
Print(PrepareModel(model=model, me="ci"))

model <- list("+",
              list("$", var=2, s=4, list(model="whi", kappa=5)),
              list("$", var=3, s=0, list(model="whi", kappa=1)) )
Print(PrepareModel(model=model, me="ci"))

} 

