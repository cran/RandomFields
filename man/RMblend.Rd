\name{RMblend}
\alias{RMblend}
\title{Scale model for a few areas of different scales and/or
  differentiabilities} 
\description{
  Let \eqn{Z=(Z_1, \ldots Z_k)} be an \eqn{k}-variate random field
  and \eqn{A_1,\ldots, A_k} a partition of the space.
  Then
  \deqn{Y(x) = \sum_{i=1}^k Z_i * 1(x \in A_i)}
  i.e. the model blends the components of \eqn{Z} to a new, univariate
  model \eqn{Y}.
}
\usage{
RMblend(multi, blend, thresholds, var, scale, Aniso, proj)
}
\arguments{
 \item{multi}{a multivariate covariance function}
 \item{blend,thresholds}{The \code{threshold} is a vector of increasing
   values. If the value of \code{blend} is below all thresholds up to
   the \eqn{k}-th threshold, then the \eqn{k}-th component of the field
   given by \code{multi} is taken. If necessary the components are recycled.

   Default: \code{threshold = 0.5}, useful for blending a bivariate
   field if \code{blend} takes only the values \eqn{0} and {1}.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 \command{\link{RMblend}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}

\references{
  \itemize{
    \item Bonat, W.H. , Ribeiro, P. Jr. and Schlather, M. (2019)
    Modelling non-stationarity in scale. In preparation.
    \item Genton, Apanovich Biometrika.
  }
}


\me

\seealso{
  \command{\link{RMSadvanced}},
  \command{\link{RMbubble}},
  \command{\link{RMscale}},
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- seq(0,1, if (interactive()) 0.01 else 0.5)
len <- length(x)
m <- matrix(1:len, nc=len, nr=len)
m <- m > t(m)
image(m) # two areas separated by the first bisector

biwm <- RMbiwm(nudiag=c(0.3, 1), nured=1, rhored=1, cdiag=c(1, 1), 
                s=c(1, 1, 0.5))
model <- RMblend(multi=biwm, blend=RMcovariate(data = as.double(m), raw=TRUE))
plot(z <- RFsimulate(model, x, x)) ## takes a while ...
\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}

