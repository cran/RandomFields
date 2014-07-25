\name{RMdelay}
\alias{RMdelay}
\title{Bivariate Delay Effect}
\description{
 \command{\link{RMdelay}} is a \eqn{(m+1)}-variate stationary covariance model.
 which depends on a univariate stationary covariance model \eqn{C_0}.
 
 The corresponding covariance function only depends on the difference
 \eqn{h \in {\bf R}^d}{h} between two points in \eqn{d}-dimensional
 space and is given by
 \deqn{C(h)=(C_0(h - s_i +s_j))_{i,j=0,\ldots,m}}
 where \eqn{s \in {\bf R}^{d\times m}}{h \in R^{d x m}} and \eqn{s_0=0}
}

\usage{
RMdelay(phi,s,var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a univariate stationary covariance model, that means an
 \command{\link{RMmodel}} whose \command{vdim} equals 1.}
 \item{s}{a \eqn{d\times m}{d x m}-dimensional
 shift matrix, where \eqn{d} is the
 dimension of the space, giving the components \eqn{s=(s_1,\ldots,
 s_m)}{s=(s_1,..., s_m)} where the \eqn{s_i} are vectors.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{%See Wackernagel, H. (2003), p.?)
 Here, a multivariate random field is obtained from single univariate
 random field, by shifting it by fixed value.
}
\value{
 \command{\link{RMdelay}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Schlather, M., Malinowski, A., Menck, P.J., Oesting, M., Strokorb,
 K. (2013). \emph{Submitted to J. Statist. Software}
 
 \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 3nd edition.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

x <- y <- seq(-10,10,0.2)
model <- RMdelay(RMstable(alpha=1.9, scale=2), s=c(4,4))
simu <- RFsimulate(model, x, y, grid=TRUE)
plot(simu, zlim="joint")

\dontshow{FinalizeExample()}
}


\keyword{spatial}
\keyword{models}


