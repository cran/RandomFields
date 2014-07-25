\name{RMfractgauss}
\alias{RMfractgauss}
\title{Fractal Gaussian Model Family}
\description{
 \command{\link{RMfractgauss}} is a stationary isotropic covariance model.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = 0.5 ((r+1)^{\alpha}-2r^{\alpha}+|r-1|^{\alpha})}{C(r) = 0.5 ((r+1)^(\alpha)-2r^(\alpha)+|r-1|^(\alpha))}
 with \eqn{0 < \alpha \le 2}. It can only be used for one-dimensional random fields.
}
\usage{
RMfractgauss(alpha,var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{ \eqn{0 < \alpha \le 2}}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{The model is only valid for dimension \eqn{d = 1}{d = 1 }. It
 is the covariance function for the fractional Gaussian noise with
 self-affinity index (Hurst parameter) \eqn{H=\alpha /2}{H=\alpha /2} with \eqn{0 < \alpha \le 2}{0 < \alpha \le 2}.
}
\value{
 \command{\link{RMfractgauss}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}
\references{
 \itemize{
 \item Gneiting, T. and Schlather, M. (2004)
 Stochastic models which separate fractal dimension and Hurst effect.
 \emph{SIAM review} \bold{46}, 269--282.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMfractgauss(alpha=0.5, scale=0.2)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(-0.3, 1))
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}