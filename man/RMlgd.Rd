\name{RMlgd}
\alias{RMlgd}
\title{Local-Global Distinguisher Family Covariance Model}
\description{
 \command{\link{RMlgd}} is a stationary isotropic covariance model, which is valid only for dimensions 
 \eqn{d =1,2}{d =1,2}.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) =1 - \beta^{-1}(\alpha + \beta)r^{\alpha} 1_{[0,1]}(r) + \alpha^{-1}(\alpha + \beta)r^{-\beta} 1_{r>1}(r) }{C(r) =1 - \beta^(-1)(\alpha + \beta)r^(\alpha) 1_{[0,1]}(r) + \alpha^(-1)(\alpha + \beta)r^(-\beta) 1_{r>1}(r) }
 where \eqn{\beta >0} and \eqn{0 < \alpha \le (3-d)/2}{0 < \alpha \le (3-d)/2},
with \eqn{d}{d} denoting the dimension of the random field.
}
\usage{
RMlgd(alpha, beta, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{argument whose range depends on the dimension of the random field: \eqn{0< \alpha \le (3-d)/2}{0< \alpha \le (3-d)/2}.}
 \item{beta}{beta > 0.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
The model is only valid for dimension \eqn{d=1,2}{d=1,2}.

This model admits simulating random fields where fractal dimension
\emph{D} of the Gaussian sample and Hurst coefficient \emph{H}
can be chosen independently (compare also \command{\link{RMgencauchy}}.):

Here, the random field has fractal dimension \deqn{D = d+1 - \alpha/2}{D = d+1
 - \alpha/2} and Hurst coefficient \deqn{H = 1-\beta/2}{H = 1-\beta/2} for \eqn{0< \beta \le 1}{0< \beta \le 1}.
}
\value{
 \command{\link{RMlgd}} returns an object of class
 \command{\link{RMmodel}}

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
model <- RMlgd(alpha=0.7, beta=4, scale=0.5)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}