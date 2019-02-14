\name{RMlsfbm}
\alias{RMlsfbm}
\title{Locally Positive Definite Function Given by the Fractal Brownian Motion}
\description{
  \command{\link{RMlsfbm}} is a positive definite function on the
  unit ball in \eqn{R^d} centred at the origin,
  \deqn{C(r) = c - r^\alpha}{C(r) = c - r^\alpha}
  with \eqn{r = \|x- y\|\in [0,1]}{0 <= r = || x - y || <= 1}.
}
\usage{
RMlsfbm(alpha, const, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{numeric in \eqn{(0,2)}; refers to the fractal dimension of the
   process.}
 \item{const}{the \code{const}ant \eqn{c} is given by the
   formula
   \deqn{
     c = 2^{-\alpha} \Gamma(d / 2 + \alpha/2) \Gamma(1 - \alpha/2) /
     \Gamma(d / 2)
   }
   and should not be changed by the user in order to ensure positive
   definiteness.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.}
}

\value{
 \command{\link{RMlsfbm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
   \item Martini, J., Schlather, M., Simianer, H. (In preparation.)
  }
}

\me

\seealso{
 \command{\link{RMbcw}} generalizes \command{RMlsfbm} in case that \eqn{c}
 is given,
 \command{\link{RMfbm}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMlsfbm(alpha=1, scale=10)
x <- seq(0, 10, 0.02)
plot(model, xlim=c(0,10))
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}}