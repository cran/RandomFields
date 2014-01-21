\name{RMdagum}
\alias{RMdagum}
\title{Dagum Covariance Model Family}
\description{
 \command{\link{RMdagum}} is a stationary isotropic covariance model.
 The corresponding covariance function only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r) = 1-(1+r^{-\beta})^{\frac{-\gamma}{\beta}}.}{C(r) = 1-(1+r^(-\beta))^(-\gamma / \beta).}
 The parameters \eqn{\beta}{\beta} and \eqn{\gamma}{\gamma} can be varied in the intervals \eqn{(0,1]}{(0,1]}
and \eqn{(0,1)}{(0,1)}, respectively. 
}
\usage{
RMdagum(beta, gamma, var, scale, Aniso, proj)
}
\arguments{
 \item{beta}{a parameter in \eqn{(0,1]}{(0,1]}}
 \item{gamma}{a parameter in \eqn{(0,1)}{(0,1)}}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{Like the generalized Cauchy model the Dagum family can be used
 to model RFfractaldimension and Hurst effect. For a comparison of these
 see Berg, C. and Mateau, J. and Porcu, E. (2008). This paper also
 establishes valid parameter choices for the Dagum family, but be
 careful because therein the model is parameterized differently.
}
\value{
 \command{\link{RMdagum}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Berg, C. and Mateau, J. and Porcu, E. (2008)
 The dagum family of isotropic correlation functions.
 \emph{Bernoulli} \bold{14}(4), 1134--1149.
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
RFoptions(seed=0)
model <- RMdagum(beta=0.5, gamma=0.5, scale=0.2)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{RFoptions(seed=NA)}
}