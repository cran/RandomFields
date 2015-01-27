\name{RMiaco}
\alias{RMiaco}
\title{Iaco-Cesare model}
\description{
 The space-time covariance function is
 \deqn{
 C(r,t) = (1.0 + r^\nu + t^\lambda)^\delta
 }
 
 }
 \usage{
RMiaco(nu, lambda, delta, var, scale, Aniso, proj)
}

\arguments{
 \item{nu,lambda}{number in \eqn{(0,2]}}
 \item{delta}{positive number}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}


\value{
 \command{\link{RMiaco}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
\itemize{
 \item de Cesare, L., Myers, D.E., and Posa, D. (2002)
 FORPRAN programs for space-time modeling. \emph{Computers \&
 Geosciences} \bold{28}, 205-212.
 \item de Iaco, S.. Myers, D.E., and Posa, D. (2002)
 Nonseparable space-time covariance models: some parameteric
 families. \emph{Math. Geol.} \bold{34}, 23-42.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 }


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMiaco(nu=1, lambda=1.5, delta=0.5)
plot(model, dim=2)

x <- seq(0, 10, if (interactive()) 0.1 else 5) 
plot(RFsimulate(model, x=x, y=x))
\dontshow{FinalizeExample()}
}