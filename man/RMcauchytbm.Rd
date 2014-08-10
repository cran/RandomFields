\name{RMcauchytbm}
\alias{RMcauchytbm}
\title{Modifications of the Cauchy Family Covariance Model}
\description{
 \code{\link{RMcauchytbm}()} is a shortcut of
 \command{\link{RMtbm}(\link{RMgencauchy}())} and is 
 given here for downwards compatibility
}
\usage{
RMcauchytbm(alpha, beta, gamma, var, scale, Aniso, proj)
}

\arguments{
 \item{alpha,beta}{see \command{\link{RMgencauchy}}.}
 \item{gamma}{is the same as \code{fulldim} in \command{\link{RMtbm}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}

 
\value{
 \command{\link{RMcauchytbm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

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
 \command{\link{RMcauchy}},
 \command{\link{RMgencauchy}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMcauchytbm(alpha=1, beta=1, gamma=3)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}