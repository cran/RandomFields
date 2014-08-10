\name{RMintrinsic}
\alias{RMintrinsic}
\title{Intrinsic Embedding Covariance Model}
\description{
 \command{\link{RMintrinsic}} is a univariate stationary isotropic covariance
 model which depends on a univariate stationary isotropic covariance model.
 
 The corresponding covariance function C of the model
 only depends on the distance \eqn{r \ge 0}{r \ge 0} between
 two points and is given by
 \deqn{C(r)=a_0 + a_2 r^2 + \phi(r), 0\le r \le diameter}
 \deqn{C(r)=b_0 (rawR D - r)^3/(r), diameter \le r \le rawR * diameter}
 \deqn{C(r) = 0, rawR * diameter \le r}
}
\usage{
RMintrinsic(phi, diameter, rawR, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{an \command{\link{RMmodel}}; has to be stationary and isotropic}
 \item{diameter}{a numerical value; positive; should be the diameter of
 the domain on which simulation is done}
 \item{rawR}{a numerical value; greater or equal to 1}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 The parameters \eqn{a_0}, \eqn{a_2} and \eqn{b_0}
 are chosen internally such that \eqn{C} becomes a smooth function.
 See formulas (3.8)-(3.10) in Gneiting et alii (2006).
 This model corresponds to the method Intrinsic Embedding.
 See also \code{\link{RPintrinsic}}.

 NOTE: The algorithm that checks the given parameters knows
 only about some few necessary conditions.
 Hence it is not ensured that
 the Stein-model is a valid covariance function for any
 choice of \eqn{\phi} and the parameters. 
 
 For certain models \eqn{\phi}{phi}, i.e. \code{stable},
 \code{whittle}, \code{gencauchy}, and the variogram
 model \code{fractalB} some sufficient conditions are known.
}
\value{
 \command{\link{RMintrinsic}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Gneiting, T., Sevecikova, H, Percival, D.B., Schlather M.,
 Jiang Y. (2006) Fast and Exact Simulation of Large {G}aussian
 Lattice Systems in {$R^2$}: Exploring the Limits.
 \emph{J. Comput. Graph. Stat.} \bold{15}, 483--501.
 \item Stein, M.L. (2002) Fast and exact simulation of fractional
 Brownian surfaces. \emph{J. Comput. Graph. Statist.} \bold{11}, 587--599
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RPintrinsic}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
x.max <- 10
model <- RMintrinsic(RMfbm(alpha=1), diameter=x.max)
x <- seq(0, x.max, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}