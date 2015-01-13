\name{RMcutoff}
\alias{RMcutoff}
\title{Gneiting's modification towards finite range}
\description{
 \command{\link{RMcutoff}} is a functional on univariate stationary
 isotropic covariance functions \eqn{\phi}{phi}. 

 The corresponding function \eqn{C} (which is not necessarily a
 covariance function, 
 see details) only depends on the distance \eqn{r}{r} between two
 points in \eqn{d}-dimensional space and is given by

 \deqn{C(r)=\phi(r), 0\le r \le d}
 \deqn{C(r) = b_0 ((dR)^a - r^a)^{2 a}, d \le r \le dR}
 \deqn{C(r) = 0, dR \le r}
 
 The parameters \eqn{R} and \eqn{b_0}
 are chosen internally such that \eqn{C} is a smooth function.
}

\usage{
RMcutoff(phi, diameter, a, var, scale, Aniso, proj)
}
\arguments{
  \item{phi}{a univariate stationary isotropic covariance model.    
    See, for instance,

    \code{RFgetModelNames(type="positive definite",
      domain="single variable", isotropy="isotropic",  vdim=1)}.
 }
 \item{diameter}{a numerical value; should be greater than 0; the
   diameter of the domain on which the simulation is done}
 \item{a}{a numerical value; should be greater than 0; has been shown to be
   optimal for \eqn{a = 1/2} or \eqn{a =1}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.}
}
\details{The algorithm that checks the given parameters knows
 only about some few necessary conditions.
 Hence it is not ensured that
 the cutoff-model is a valid covariance function for any
 choice of \eqn{\phi} and the parameters. 
 
 For certain models \eqn{\phi}{phi}, e.g. \command{\link{RMstable}},
 \command{\link{RMwhittle}} and \command{\link{RMgencauchy}}, some
 sufficient conditions 
 are known (cf. Gneiting et al. (2006)).
}
\value{
 \command{\link{RMcutoff}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
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
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
 ## For examples see the help page of 'RFsimulateAdvanced' ##
model <- RMexp()
plot(model, model.cutoff=RMcutoff(model, diameter=1), xlim=c(0, 4))

\dontshow{FinalizeExample()}
}