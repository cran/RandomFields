\name{RMexponential}
\alias{RMexponential}
\alias{RMexponential}
\title{Exponential operator}
\description{
 \command{\link{RMexponential}} yields a covariance model 
 from a given variogram or covariance model.
 The covariance \eqn{C} is given as
 \deqn{
 C(h) = \frac{\exp(\phi(h)) -\sum_{k=0}^n \phi^k(h)/k!}{\exp(\phi(0))
 -\sum_{k=0}^n \phi^k(0)/k!}}{ 
 C(h) = (\exp(\phi(h)) -\sum_{k=0}^n \phi^k(h)/k!) / (\exp(\phi(0))
 -\sum_{k=0}^n \phi^k(0)/k!)}
 if \eqn{\phi} is a covariance model, and as
 \deqn{ C(h) = \exp(-\phi(h))}
 if \eqn{\phi} is a variogram model.
}
\usage{
RMexponential(phi, n, standardised, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a valid \command{\link{RMmodel}}; either a variogram model or a
 covariance model}
 \item{n}{integer, see formula above. Default is -1.; if the
 multivariate dimension of the submodel is greater than 1 then only the
 default value is valid.}
\item{standardised}{logical. If \code{TRUE} then the above formula
  holds. If \code{FALSE} then only the nominator of the above formula is
  returned. Default value is \code{TRUE}.
}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any \command{\link{RMmodel}}. If not passed, the above covariance function remains unmodified.}
}
\details{
 If \eqn{\gamma} is a variogram, then \eqn{\exp(-\gamma)} is a valid
 covariance. 
}
\value{
 \command{\link{RMexponential}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 See, for instance,

 \itemize{
   \item
   Berg, C., Christensen, J. P. R., Ressel, P. (1984)
   Harmonic Analysis on Semigroups. \emph{Theory of Positive Definite and
     Related Functions.} Springer, New York.

   \item
   Sasvari, Z. (2013) \emph{Multivariate Characteristic and Correlation
     Functions}. de Gruyter, Berlin.
   
   \item Schlather, M. (2010) \emph{Some covariance models based on
     normal scale mixtures}, \emph{Bernoulli} \bold{16}, 780-797.
   
   \item Schlather, M. (2012) Construction of covariance functions and
   unconditional simulation of random fields. In
   Porcu, E., Montero, J. M., Schlather, M. 
   \emph{Advances and Challenges in Space-time Modelling of Natural Events},
   Springer, New York.
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
model <- RMexponential(RMfbm(alpha=1))  ## identical to RMexp()
plot(RMexp(), model=model, type=c("p", "l"), pch=20) 

\dontshow{FinalizeExample()}
}