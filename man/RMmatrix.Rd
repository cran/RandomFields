\name{RMmatrix}
\alias{RMmatrix}
\title{Matrix operator}
\description{
 \command{\link{RMmatrix}} is a multivariate covariance model
 depending on a multivariate covariance model \eqn{\phi}{phi}.
 The corresponding covariance function is given by
 \deqn{ C(h) = M \phi(h) M^t}{ C(h) = M phi(h) M^t}
}
\usage{
RMmatrix(phi, M, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a k-variate covariance \command{\link{RMmodel}}.}
 \item{M}{a k times k matrix}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
%\details{}
\value{
 \command{\link{RMmatrix}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
  Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and
    Strokorb, K. (2015) 
    Analysis, simulation and prediction of multivariate
    random fields with package \pkg{RandomFields}. \emph{
      Journal of Statistical Software}, \bold{63} (8), 1-25,
    url          = \sQuote{http://www.jstatsoft.org/v63/i08/}
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
# bivariate Linear Model of Coregionalisation
model <- RMmatrix(M = c(0.9, 0.43), RMwhittle(nu = 0.3)) + 
         RMmatrix(M = c(0.6, 0.8), RMwhittle(nu = 2))
x <- y <- seq(-10, 10, 0.2)
simu <- RFsimulate(model, x, y)
plot(simu)

\dontshow{FinalizeExample()}

}