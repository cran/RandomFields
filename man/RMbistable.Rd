\name{RMbistable}
\alias{RMbistable}
\title{Bivariate stableModel}
\description{
 \command{\link{RMbistable}} is a bivariate stationary isotropic covariance model
 whose corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between
 two points
% and is given for \eqn{i,j \in \{1,2\}}{i,j = 1,2} by
% \deqn{C_{ij}(r)=[to be done].}

 For constraints on the constants see details.
}

\usage{
RMbistable(alpha, s, cdiag, rho, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{[to be done]}
 \item{s}{a vector of length 3 of numerical values; each entry
 positive; the vector \eqn{(s_{11},s_{21},s_{22})}}
 \item{cdiag}{[to be done]}
 \item{rho}{[to be done]
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 Constraints on the constants: [to be done] 
 }
\value{
 \command{\link{RMbistable}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Moreva, O., Schlather, M. (2016) Modelling and simulation of bivariate Gaussian random fields.
 \emph{arXiv 1412.1914}
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMstable}},
 \link{Multivariate RMmodels},
}


\keyword{spatial}
\keyword{models}



\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
## todo

\dontshow{FinalizeExample()}
}