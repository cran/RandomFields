\name{RMmppplus}
\alias{RMmppplus}
\alias{++}
\title{Mixture of shape functions}
\description{
 \command{\link{RMmppplus}} is a multivariate covariance model which depends on
 up to 10 submodels \eqn{C_0, C_1, ..., C_10}.

 It is used to gether with \command{\link{RPsmith}} it allowed for mixed
 moving maxima with a finite number of shape functions.
 
}
\usage{
RMmppplus(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, p)%, var,scale,Aniso,proj)
}
\arguments{
 \item{C0}{an \command{\link{RMmodel}}.}
 \item{C1,C2,C3,C4,C5,C6,C7,C8,C9}{optional; each an \command{\link{RMmodel}}.}
 \item{p}{vector of probabilities for the shape functions.
   The probabilities should add up to 1. The length
   of the vector equals the number of given submodels}
% \item{var,scale,Aniso,proj}{optional arguments; same meaning for any \command{\link{RMmodel}}. If not passed, the above model remains unmodified.}
}
\value{
 \command{\link{RMmppplus}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMplus}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RPsmith}}.
}

\keyword{spatial}
\keyword{models}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
\dontshow{\dontrun{
model <- RPcoins(RMmppplus(p=c(0.8, 0.2), RMgauss(), RMexp()), intensity=1)
x <- seq(0, 10, 0.02)
plot(RFsimulate(model, x=x))
}}
\dontshow{FinalizeExample()}
}
