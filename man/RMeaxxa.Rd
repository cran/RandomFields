\name{RMeaxxa}
\alias{RMeaxxa}
\alias{RMetaxxa}
\title{Special models for rotation like fields}
\description{
  \command{RMeaxxa} and \command{RMetaxxa} define the auxiliary functions
 \deqn{f(h) = h^\top A A^\top h + diag(E)}{C(h) = h^T A A^T h + diag(E)}
 and
 \deqn{f(h) = h^\top A R R A^\top h + diag(E)}{C(h) = h^T A R R^R A^T h
   + diag(E),}
 respectively.
 }
 \usage{
RMeaxxa(E, A) 
RMetaxxa(E, A, alpha)
}

\arguments{
 \item{E}{m-variate vector of positive values}
 \item{A}{\eqn{m\times k}{m x k} matrix}
 \item{alpha}{angle for the rotation matrix \eqn{R}}
}

\details{
 \command{\link{RMeaxxa}} is defined in space and returns an m-variate model.

 \command{\link{RMetaxxa}} is a space-time model with two spatial dimensions.
 The matrix R is a rotation matrix with angle \eqn{\beta t} where
 \eqn{t} is the time component.
}

 \value{
 \command{\link{RMeaxxa}} and \command{\link{RMetaxxa}} return an
 object of class \code{\link[=RMmodel-class]{RMmodel}}.

}

\references{
 \itemize{
 \item Schlather, M. (2010)
 On some covariance models based on normal scale mixtures.
 \emph{Bernoulli}, \bold{16}, 780-797. 
 }
}

\me
\seealso{
 \command{\link{RMmodel}},
  \link{S10}
 }

\keyword{spatial}
\keyword{models}


\examples{\dontshow{StartExample()}
# see S10
\dontshow{FinalizeExample()}}

