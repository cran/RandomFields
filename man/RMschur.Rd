\name{RMschur}
\alias{RMschur}
\title{Schur product}
\description{
 The covariance function is

 \deqn{C(x)= M * \phi(x)}

 where \sQuote{*} denotes the Schur product, i.e. elementwise
 multiplication
 
 }
 \usage{
RMschur(phi, M,  diag, rhored, var, scale, Aniso, proj)
}

\arguments{
 \item{phi}{covariance function of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{M}{constant \eqn{n \times n}{n x n} covariance matrix of the
   same size as multivariate model \code{phi}}
 \item{diag,rhored}{alternative way of passing \code{M}:
   \code{diag} is a vector of variances, \code{rhored} is a vector
   containing the correlations of lower triangle of the \code{M}
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}

\value{
 \command{\link{RMschur}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
\itemize{
 \item ?
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RMmatrix}},
 }


\keyword{spatial}
\keyword{models}




\examples{
 ## For an example see the help page of 'RMparswmX' ##
\dontshow{FinalizeExample()}
}