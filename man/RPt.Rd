\name{RPt}
\alias{RPt}

\title{Simulation of T Random Fields}
\description{
  \command{RPt} defines a t field.
}

\usage{

RPt(phi, boxcox, nu)

}

\arguments{
 \item{phi}{the \command{\link{RMmodel}}. If a model for the
 distribution is not specified, \command{\link{RPgauss}} is used as
 default and a covariance model is expected.}
\item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  See \command{\link{RFboxcox}} for details.
 }
 \item{nu}{non-negative number. Degree of freedom.}
}

\value{
 The function returns an object of class \code{\link[=RMmodel]{RMmodel}}.
}

\me

\seealso{
  \command{\link{Auxiliary RMmodels}},
  \command{\link{RP}},
  \command{\link{RPgauss}}.
}

\references{
  Related to the extremal t process
  \itemize{
    \item
    T. Opitz (2012) A spectral construction of the extremal t process.
    \emph{arxiv} \bold{1207.2296}.
  }
}


\keyword{spatial}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RPt(RMexp(), nu=2)
x <- seq(0, 10, 0.1)
z <- RFsimulate(model, x, x, n=4)
plot(z)

\dontshow{FinalizeExample()}}
