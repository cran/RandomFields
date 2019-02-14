\name{Smith}
\alias{RPsmith}
\alias{mixed moving maxima}
\alias{moving maxima}
\alias{M2}
\alias{M3}

\title{(Mixed) Moving Maxima}
\description{
  \command{RPsmith} defines a moving maximum process or a mixed moving
  maximum process with finite number of shape functions.
}

\usage{
RPsmith(shape, tcf, xi, mu, s)
}

\arguments{
 \item{shape}{an \command{\link{RMmodel}} giving the spectral function}
 \item{tcf}{an \command{\link{RMmodel}} specifying the
   extremal correlation function; either \code{shape} or \code{tcf} must
   be given. If \code{tcf} is given a shape function is tried to be 
   constructed via the \command{\link{RMm2r}} construction of
   deterministic, monotone functions.
 }
 \item{xi,mu,s}{the extreme value index, the location parameter and the
   scale parameter, respectively, of the generalized extreme value
   distribution. See Details.
 }

}

\note{
  IMPORTANT: For consistency reasons with the geostatistical definitions
  in this package the scale argument differs froms the original
  definition of the Smith model! See the example below.
  
  \command{RPsmith} depends on \command{\link{RRrectangular}}
  and its arguments.

  Advanced options
   are \code{maxpoints} and \code{max_gauss}, see
   \command{\link{RFoptions}}.
}

\details{
  \GEV
 
 It simulates max-stable processes \eqn{Z} that are referred to as
 \dQuote{Smith model}.
 \deqn{Z(x) = \max_{i=1}^\infty X_i Y_i(x-W_i),
 }{Z(x) = max_{i=1, 2, ...} X_i * Y_i(x - W_i),}
 where \eqn{(W_i, X_i)} are the points of a Poisson point process on 
 \eqn{\R^d \times (0, \infty)}{R^d x (0, \infty)} with intensity
 \eqn{dw * c/x^2 dx} and \eqn{Y_i \sim Y}{Y_i ~ Y} are iid measurable
 random functions with
 \eqn{E[\int \max(0, Y(x)) dx] < \infty}{E[int max(0, Y(x)) dx ] < \infty}. 
 The constant \eqn{c} is chosen such that \eqn{Z} has standard Frechet
 margins.
 }
 
\me

\references{
  \itemize{
    \item Haan, L. (1984) 
    A spectral representation for max-stable processes.
    \emph{Ann. Probab.}, \bold{12},
    1194-1204.
    
    \item Smith, R.L. (1990) Max-stable processes and spatial extremes
    Unpublished Manuscript.
    }
}

\seealso{
  \command{\link{Advanced RMmodels}},
  \command{\link{Auxiliary RMmodels}},
  \command{\link{RMmodel}},
  \command{\link{RPbernoulli}},
  \command{\link{RPgauss}},
  \link{maxstable},
  \command{\link{maxstableAdvanced}}.
}

\keyword{spatial}

 
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMball()
x <- seq(0, 1000, 0.2)
z <- RFsimulate(RPsmith(model, xi=0), x)
plot(z)
hist(z@data$variable1, 50, freq=FALSE)
curve(exp(-x) * exp(-exp(-x)), from=-3, to=8, add=TRUE)

## 2-dim
x <- seq(0, 10, 0.1) 
z <- RFsimulate(RPsmith(model, xi=0), x, x)
plot(z)

## original Smith model
x <- seq(0, 10, 0.05)
model <- RMgauss(scale = sqrt(2)) # !! cf. definition of RMgauss
z <- RFsimulate(RPsmith(model, xi=0), x, x)
plot(z)

\dontshow{\dontrun{
## original Smith model at arbitrary locations
n <- 200
model <- RMgauss(scale = sqrt(2)) # !! cf. definition of RMgauss
z <- RFsimulate(RPsmith(model, xi=0), runif(n, 0, 10), runif(0, 0, 10))
plot(z)
}}


## for some more sophisticated models see 'maxstableAdvanced'

\dontshow{FinalizeExample()}}

