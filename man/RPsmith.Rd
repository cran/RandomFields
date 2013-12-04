\name{Smith}
\alias{RPsmith}
\alias{mixed movin maxima}
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
   constructed via the \link[=RMstrokorbMono]{Strokorb} construction of
   deterministic, monotone functions.
 }
 \item{xi,mu,s}{the extreme value index, the location parameter and the
   scale parameter, respectively, of the generalized extreme value
   distribution. See Details.
 }

}

\note{
  \command{RPsmith} depends on \command{\link{RRrectangular}}
  and its parameters.

  Advanced options
   are \code{maxpoints} and \code{max_gauss}, see
   \command{\link{RFoptions}}.
}

\details{
  The parameter \code{xi}
  is always a number, i.e. \eqn{\xi} is constant in
 space. In contrast, \eqn{\mu} and \eqn{s} might be constant
 numerical value or given a \code{\link{RMmodel}}, in particular by a
 \code{\link{RMtrend}} model.
 

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

 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}



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
 
\examples{
set.seed(0)

model <- RMball()
x <- if (interactive()) seq(0, 1000, 0.02) else seq(0, 100, 10)
z <- RFsimulate(RPsmith(model, xi=0), x)
plot(z)
hist(z@data$variable1, 50, freq=FALSE)
curve(exp(-x) * exp(-exp(-x)), from=-3, to=8, add=TRUE)

\dontshow{if (!.C("isAuthor", a=integer(1))$a) {
}}


## for some more sophisticated models see 'maxstableAdvanced'
}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RPbernoulli}},
 \command{\link{RPgauss}},
 \link{maxstable}
 \command{\link{maxstableAdvanced}}
}




\keyword{spatial}


