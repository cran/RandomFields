\name{RMchoquet}
\alias{RMchoquet}
\alias{Choquet's representation}
\alias{Schoenberg's representation}
\title{Schoenberg's representation for the classes psi_d and \eqn{psi_{\infty}} in d=2}
\description{
  \command{\link{RMchoquet}} is an isotropic covariance model. The
  corresponding covariance function only depends on the angle
  \eqn{0 \le \theta \le \pi}{0 \le \theta \le \pi}
  between two points on the sphere and is given for d=2 by 
 \deqn{\psi(\theta) = \sum_{n=0}^{\infty} b_{n,2}/(n+1)*P_n(cos(\theta)),}
 where \deqn{\sum_{n=0}^{\infty} b_{n,d}=1}
   and \eqn{P_n} is the Legendre Polynomial of integer order \eqn{n >= 0}{n >= 0}.
}

\usage{
RMchoquet(b)
}

\arguments{
  \item{b}{a numerical vector of weights in \eqn{(0,1)}, such that sum(b)=1.}
}

\details{
  By the results (cf. Gneiting, T. (2013), p.1333) of Schoenberg and
  others like Menegatto, Chen, Sun,  Oliveira and Peron, the class
  \eqn{psi_d} of all real valued funcions on \eqn{[0,\pi]}, with
  \eqn{\psi(0)=1} and such that the associated isotropic function
  \deqn{h(x,y)=\psi(theta)     with     cos(\theta)=<x,y>}
  \deqn{for x,y in  {x in R^d: ||x|| = 1}}
  is (strictly) positive definite is represented by this
  covariance model. The model can be interpreted as Choquet
  representation in terms of extremal members, which are non-strictly
  positive definite.

  Special cases are the multiquadric family (see
  \command{\link{RMmultiquad}}) and the model of the sine power function (see
  \command{\link{RMsinepower}}).

}


\value{
 \command{\link{RMchoquet}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
  \itemize{
    \item Gneiting, T. (2013)
    \emph{Strictly and non-strictly positive definite functions on
      spheres.} Bernoulli, \bold{19}(4), 1327-1349.

  \item Schoenberg, I.J. (1942)
  \emph{Positive definite functions on spheres.} Duke Math.J.,\bold{9}, 96-108.

  \item Menegatto, V.A. (1994)
  \emph{Strictly positive definite kernels on the Hilbert sphere.}
  Appl. Anal., \bold{55}, 91-101.

  \item Chen, D., Menegatto, V.A., and Sun, X. (2003)
  \emph{A necessary and sufficient condition for strictly positive
    definite functions on spheres.}
  Proc. Amer. Math. Soc.,\bold{131}, 2733-2740.

  \item Menegatto, V.A., Oliveira, C.P. and Peron, A.P. (2006)
  \emph{Strictly positive definite kernels on subsets of the complex
    plane.} Comput. Math. Appl., \bold{51}, 1233-1250.
  }
}

\author{Christoph Berreth; \martin}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{spherical models}},
 \command{\link{RMmultiquad}},
 \command{\link{RMsinepower}}

}

\keyword{spatial}
\keyword{models}


\examples{\dontshow{StartExample()}

RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## to do

\dontshow{\dontrun{
#b = 
#model <- RMchoquet(b=b)
#x <- seq(0, 10, 0.02)
#plot(model)
#plot(RFsimulate(model, x=x))
}}
\dontshow{FinalizeExample()}}