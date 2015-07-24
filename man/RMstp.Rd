\name{RMstp}
\alias{RMstp}
\alias{RMstp}
\title{Single temporal process}
\description{
  \command{\link{RMstp}} is a univariate covariance model which depends on
  a normal mixture submodel \eqn{\phi}{phi}.
  The covariance is given by
  \deqn{ 
    C(x,y) = |S_x|^{1/4} |S_y|^{1/4} |A|^{-1/2}
    \phi(Q(x,y)^{1/2})
  }
  where
  \deqn{
    Q(x,y) = c^2 - m^2 + h^t (S_x + 2(m + c)M) A^{-1} (S_y + 2
    (m-c)M)h,
  }
  \deqn{
    c = -z^t h + \xi_2(x) - \xi_2(y),
  }
  \deqn{
    A = S_x + S_y + 4 M h h^t M 
  }
  \deqn{
    m = h^t M h
  }
  \deqn{
    h = x - y
  }
}
\usage{
RMstp(xi, phi, S, z, M, var, scale, Aniso, proj) 
}
\arguments{
 \item{xi}{arbitrary univariate function on \eqn{R^d}}
 \item{phi}{an \command{\link{RMmodel}} that is a normal mixture model,
   cf.\cr
    \code{RFgetModelNames(monotone="normal mixture")}
 }
 \item{S}{functions that returns strictly positive definite
   \eqn{d\times d}{d x d}}
 \item{z}{arbitrary vector, \eqn{z \in R^d}}
 \item{M}{an arbitrary, symmetric \eqn{d \times d}{d x d} matrix}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any \command{\link{RMmodel}}. If not passed, the above covariance function remains unmodified.}
}
\details{
 See Schlather (2008) formula (13).
 The model allows for mimicking cyclonic behaviour.
}
\value{
 \command{\link{RMstp}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
  \itemize{
 
    \item Paciorek C.J., and Schervish, M.J. (2006)
    Spatial modelling using a new class of nonstationary covariance functions,
    \emph{Environmetrics} \bold{17}, 483-506.
    
    \item Schlather, M. (2010)
    Some covariance models based on normal scale mixtures.
    \emph{Bernoulli}, \bold{16}, 780-797.
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
\dontshow{StartExample()}
model <- RMstp(xi = RMrotat(phi= -2 * pi, speed=1),
               phi = RMwhittle(nu = 1),
               M=matrix(nc=3, rep(0, 9)),
               S=RMetaxxa(E=rep(1, 3), alpha = -2 * pi,
                          A=t(matrix(nc=3, c(2, 0, 0, 1, 1 , 0, 0, 0, 0))))
              )
x <- seq(0, 10, 0.7)
plot(RFsimulate(model, x=x, y=x, z=x))
\dontshow{FinalizeExample()}
}
