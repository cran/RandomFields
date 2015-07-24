
\name{Square roots}
\alias{Direct}
\alias{RPdirect}
\alias{Sequential}
\alias{RPsequential}
\title{Methods relying on square roots of the covariance matrix}



\description{
 Methods relying on square roots of the covariance matrix
}

\usage{
RPdirect(phi, boxcox, root_method, svdtolerance, max_variab) 

RPsequential(phi, boxcox, max_variables, back_steps, initial)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 specifies the covariance model to be simulated.}
\item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  see \command{\link{RFboxcox}} for Details.
 }

% \item{loggauss}{optional argument; same meaning as for 
% \command{\link{RPgauss}}.}

\item{root_method}{Decomposition of the covariance matrix:
  \itemize{
    \item  \code{root_method=0} Cholesky decomposition
    \item  \code{root_method=2} SVD
    }
    If Cholesky decomposition failes then alway SVD is tried.

    Default: \code{0} .
    }
    \item{svdtolerance}{ If SVD decomposition is used for calculating the square root of
 the covariance matrix then the absolute componentwise difference between
 the covariance matrix and square of the square root must be less
 than \code{svdtolerance}. No check is performed if
 \code{svdtolerance} is negative.
 
 Default: \code{1e-12} .
}
 \item{max_variab}{If the number of variables to generate is
 greater than \code{maxvariables}, then any matrix decomposition
 method is rejected. It is important that this option is set
 conveniently to avoid great losses of time during the automatic
 search of a simulation method (\code{method="any method"}). 

 Default: \code{8192} 
}

 \item{max_variables}{The maximum size of the conditional covariance matrix
 (default to 5000)}
 \item{back_steps}{
   Number of previous instances on which
   the algorithm should condition.
   If less than one then the number of previous instances
   equals \code{max} / (number of spatial points).
   
   Default: \code{10} .
   
 }
 \item{initial}{
   First, N=(number of spatial points) * \code{back_steps}
   number of points are simulated. Then, sequentially,
   all spatial points for the next time instance
   are simulated at once, based on the previous \code{back_steps}
   instances. The distribution of the first N points
   is the correct distribution, but
   differs, in general, from the distribution of the sequentially
   simulated variables. We prefer here to have the same distribution
   all over (although only approximatively the correct one),
   hence do some initial sequential steps first.
   If \code{initial} is non-negative, then \code{initial}
   first steps are performed.
   If \code{initial} is negative, then
   \code{back_steps} - \code{initial}
 initial steps are performed. The latter ensures that
 none of the very first N variables are returned.
 
 Default: \code{-10} .
 }
}


\details{
 \command{RPdirect}
 is based on the well-known method for simulating
 any multivariate Gaussian distribution, using the square root of the
 covariance matrix. The method is pretty slow and limited to
 about 8000 points, i.e. a 20x20x20 grid in three dimensions. 
 This implementation can use the Cholesky decomposition and
 the singular value decomposition. 
 It allows for arbitrary points and arbitrary grids.
 
 \command{RPsequential}
 is programmed for spatio-temporal models
 where the field is modelled sequentially in the time direction
 conditioned on the previous \eqn{k} instances.
 For \eqn{k=5} the method has its limits for about 1000 spatial
 points. It is an approximative method. The larger \eqn{k} the
 better.
 It also works for certain grids where the last dimension should
 contain the highest number of grid points. 
}

\value{
 \command{\link{RPsequential}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

}

\references{
 \itemize{
 \item
 Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University. 

 }}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\seealso{
  \link{RP},
 \command{\link{RPcoins}},
 \command{\link{RPhyperplane}},
 \command{\link{RPspectral}},
 \command{\link{RPtbm}}.
}

\keyword{methods}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMgauss(var=10, s=10) + RMnugget(var=0.01)
plot(model, xlim=c(-25, 25))

z <- RFsimulate(model=RPdirect(model), 0:10, 0:10, n=4)
plot(z)

\dontshow{FinalizeExample()}
}
