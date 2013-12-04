\name{RFcov}
\alias{RFcov}
\alias{RFvariogram}
\alias{RFpseudovariogram}
\alias{RFcovmatrix}
\alias{RFfctn}
\title{Evaluate Covariance and Variogram Functions}
\description{
 \command{\link{RFcov}} returns the values of a covariance function;
 \command{\link{RFvariogram}} returns the values of a variogram;
 \command{\link{RFpseudovariogram}} returns the values of a pseudovariogram;
 \command{\link{RFcovmatrix}} returns the covariance matrix for a set of points
 \command{\link{RFfctn}} returns the values of a shape function; 
}
\usage{
RFcov(model, x, y = NULL, z = NULL, T = NULL, grid, distances, dim,
      fctcall = c("Cov", "CovMatrix", "Variogram", "Pseudovariogram",
                  "Fctn"),
       ...)

RFvariogram(model, x, y = NULL, ...)

RFpseudovariogram(model, x, y = NULL, ...) 

RFcovmatrix(...)

RFfctn(...)

}
\arguments{
 \item{model}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 the covariance or variogram model, which is to be evaluated}
 \item{x}{vector or \eqn{(n \times \code{dim})}{(n x
 \code{dim})}-matrix, where \eqn{n} is the number of points at
 which the covariance function is to be evaluated;
 in particular,
 if the model is isotropic or \code{dim=1} then \code{x}
 is a vector. \code{x}}
 \item{y}{second vector or matrix for non-stationary covariance
 functions}
 \item{z}{z-component of point if xyzT-specification of points is used}
 \item{T}{T-component of point if xyzT-specification of points is used}
 \item{grid}{boolean; whether xyzT specify a grid}
 \item{distances}{vector; only if the function \command{\link{RFcovmatrix}} is
 used; the lower triangular part of the distance matrix column-wise;
 equivalently the upper triangular part of the distance matrix row-wise;
 either \code{x} or \code{distances} must be missing}
 \item{dim}{dimension of the coordinate space in which the model is applied}
 \item{fctcall}{internal. This parameter should not be considered by
 the user}
 \item{...}{arguments passed to \command{\link{RFcov}} (\command{\link{RFcovmatrix}})
 and arguments passed to \command{\link{RFoptions}}}
}
\details{
 \command{\link{RFcovmatrix}} returns a covariance matrix. Here a matrix of
 of coordinates (\code{x}) or a vector or a matrix of \code{distances}
 is expected.
 \command{\link{RFcovmatrix}} allows also for variogram models. Then the negative of the
 variogram matrix is returned.
}

\value{
 \command{\link{RFcov}} returns a vector of values of the covariance
 function.
 
 \command{\link{RFvariogram}} returns a vector of values of the variogram model.

 \command{\link{RFpseudovariogram}} returns a vector of values of the variogram model.

 \command{\link{RFcovmatrix}} returns a covariance matrix. 
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\seealso{	
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}

\examples{
set.seed(0)

# locations:
X <- matrix(runif(15), ncol=3)
 # coordinate matrix of 5 arbitrary points
 # p1, p2, p3, p4, p5 in 3-dimensional space
Y <- matrix(runif(15),ncol=3)
 # coordinate matrix of 5 arbitrary points
 # q1, q2, q3, q4, q5 in 3-dimensional space


#####################################################
# evaluate covariance function C(x,y) #
# or variogram gamma(x,y) #
# at given locations #
#####################################################

RFgetModelNames() # returns available models


#####################################################
# a stationary covariance model

model <- RMexp()

# covariance only depends on differences
# of locations hi=pi-qi in 3-dimensional space
# therefore, the following 2 commands yield the same

RFcov(model=model, x=X, y=Y)
RFcov(model=model, x=X-Y)
 # yields 5 values C(hi)=C(pi,qi) for i=1,2,3,4,5

#####################################################
# get covariance matrix C(x_i,x_j) #
# at given locations x_i, i=1,...,n #
#####################################################


#####################################################
# in particular compare the following two evaluations

x <- c(1,2,1)
y <- c(4,5,6)
RFcov(model=model,
 x=as.matrix(x), y=as.matrix(y), grid=FALSE)
 # coord space 1-dim'l
 # evaluate at 3 points
RFcov(model=model,
 x=t(x),y=t(y)) # coord space 3-dim'l
 # evaluate at 1 point

#####################################################
# an isotropic stationary covariance model

model <- RMexp()
# the covariance function C(x,y)=C(r) of this model
# depends only on the distance r between x and y

RFcovmatrix(model=model, distances=c(5,3,2.5,4,2.5,2.5),
 dim=4)

## yields a 4 times 4 covariance matrix of the form
## C(0) C(5) C(3) C(2.5)
## C(5) C(0) C(4) C(2.5)
## C(3) C(4) C(0) C(2.5)
## C(2.5) C(2.5) C(2.5) C(0)

#####################################################
# a variogram model

model <- RMfbm(alpha=1)
 # fractal Brownian Motion is a variogram model

RFvariogram(model=model, x=X, y=Y)
 # yields 5 values gamma(pi,qi) for i=1,2,3,4,5
RFvariogram(model=model, x=X-Y)
 # gives same output

\dontshow{\dontrun{
#####################################################
# a non-isotropic covariance model

model <- RMnonstwm(nu=RMexp())
C <- RFcovmatrix(model=model, x=X)
# yields covariance matrix C(xi,xj) for i,j=1,2,3,4,5
C
}}

}
\keyword{spatial}






