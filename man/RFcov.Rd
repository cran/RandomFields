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
 %     fctcall = c("Cov", "CovMatrix", "Variogram", "Pseudovariogram", "Fctn"),
       ...)

RFvariogram(model, x, y = NULL, z = NULL, T = NULL, grid, distances, dim, ...)

RFpseudovariogram(model, x, y = NULL, z = NULL, T = NULL, grid, distances, dim, ...) 

RFcovmatrix(model, x, y = NULL, z = NULL, T = NULL, grid, distances, dim,...)

RFfctn(model, x, y = NULL, z = NULL, T = NULL, grid, distances, dim,...)

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
% \item{fctcall}{internal. This argument should not be considered by
% the user}
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
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\seealso{	
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
 \command{\link{RFcalc}}.
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

# locations:
x <- matrix(runif(15), ncol=3)
 # coordinate matrix of 5 arbitrary points
 # p1, p2, p3, p4, p5 in 3-dimensional space
y <- matrix(runif(15),ncol=3)
 # coordinate matrix of 5 arbitrary points
 # q1, q2, q3, q4, q5 in 3-dimensional space


#####################################################
# get available models
RFgetModelNames(type=c("tail correlation", "positive definite",
               "variogram", "shape function", "trend"))

#####################################################
# Example 1: a stationary covariance model

model <- RMexp()

# covariance only depends on differences
# of locations hi=pi-qi in 3-dimensional space
# therefore, the following 2 commands yield the same

RFcov(model=model, x=x, y=y)
RFcov(model=model, x=x-y)
 # yields 5 values C(hi)=C(pi,qi) for i=1,2,3,4,5



##################################################
# Example 2: get covariance matrix C(x_i,x_j)    
# at given locations x_i, i=1,...,n 
#
# here for an isotropic stationary covariance model
# yields a 4 times 4 covariance matrix of the form
# C(0)   C(5)   C(3)   C(2.5)
# C(5)   C(0)   C(4)   C(2.5)
# C(3)   C(4)   C(0)   C(2.5)
# C(2.5) C(2.5) C(2.5) C(0)

model <- RMexp() # the covariance function C(x,y)=C(r) of this model
#                  depends only on the distance r between x and y
RFcovmatrix(model=model, distances=c(5,3,2.5,4,2.5,2.5), dim=4)



######################################################
# Example 3: distinguish the different uses of x and y

x <- c(1,2,1)
y <- c(4,5,6)

# coordinate space 1-dimensional, evaluated at 3 points:
RFcov(model=model, x=as.matrix(x), y=as.matrix(y))
# coordinate space is 3-dimensional, evaluated at a pair of points
RFcov(model=model, x=t(x), y=t(y))



\dontshow{FinalizeExample()}

}
\keyword{spatial}






