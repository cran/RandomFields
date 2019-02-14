\name{RFcovmatrix}
\alias{RFcovmatrix}
\title{Covariance matrix}
\description{
 \command{\link{RFcovmatrix}} returns the covariance matrix for a set of points;
 }
\usage{
RFcovmatrix(model, x, y = NULL, z = NULL, T = NULL, grid, params,
            distances, dim,...)

}
\arguments{
  \item{model,params}{\argModel }
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{distances,dim}{\argDistances}
 \item{...}{\argDots}
}
\details{
 \command{RFcovmatrix} returns a covariance matrix. Here, a matrix of coordinates (\code{x}) or a vector or a matrix of \code{distances}
 is expected.
 
 \command{RFcovmatrix} also allows for variogram models. Then the negative of the
 variogram matrix is returned.
}

\value{

 \command{RFcovmatrix} returns a covariance matrix. 
}

\me

\seealso{	
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RFfctn}},
 \command{\link{RFcalc}},
 \command{\link{RFcov}},
 \command{\link{RFpseudovariogram}},
 \command{\link{RFvariogram}}.
}

\examples{\dontshow{StartExample()}

##################################################
# Example: get covariance matrix C(x_i,x_j)    
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


\dontshow{FinalizeExample()}}
\keyword{spatial}






