\name{RFfctn}
\alias{RFfctn}
\title{Evaluate Covariance and Variogram Functions}
\description{
 \command{\link{RFfctn}} returns the values of a shape function. 
}
\usage{
RFfctn(model, x, y = NULL, z = NULL, T = NULL, grid, params, distances, dim,...)

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
 \command{\link{RFcovmatrix}} also allows for variogram models. Then the negative of the
 variogram matrix is returned.
}

\value{
\command{\link{RFfctn}} returns a vector. 
}

\me

\seealso{	
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}},
 \command{\link{RFcalc}},
 \command{\link{RFcov}},
 \command{\link{RFcovmatrix}},
 \command{\link{RFpseudovariogram}},
 \command{\link{RFvariogram}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMexp() - 1
RFfctn(model, 1:10)


\dontshow{FinalizeExample()}}
\keyword{spatial}






