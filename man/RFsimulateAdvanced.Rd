\name{RFsimulateAdvanced}
\alias{RFsimulateAdvanced}
\title{Simulation of Random Fields -- Advanced}


\description{
 This function simulates \bold{unconditional} random fields:
  \itemize{
    \item univariate and multivariate,
    spatial and spatio-temporal Gaussian random fields
    \item \link[=RPpoisson]{stationary Poisson fields}
    \item \link[=RPchi2]{Chi2 fields}
    \item \link[=RPt]{t fields}
    \item \link[=RPbernoulli]{Binary fields}
    \item \link[=RPmaxstable]{stationary max-stable random fields}.
  }

  It also simulates \bold{conditional} random fields for
  \itemize{
    \item univariate and multivariate,
    spatial and spatio-temporal Gaussian random fields.
  }

  For basic simulation of Gaussian random fields, see \link{RFsimulate}.
  See  \link{RFsimulate.more.examples} and
  \link{RFsimulate.sophisticated.examples} for further examples.

}


\arguments{
 \item{model}{object of class \code{\link[=RMmodel-class]{RMmodel}},
 \command{\link{RFformula}} or \command{\link[stats]{formula}};
 specifies the model to be simulated

 \itemize{
 \item if of class \code{\link[=RMmodel-class]{RMmodel}}, \code{model}
 specifies
 
 \itemize{
   \item the type of random field by using \code{RP}functions,
   e.g.,
   \itemize{ 
     \item \command{\link{RPgauss}}: Gaussian random field (default
     if none of the functions in the list is given)
     %	 \item \command{\link{RPbernoulli}}:
     %	 indicator field, based on a Gaussian field
     %	 \item \command{\link{RPpoisson}}: shot noise model
     %	 \item \command{\link{RPbrownresnick}}: Brown-Resnick processes
     %	 \item \command{\link{RPopitz}}: Extremal t process
     \item \command{\link{RPsmith}}: Smith model
     %	 \item \command{\link{RPschlather}}: Extremal Gaussian process
     %	 (also called Schlather model)
     See\code{ \link{RP}} for an overview.
   }
   
   \item the covariance or variogram model in case of a Gaussian random
   field (\command{\link{RPgauss}}) and for fields based on Gaussian fields
   (e.g. \command{\link{RPbernoulli}});
   type \code{\link{RFgetModelNames}(type="variogram")}
   for a list of available models; see also \command{\link{RMmodel}}.
   
   \item the shape function in case of a shot noise process; type
   \code{\link{RFgetModelNames}(type='shape')} for a list of
   available models.
   
%   \item the main structural element of a max-stable field
%   (\command{\link{RPbrownresnick}}, \command{\link{RPopitz}},
%   \command{\link{RPsmith}} and \command{\link{RPschlather}})	
 }
 
 \item if of class \code{\link{RFformula}} or \command{\link[stats]{formula}},
 \code{submodel} specifies a linear mixed model where random
 effects can be modelled by Gaussian random fields;
 see \command{\link{RFformula}} for details on model
 specification.
 }
 }
 \item{x}{matrix of coordinates, or vector of x coordinates, or object
 of class \code{\link[sp:GridTopology-class]{GridTopology}} or
   \code{\link[raster]{raster}}; 
 if matrix, \code{ncol(x)} is the dimension
 of the index space; matrix notation is required in case of more than 3 space
 dimensions; in this case, if \code{grid=FALSE}, \code{x_ij} is the
 i-th coordinate in the j-th dimension; otherwise, if
 \code{grid=TRUE}, the columns of \code{x} are interpreted as
 gridtriples (see \code{grid}); if of class \code{\link[sp:GridTopology-class]{GridTopology}},
 \code{x} is interpreted as grid definition and \code{grid}
 is automatically set to \code{TRUE}. 
 }
 \item{y}{optional vector of y coordinates, ignored if \code{x}
 is a matrix}
 \item{z}{optional vector of z coordinates, ignored if \code{x}
 is a matrix}
 \item{T}{optional vector of time coordinates, 
 \code{T} must always be an equidistant vector or given in a
 gridtriple format (see argument \code{grid}); for each component
 of \code{T}, the random field is simulated at all location points.
 }
 \item{grid}{logical; determines whether the vectors \code{x},
 \code{y}, and \code{z} or the columns of \code{x} should be
 interpreted as a grid definition (see Details). If \code{grid=TRUE},
 either \code{x}, \code{y}, and \code{z} must
 be equidistant vectors in ascending order or the columns of \code{x}
 must be given in the gridtriple format:
 \code{c(from, stepsize, len)}.


 Note: If \code{grid} is not given, \command{\link{RFsimulate}} tries to
 guess what is meant.
 
 \code{c(from, stepsize, len)} (see Details)
 } 
 \item{data}{matrix, data.frame or object of class \command{\link[=RFsp-class]{RFsp}};
 coordinates and response values of
 measurements in case that conditional simulation is to
 be performed;
 if a matrix or a data.frame, the first columns are
 interpreted as coordinate vectors, and the last column(s) as
 (multiple) measurement(s) of the field; if \code{x} is missing,
 \code{data} may contain \code{NA}s, which are then replaced by
 conditionally simulated values; if \code{data} is missing, unconditional
 simulation is performed;
 for details on matching of variable names see Details; if of class
 \command{\link[=RFsp-class]{RFsp}}
 }
 \item{err.model}{same as \code{model}; gives the model of the
 measurement errors for the measured \code{data} (which must be given
 in this case!), see Details.
 \code{err.model=NULL} (default) corresponds to error-free
 measurements, the most common alternative is
 \code{err.model=\link{RMnugget}()};
 ignored if \code{data} is missing.
 }
 \item{distances}{object of class \code{\link{dist}} representing
 the upper triangular part of the matrix of Euclidean distances
 between the points at which the field is to be simulated; only
 applicable for stationary and isotropic models; if not \code{NULL},
 \code{dim} must be given and \code{x}, \code{y}, \code{z} and
 \code{T} must be missing or \code{NULL}.

 If \code{distances} are given, the current value of \code{spConform}, see
 \command{\link{RFoptions}}, is ignored and instead
 \code{spConform=FALSE} is used. (This fact may change in future.) % To Do
 
 }
 \item{dim}{
 integer; space or space-time dimension of the field
 }
 \item{n}{number of realizations to generate}
 \item{...}{further options and control arguments for the simulation
 that are passed to and processed by \command{\link{RFoptions}}
 }
}


\details{
 \command{\link{RFsimulate}} simulates different classes of random fields,
 controlled by the wrapping model.

 If the wrapping function of the \code{model} argument is a covariance
 or variogram model, i.e., one of the list obtained by
 \code{\link{RFgetModelNames}(type="variogram",
 group.by="type")}, by default, a Gaussian field 
 with the corresponding covariance structure is simulated.
 By default, the simulation method is chosen automatically through internal algorithms.
 The simulation method can be set explicitly by enclosing the
 covariance function with a \link[=Gaussian]{method specification}.

 If other than Gaussian fields are to be simulated, the \code{model}
 argument must be enclosed by a function specifying the type of the
 random field.


 There are different possibilities of passing the locations at which
 the field is to be simulated. If \code{grid=FALSE}, all coordinate
 vectors (except for the time component \eqn{T}) must have the same
 length and the field is only simulated at the locations given by the
 rows of \eqn{x} or of \code{cbind(x, y, z)}. If \eqn{T} is not
 missing, the field is simulated for all combinations
 \eqn{(x[i, ], T[k])} or \eqn{(x[i], y[i], z[i], T[k])},
 \eqn{i=1, ..., }\code{nrow(x)}, \eqn{k=1, ..., }\code{length(T)},
 even if \code{model} is not explicitly a space-time model. 
 \cr
 If \code{grid=TRUE}, the vectors \code{x},
 \code{y}, \code{z} and \code{T} or the columns of \code{x} and
 \code{T} are
 interpreted as a grid definition, i.e. the field is simulated at all
 locations \eqn{(x_i, y_j, z_k, T_l)}, as given by
 \code{expand.grid(x, y, z, T)}. 
 Here, \dQuote{grid} means \dQuote{equidistant in each direction}, i.e. all
 vectors must be equidistant and in ascending order.
 In case of more than 3 space dimensions, the coordinates must be
 given in matrix notation. To enable different grid lengths for each
 direction in combination with the matrix notation, the
 \dQuote{gridtriple} notation \code{c(from, stepsize, len)} is used:
 If \code{x}, \code{y},
 \code{z}, \code{T} or the columns of \code{x} are of length 3, they
 are internally replaced by \code{seq(from=from,
 to=from+(len-1)*stepsize, by=stepsize)} , i.e. the field
 is simulated at all locations\cr
 \code{expand.grid(seq(x$from, length.out=x$len, by=x$stepsize),
 seq(y$from, length.out=y$len, by=y$stepsize),
 seq(z$from, length.out=z$len, by=z$stepsize),
 seq(T$from, length.out=T$len, by=T$stepsize))}
 

 If \code{data} is passed, conditional simulation is performed. 
 \itemize{
 \item
 If of class \code{\link[=RFsp-class]{RFsp}},
 \code{ncol(data@coords)} must equal the dimension of the index
 space. If \code{data@data} contains only a single variable,
 variable names are optional. If \code{data@data} contains
 more than one variable, variables must be named and \code{model}
 must be given in the tilde notation \code{resp ~ ... } (see
 \command{\link{RFformula}}) and \code{"resp"} must be contained
 in \code{names(data@data)}.
 \item
 % Beschreibung hier stimmt nicht so ganz mit Examples unten ueberein
 If \code{data} is a matrix or a data.frame, either \code{ncol(data)}
 equals \eqn{(dimension of index space + 1)} and the order of the
 columns is (x, y, z, T, response) or, if \code{data} contains
 more than one 
 response variable (i.e. \code{ncol(data) > (dimension of index
 space + 1)}), \code{colnames(data)} must contain
 \code{colnames(x)} or those of \code{"x", "y", "z", "T"} that
 are not missing. The response variable name is matched with
 \code{model}, which must be given in the tilde notation. If
 \code{"x", "y", "z", "T"} are missing and \code{data} contains
 \code{NA}s, \code{colnames(data)} must contain an element which starts
 with \sQuote{data}; the corresponding column and those behind it are
 interpreted as the given data and those before the corresponding
 column are interpreted as the coordinates.
 \item
 If \code{x} is missing, \command{\link{RFsimulate}} searches for
 \code{NA}s in the data and performs a conditional simulation
 for them.
 }

 Specification of \code{err.model}:
 In geostatistics we have two different interpretations of a nugget
 effect: small scale variability and measurement error.
 The result of conditional simulation usually does not include the
 measurement error. Hence the measurement error \code{err.model}
 must be given separately. For sake of generality, any model (and not
 only the nugget effect) is allowed.
 Consequently, \code{err.model} is ignored
 when unconditional simulation is performed.
 
 }

\value{By default, 
 an object of the virtual class
 \command{\link[=RFsp-class]{RFsp}};
 result is of class
 \command{\link[=RFspatialGridDataFrame]{RFspatialGridDataFrame}}
 if \eqn{[space-time-dimension > 1]} and the coordinates are on a grid,
 result is of class 
 \command{\link[=RFgridDataFrame]{RFgridDataFrame}}
 if \eqn{[space-time-dimension = 1]} and the coordinates are on a grid,
 result is of class
 \command{\link[=RFspatialPointsDataFrame]{RFspatialPointsDataFrame}}
 if \eqn{[space-time-dimension > 1]} and the coordinates are not on a grid,
 result is of class 
 \command{\link[=RFpointsDataFrame]{RFpointsDataFrame}}
 if \eqn{[space-time-dimension = 1]} and the coordinates are not on a
 grid.

 The output format can be switched to the "old" array format using
 \code{\link{RFoptions}}, either by globally setting
 \code{\link{RFoptions}(spConform=FALSE)} or by passing \code{spConform=FALSE}
 in the call of \command{\link{RFsimulate}}.
 Then the object returned by \command{\link{RFsimulate}}
 depends on the arguments \code{n} and \code{grid} in the following way:\cr

 If \code{vdim > 1} the \code{vdim}-variate vector makes the first dimension.

 If \code{grid=TRUE} an array of the dimension of the
 random field makes the next dimensions. Here, the dimensions
 are ordered in the sequence \code{x}, \code{y}, \code{z}, \code{T}
 (if given).
 
 Else if no time component is given, then the values are passed as a
 single vector. Else if the time component is given the next 2
 dimensions give the space and the time, respectively.
 
 If \code{n > 1} the repetitions make the last dimension.

 Note: Conversion between the \pkg{sp}
 format and the conventional format can be
 done using the method \command{RFspDataFrame2conventional} and the
 function \command{conventional2RFspDataFrame}. 
 
 \code{InitRFsimulate} returns 0 if no error has occurred and a positive value
 if failed.\cr

}

\note{Advanced options are
  \itemize{
    \item \code{spConform} (suppressed return of S4 objects)
    \item \code{practicalrange} (forces range of covariances to be one)
    \item \code{exactness} (chooses the simulation method by precision)
    \item \code{seed} (sets \code{\link[base]{.Random.seed}} locally or
    globally)
  }
  See \command{\link{RFoptions}} for further options.
}

\references{
% Gneiting, T. and Schlather, M. (2004)
% Statistical modeling with covariance functions.
% \emph{In preparation.}
General
\itemize{
  \item
 Lantuejoul, Ch. (2002) \emph{Geostatistical simulation.}
 \bold{New York:} Springer.

 \item
 Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University.
 }

 Original work:
 \itemize{
 \item Circulant embedding:
 
 Chan, G. and Wood, A.T.A. (1997)
 An algorithm for simulating stationary Gaussian random fields.
 \emph{J. R. Stat. Soc., Ser. C} \bold{46}, 171-181.
 
 Dietrich, C.R. and Newsam, G.N. (1993)
 A fast and exact method for multidimensional Gaussian
 stochastic simulations.
 \emph{Water Resour. Res.} \bold{29}, 2861-2869. 
 
 Dietrich, C.R. and Newsam, G.N. (1996)
 A fast and exact method for multidimensional {G}aussian stochastic
 simulations: Extensions to realizations conditioned on direct and
 indirect measurement 
 \emph{Water Resour. Res.} \bold{32}, 1643-1652. 

 Wood, A.T.A. and Chan, G. (1994)
 Simulation of stationary Gaussian processes in \eqn{[0,1]^d}{[0,1]^d}
 \emph{J. Comput. Graph. Stat.} \bold{3}, 409-432.

 The code used in \cite{RandomFields} is based on
 Dietrich and Newsam (1996).
 
 \item Intrinsic embedding and Cutoff embedding:
 
 Stein, M.L. (2002)
 Fast and exact simulation of fractional Brownian surfaces.
 \emph{J. Comput. Graph. Statist.} \bold{11}, 587--599.
 
 
 Gneiting, T., Sevcikova, H., Percival, D.B., Schlather, M. and
 Jiang, Y. (2005)
 Fast and Exact Simulation of Large Gaussian Lattice Systems in
 \eqn{R^2}: Exploring the Limits
 \emph{J. Comput. Graph. Statist.} Submitted.

 \item Markov Gaussian Random Field:
 
 Rue, H. (2001) Fast sampling of Gaussian Markov random fields.
 \emph{J. R. Statist. Soc., Ser. B}, \bold{63} (2), 325-338.
 
 Rue, H., Held, L. (2005) \emph{Gaussian Markov Random Fields:
 Theory and Applications.}
 Monographs on Statistics and Applied Probability, no \bold{104},
 Chapman \& Hall.


 \item Turning bands method (TBM), turning layers:
 
 Dietrich, C.R. (1995) A simple and efficient space domain implementation
 of the turning bands method. \emph{Water Resour. Res.} \bold{31},
 147-156.
 
 Mantoglou, A. and Wilson, J.L. (1982) The turning bands method for
 simulation of random fields using line generation by a spectral
 method. \emph{Water. Resour. Res.} \bold{18}, 1379-1394.

 Matheron, G. (1973)
 The intrinsic random functions and their applications.
 \emph{Adv. Appl. Probab.} \bold{5}, 439-468.

 Schlather, M. (2004)
 Turning layers: A space-time extension of turning bands.
 \emph{Submitted}

 \item Random coins:
 
 Matheron, G. (1967) \emph{Elements pour une Theorie des Milieux
 Poreux}. Paris: Masson.
 }
}

\me

\seealso{
 \command{\link{RFoptions}},
 \command{\link{RMmodel}},
 \command{\link{RFgui}},
 \command{\link[=Gaussian]{methods for simulating Gaussian random fields}},
 \command{\link{RFfit}},
 \command{\link{RFvariogram}},
 \link{RFsimulate.more.examples},
 \link{RFsimulate.sophisticated.examples},
 \command{\link{RPgauss}}.
}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

\dontshow{\dontrun{

#############################################################
##                                                         ##
## Example 1: Specification of simulation method           ##
##                                                         ##
#############################################################

 ## usage of a specific method
 ## -- the complete list is obtained by RFgetMethodNames()
 model <- RMstable(alpha=1.5)
 x <- runif(100, max=20) 
 y <- runif(100, max=20) # 100 points in 2 dimensional space
 simulated <- RFsimulate(model = RPdirect(model), x=x, y=y) # cholesky
 plot(simulated)



#############################################################
##                                                         ##
## Example 2: Turnung band with different number of lines  ##
##                                                         ##
#############################################################
model <- RMstable(alpha=1.5)
x <- seq(0, 10, 0.01)
z <- RFsimulate(model = RPtbm(model), x=x, y=x)
plot(z)



#############################################################
##                                                         ##
## Example 3: Shot noise fields (random coins)             ##
##                                                         ##
#############################################################

 x <- GridTopology(0, .1, 500)

 z <- RFsimulate(model=RPpoisson(RMgauss()), x=x, mpp.intensity = 100)

 plot(z)
 par(mfcol=c(2,1))
 plot(z@data[,1:min(length(z@data), 1000)], type="l")
 hist(z@data[,1])
 
 
 z <- RFsimulate(x=x, model=RPpoisson(RMball()), mpp.intensity = 0.1)
 
 plot(z)
 par(mfcol=c(2,1))
 plot(z@data[,1:min(length(z@data), 1000)], type="l")
 hist(z@data[,1])
 



 #############################################################
 ##                                                         ##
 ## Example 4: a 2d random field based on                   ##
 ## covariance functions valid in 1d only                   ##
 ##                                                         ##
 #############################################################

x <- seq(0, 2, 0.1)
model <- RMfbm(alpha=0.5, Aniso=matrix(nrow=1, c(1, 0))) + 
         RMfbm(alpha=0.9, Aniso=matrix(nrow=1, c(0, 1)))
z <- RFsimulate(x, x, model=model)
plot(z)



#############################################################
##                                                         ##
## Example 5 : Brownian sheet                              ##
## (using Stein's method)                                  ##
##                                                         ##
#############################################################
 
# 2d
step <- 0.3 ## nicer, but also time consuming if step = 0.1
x <- seq(0, 5, step)
alpha <- 1 # in [0,2)
z <- RFsimulate(x=x, y=x, model=RMfbm(alpha=alpha))
plot(z)


# 3d
z <- RFsimulate(x=x, y=x, z=x,
 model=RMfbm(alpha=alpha))


#############################################################
##                                                         ##
## Example 5 : Non-Geometric anisotropy                    ##
##                                                         ##
#############################################################

x <- seq(0.1, 6, 0.12)
Aniso <- R.c(R.p(1)^2, R.p(2)^1.5)
z <- RFsimulate(RMexp(Aniso = Aniso) + 10, x, x)
plot(z) 
}}

\dontshow{FinalizeExample()}}




\keyword{spatial}


