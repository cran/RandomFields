\name{RFinterpolate}
\alias{RFinterpolate}
\alias{kriging}
\title{Interpolation methods}
\description{
  The function allows for different methods of interpolation.
  Currently only various kinds of kriging are installed.
}
\usage{
RFinterpolate(model, x, y = NULL, z = NULL, T = NULL, grid, data,
              distances, dim, err.model, method = "ml", ...)
}
\arguments{
 \item{model}{string; covariance model, see \command{\link{RMmodel}}, or
 type \code{\link{RFgetModelNames}()} to get all options.}
 
 \item{x}{\eqn{(n \times d)}{(n x d)} matrix or vector of x
 coordinates, or object of class
 \code{\link[sp:GridTopology-class]{GridTopology}} or
   \code{\link[raster]{raster}};
 coordinates of \eqn{n} points to be kriged.
 For more options see \link{RFsimulateAdvanced}.
 }
 \item{y}{optional vector of y coordinates}
 \item{z}{optional vector of z coordinates}
 \item{T}{optional vector of time coordinates, 
 \code{T} must always be an equidistant vector.
 Instead of \code{T=seq(from=From, by=By, len=Len)} one may also write
 \code{T=c(From, By, Len)}.
 }
 \item{grid}{logical; determines whether the vectors \code{x},
 \code{y}, and \code{z} should be
 interpreted as a grid definition; \code{RandomFields} can
 find itself the correct value in nearly all cases.
 See also \link{RFsimulateAdvanced}.
 }
 \item{data}{
 Matrix, data.frame or object of class \command{\link[=RFsp-class]{RFsp}};
 coordinates and response values of measurements;
 If a matrix or a data.frame, the first columns are
 interpreted as coordinate vectors, and the last column(s) as
 (multiple) measurement(s) of the field which are kriged 
 \emph{separately}; if the argument \code{x} is missing,
 \code{data} may contain \code{NA}s, which are then replaced by
 the kriged values (imputing); 
 for details on matching of variable names see 
 \link{RFsimulateAdvanced}; if of class \command{\link[=RFsp-class]{RFsp}}
 }
 \item{distances}{another alternative to pass the (relative)
 coordinates,
 see \link{RFsimulateAdvanced}.
 }
 \item{dim}{
 Only used if \code{distances} are given.
 }
 \item{err.model}{For conditional simulation and random imputing only.\cr
 Usually \code{err.model=RMnugget(var=var)}, or not given at all
 (error-free measurements).
 }
 \item{method}{
   character. A single method out of \code{methods} or
   \code{sub.methods},
   see \command{\link{RFfit}}.
 }
 \item{...}{
 Two intrinsic arguements can be given:
 \describe{
 \item{given}{matrix. If \code{data} is a matrix, then the coordinates
 can be given separately, namely by \code{given} where, in each row,
 a single location is given}
 \item{MARGIN}{character or integer.
	This argument is only given when the mean is kriged, i.e.
	\code{method='M'} (see \command{\link{RFoptions}}).
	Then it has a similar meaning as in \command{\link{apply}}:
	e.g. if \code{MARGIN="T"} in a space-time model, then spatial
	means are calculated for each temporal instance.
 }
 }
 For further optional arguments, see \command{\link{RFoptions}}.
 }
 % \item{trend}{only used for universal and intrinsic kriging. In case of
 % universal kriging \code{trend} is a non-negative integer (monomials
 % up to order k as trend functions), a list of functions or a formula (the
 % summands are the trend functions); you have the choice of using either
 % x, y, z or X1, X2, X3,... as spatial variables;
 % in case of intrinsic kriging trend should be a nonnegative integer which
 % is the order of the underlying model.
 % }
 % \item{pch}{Kriging procedures are quite time consuming in general. 
 % The character \code{pch} is printed after roughly
 % each 80th part of calculation.}
 % \item{return.variance}{logical. If \code{FALSE} the kriged field is
 % returned. If \code{TRUE} a list of two elements, \code{estim} and
 % \code{var}, i.e. the kriged field and the kriging variances,
 % is returned.}
 % \item{allowdistanceZero}{if \code{TRUE} then
 % identical locations are slightly scattered}
 % \item{cholesky}{
 % if \code{TRUE} cholesky decomposition is used instead of LU.
 %}
}

\note{Important options are
  \itemize{
    \item \code{method} (overwriting the automatically detected variant
    of kriging)
    \item \code{return_variance} (returning also the kriging variance)
    \item \code{locmaxm} (maximum number of conditional values before
    neighbourhood kriging is performed)
    \item \code{fillall} imputing estimates location by default
    \item \code{varnames} and \code{coordnames} in case
    \code{data.frame}s are used to tell which column contains the data
    and the coordinates, respectively.
}}

\details{
 In case of intrinsic cokriging (intrinsic kriging for a multivariate
 random fields) the pseudo-cross-variogram is used (cf. Ver Hoef and
 Cressie, 1991).

}
\value{
 The value depends on the additional argument \code{variance.return},
 see \command{\link{RFoptions}}. 
 
 If \code{variance.return=FALSE} (default), \code{Kriging} returns a
 vector or matrix of kriged values corresponding to the
 specification of \code{x}, \code{y}, \code{z}, and
 \code{grid}, and \code{data}.\cr
 
 \code{data}: a vector or matrix with \emph{one} column\cr
 * \code{grid=FALSE}. A vector of simulated values is
 returned (independent of the dimension of the random field)\cr
 * \code{grid=TRUE}. An array of the dimension of the
 random field is returned (according to the specification
 of \code{x}, \code{y}, and \code{z}).\cr
 
 \code{data}: a matrix with \emph{at least two} columns\cr
 * \code{grid=FALSE}. A matrix with the \code{ncol(data)} columns
 is returned.\cr
 * \code{grid=TRUE}. An array of dimension
 \eqn{d+1}{d+1}, where \eqn{d}{d} is the dimension of
 the random field, is returned (according to the specification
 of \code{x}, \code{y}, and \code{z}). The last
 dimension contains the realisations.

 If \code{variance.return=TRUE}, a list of two elements, \code{estim} and
 \code{var}, i.e. the kriged field and the kriging variances,
 is returned. The format of \code{estim} is the same as described
 above. The format of \code{var} is accordingly.
 }
\references{
 Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.

 Cressie, N.A.C. (1993)
 \emph{Statistics for Spatial Data.}
 New York: Wiley.
 
 Goovaerts, P. (1997) \emph{Geostatistics for Natural Resources
 Evaluation.} New York: Oxford University Press.
 
 Ver Hoef, J.M. and Cressie, N.A.C. (1993)
 Multivariate Spatial Prediction.
 \emph{Mathematical Geology} \bold{25}(2), 219-240.
 
 Wackernagel, H. (1998) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 2nd edition. 
}

\author{
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
 
 Marco Oesting, \email{oesting@math.uni-mannheim.de}

 \subsection{Author(s) of the code:}{
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}

 Alexander Malinowski, \email{malinowski@math.uni-mannheim.de}
 
 Marco Oesting, \email{oesting@math.uni-mannheim.de}
 }
}

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFempiricalvariogram}},
 \code{\link[=RandomFields-package]{RandomFields}},
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## Preparation of graphics
if (interactive()) dev.new(height=7, width=16)
RFoptions(always_close_screen=FALSE)

## creating random variables first
## here, a grid is chosen, but does not matter
p <- 3:8
points <- as.matrix(expand.grid(p,p))
model <- RMexp() + RMtrend(mean=1)
data <- RFsimulate(model, x=points)
plot(data)
x <- seq(0, 9, 0.25)

## Simple kriging with the exponential covariance model
model <- RMexp()
z <- RFinterpolate(model, x=x, y=x, data=data)
plot(z, data)

## Simple kriging with mean=4 and scaled covariance
model <- RMexp(scale=2) + RMtrend(mean=4)
z <- RFinterpolate(model, x=x, y=x, data=data)
plot(z, data)


## Ordinary kriging
model <- RMexp() + RMtrend(mean=NA)
z <- RFinterpolate(model, x=x, y=x, data=data)
plot(z, data)


\dontshow{\dontrun{

## alternatively 

## Intrinsic kriging
model <- RMfbm(a=1)
z <- RFinterpolate(krige.meth="U", model, x, x, data=data)
screen(scr <- scr+1); plot(z, data)


## Interpolation nicht korrekt
## Intrinsic kriging with Polynomial Trend
model <- RMfbm(a=1) + RMtrend(polydeg=2)
z <- RFinterpolate(model, x, x, data=data)
screen(scr <- scr+1); plot(z, data)
}}

\dontshow{\dontrun{
## Universal kriging with trend as formula
model <- RMexp() + RMtrend(arbit=function(X1,X2) sin(X1+X2)) +
 RMtrend(mean=1)
z <- RFinterpolate(model, x=x, y=x, data=data)
screen(scr <- scr+1); plot(z, data)

## Universal kriging with several arbitrary functions
model <- RMexp() + RMtrend(arbit=function(x,y) x^2 + y^2) +
 RMtrend(arbit=function(x,y) (x^2 + y^2)^0.5) + RMtrend(mean=1)
z <- RFinterpolate(model, x=x, y=x, data=data)
screen(scr <- scr+1); plot(z, data)
}}

% folgender Befehl muss unbedingt drin bleiben
close.screen(all = TRUE)

\dontshow{RFoptions(always_close_screen=TRUE);
          close.screen(all.screens=TRUE);
          while (length(dev.list()) >= 2) dev.off()}
\dontshow{FinalizeExample()}

}


\keyword{spatial}%-- one or more ...
