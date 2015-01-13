\name{plot-method}
\docType{methods}
\alias{plot-method}
\alias{plot,RFgridDataFrame,missing-method}
\alias{plot,RFpointsDataFrame,missing-method}
\alias{plot,RFspatialDataFrame,missing-method}
\alias{plot,RFdataFrame,missing-method}
\alias{plot,RFdataFrame,matrix-method}
\alias{plot,RFdataFrame,data.frame-method}
\alias{plot,RFdataFrame,RFdataFrame-method}
\alias{persp,RFspatialGridDataFrame-method}
\alias{plot,RFspatialGridDataFrame,missing-method}
\alias{plot,RFspatialDataFrame,RFspatialPointsDataFrame-method}
\alias{plot,RFspatialDataFrame,matrix-method}
\alias{plot,RFspatialDataFrame,data.frame-method}
\alias{plot,RFspatialGridDataFrame,RFspatialPointsDataFrame-method}
\alias{plot,RFspatialDataFrame,RFspatialGridDataFrame-method}
\alias{plot,RFspatialPointsDataFrame,RFspatialGridDataFrame-method}
\alias{plot,RFspatialPointsDataFrame,missing-method}
\alias{plot,RFspatialPointsDataFrame,RFspatialPointsDataFrame-method}
\alias{plot,RFempVariog,missing-method}
\alias{plot,RFfit,missing-method}
\alias{plot,RMmodel,missing-method}
\alias{points.RMmodel}
\alias{lines.RMmodel}
\title{Methods for function \code{plot} in package \pkg{RandomFields}}
\description{
Plot methods are implemented for simulated random fields (objects of
class \code{\link[=RFsp-class]{RFsp}}), explicit covariance models (objects of class \code{\link[=RMmodel-class]{RMmodel}}),
empirical variograms (objects of class \code{\link[=RFempVariog-class]{RFempVariog}}) and fitted models (objects of class \code{\link[=RFfit]{RFfit}}).
}

\section{Methods}{
\describe{
\item{\code{signature(x = "RFspatialGridDataFrame", y =
 "missing")}}{Generates nice image plots of simulation results for
 simulation on a grid and space-time-dimension \eqn{\ge 2}{>1}. If
 space-time-dimension \eqn{\ge 3}{>2}, plots are on 2-dimensional
 subspaces. 
 Handles multivariate random fields (\code{.RFparams$vdim>1}) as well as
 repeated iid simulations (\code{.RFparams$vdim>n}).}
\item{\code{signature(x = "RFspatialGridDataFrame", y =
 "RFspatialPointsDataFrame")}}{
 Similar to method for \code{y="missing"}, but additionally adds the
 points of \code{y}. Requires \code{MARGIN.slices=NULL} and 
 \code{all.equal(x@.RFparams, y@.RFparams)}.}
\item{\code{signature(x = "RFspatialGridDataFrame", y =
 "matrix")}}{
 Similar to method for \code{y="missing"}, but additionally adds the
 points of \code{y}. Requires \code{MARGIN.slices=NULL} and 
 \code{all.equal(x@.RFparams, y@.RFparams)}.}
\item{\code{signature(x = "RFspatialPointsDataFrame", y =
 "RFspatialGridDataFrame")}}{
 Throws an error. Probably \code{x} and \code{y} have been interchanged.}
\item{\code{signature(x = "RFspatialPointsDataFrame", y =
 "missing")}}{Similar to method for class
 \code{\link[=RFspatialGridDataFrame-class]{RFspatialGridDataFrame}}, but for non-gridded simulation results.
 Instead of a grid, only existing points are plotted with colors indicating
 the value of the random field at the respective location. Handles
 multivariate random fields (\code{.RFparams$vdim>1}) as well as
 repeated iid simulations (\code{.RFparams$vdim>n}).}
\item{\code{signature(x = "RFspatialPointsDataFrame", y =
 "RFspatialPointsDataFrame")}}{
 Similar to method for \code{y="missing"}, but additionally adds the
 points of \code{y}. Requires 
 \code{all.equal(x@.RFparams, y@.RFparams)}.}

\item{\code{signature(x = "RFgridDataFrame", y = "missing")}}{Generates
 plots of simulation results for space-time-dimension \eqn{=1}.
 Handles different values for the number of repetitions as well as
 multivariate responses.}

\item{\code{signature(x = "RFpointsDataFrame", y = "missing")}}{Similar
 to method for class \code{\link[=RFgridDataFrame-class]{RFgridDataFrame}}, but for non-gridded data.}

\item{\code{signature(x = "RFempVario", y = "missing")}}{Gives nice
 plots of the empirical variogram; handles binning in up to three
 space-dimensions and a time-dimension, where the empirical variogram
 is plotted along lines which are directed according to the
 angle-centers given in \code{x@phi.centers} and
 \code{x@theta.centers}; arbitrary theoretical model curves
 can be added to the plot by using the argument \code{model}.
 If no bins are given, i.e. (\code{x@bin=NULL}),
 \command{\link[graphics]{image}}-plots are generated.}

\item{\code{signature(x = "RFfit", y = "missing")}}{Combines the plot of
 the empirical variogram with the estimated covariance or variogram
 model (theoretical) curves; further models can be added via the
 argument \code{model}.}

\item{\code{signature(x = "RMmodel", y = "missing")}}{Generates
 covariance function or variogram function plots in one or two
 dimensions.}

}



%% 2.1.3 Documenting S4 classes and methods
%% ----------------------------------------
 
%% There are special ways to use the ? operator, namely class?TOPIC
%% and methods?TOPIC, to access documentation for S4 classes and
%% methods, respectively. This mechanism depends on conventions for the
%% topic names used in alias entries. The topic names for S4 classes
%% and methods respectively are of the form
 
%% CLASS-class
%% GENERIC,SIGNATURE_LIST-method
 
%% where SIGNATURE_LIST contains the names of the classes in the signature
%% of the method (without quotes) separated by , (without whitespace),
%% with ANY used for arguments without an explicit specification. E.g.,
%% genericFunction-class is the topic name for documentation for the S4
%% class "genericFunction", and coerce,ANY,NULL-method is the topic
%% name for documentation for the S4 method for coerce for signature
%% c("ANY", "NULL").
 
%% Skeletons of documentation for S4 classes and methods can be
%% generated by using the functions promptClass() and promptMethods()
%% from package *methods*. If it is necessary or desired to provide an
%% explicit function declaration (in a usage section) for an S4 method
%% (e.g., if it has \sQuote{surprising arguments} to be mentioned explicitly),
%% one can use the special markup
 
%% % \S4method{GENERIC}{SIGNATURE_LIST}(ARGUMENT_LIST)
 
%% %(e.g., \S4method{coerce}{ANY,NULL}(from, to)).
 
%% To allow for making full use of the potential of the on-line
%% documentation system, all user-visible S4 classes and methods in a
%% package should at least have a suitable alias entry in one of the
%% packages Rd files. If a package has methods for a function defined
%% originally somewhere else, and does not change the underlying default
%% method for the function, the package is responsible for documenting the
%% methods it creates, but not for the function itself or the default
%% method.
 
%% See help("Documentation", package = "methods") for more
%% information on using and creating on-line documentation for S4 classes
%% and methods.
 


}

\usage{
\S4method{plot}{RFspatialDataFrame,missing}(x, y,
 MARGIN=c(1,2), MARGIN.slices=NULL,
 n.slices = if (is.null(MARGIN.slices)) 1 else 10,
 nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), select.variables, zlim, legend=TRUE,
 MARGIN.movie = NULL, ...)

\S4method{plot}{RFspatialDataFrame,RFspatialGridDataFrame}(x, y,
 MARGIN=c(1,2), MARGIN.slices=NULL,
 n.slices = if (is.null(MARGIN.slices)) 1 else 10,
 nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), select.variables, zlim, legend=TRUE,
  MARGIN.movie = NULL,...)

\S4method{plot}{RFspatialDataFrame,RFspatialPointsDataFrame}(x, y,
 MARGIN=c(1,2), MARGIN.slices=NULL,
 n.slices = if (is.null(MARGIN.slices)) 1 else 10,
 nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), select.variables, zlim, legend=TRUE,
  MARGIN.movie = NULL,...)

\S4method{plot}{RFspatialDataFrame,matrix}(x, y,
 MARGIN=c(1,2), MARGIN.slices=NULL,
 n.slices = if (is.null(MARGIN.slices)) 1 else 10,
 nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), select.variables, zlim, legend=TRUE,
  MARGIN.movie = NULL,...)
\S4method{plot}{RFspatialDataFrame,data.frame}(x, y,
 MARGIN=c(1,2), MARGIN.slices=NULL,
 n.slices = if (is.null(MARGIN.slices)) 1 else 10,
 nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), select.variables, zlim, legend=TRUE,
  MARGIN.movie = NULL,...)


\S4method{persp}{RFspatialGridDataFrame}(x, y,
 MARGIN=c(1,2), MARGIN.slices=NULL,
 n.slices = if (is.null(MARGIN.slices)) 1 else 10,
 nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), select.variables, zlim, legend=TRUE, 
 MARGIN.movie = NULL,...)


\S4method{plot}{RFdataFrame,missing}(x, y, nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), legend=TRUE, ...)

\S4method{plot}{RFdataFrame,RFdataFrame}(x, y, nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), legend=TRUE, ...)

\S4method{plot}{RFdataFrame,matrix}(x, y, nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), legend=TRUE, ...)

\S4method{plot}{RFdataFrame,data.frame}(x, y, nmax=6, 
 plot.variance = (!is.null(x@.RFparams$has.variance) &&
 x@.RFparams$has.variance), legend=TRUE, ...)

\S4method{plot}{RFempVariog,missing}(x, model=NULL,
%nmax.phi=6, nmax.theta=3, nmax.T=3,
nmax.phi=NA, nmax.theta=NA, nmax.T=NA,
 plot.nbin=TRUE, plot.sd=FALSE, variogram=TRUE, boundaries = TRUE,...)

\S4method{plot}{RFfit,missing}(x, model=NULL,  fit.method="ml",
 nmax.phi=NA, nmax.theta=NA, nmax.T=NA,
 plot.nbin=TRUE, plot.sd=FALSE, variogram = TRUE, boundaries = TRUE,...)

\S4method{plot}{RMmodel,missing}(x, y, dim=1, n.points=200,
 fct.type=NULL, MARGIN, fixed.MARGIN, maxchar=15, ...)

\S3method{points}{RMmodel}(x, y, n.points=200, fct.type=NULL, ...)

\S3method{lines}{RMmodel}(x, y, n.points=200, fct.type=NULL, ...)

%\S3method{contour}{RMmodel}(x, y, n.points=200, fct.type=NULL, ...)
}


\arguments{
  \item{x}{object of class \code{\link[=RFsp-class]{RFsp}} or
    \command{\link[=RFempVariog-class]{RFempVario}} or
    \command{\link[=RFfit-class]{RFfit}} or
    \command{\link[=RMmodel-class]{RMmodel}}; in the latter case, \code{x} can
    be any sophisticated model but it must be either stationary or a
    variogram model} 
  \item{y}{ignored in most methods}
  \item{MARGIN}{vector of two; two integer values giving the coordinate
    dimensions 
    w.r.t. which the field or the covariance model
    is to be plotted; in all other directions,
    the first index is taken}
  \item{MARGIN.slices}{integer value; if \eqn{[space-time-dimension>2]},
    \code{MARGIN.slices} can specify a third dimension w.r.t. which a
    sequence of slices is plotted. Currently only works for grids.} 
  \item{fixed.MARGIN}{only for \code{class(x)=="RMmodel"} and if
    \code{dim > 2}; a vector of length \code{dim}-2 with distance values
    for the coordinates that are not displayed} 
  \item{n.slices}{integer value. 
    The number of slices to be plotted; if
    \code{n.slices>1}, \code{nmax} is set to 1. Or \code{n.slices} is a
    vector of 3 elements: start, end, length. Currently only works for grids.}
  \item{nmax}{the maximal number of the \code{x@.RFparams$n} iid copies
    of the field that are to be plotted}
  \item{MARGIN.movie}{integer. If given a sequence of figures is shown for this
    direction. This option is in an experimental stage.
    It works only for grids.
  }
  \item{...}{arguments to be passed to methods; mainly graphical
    arguments, or further models in case of class \code{'RMmodel'},
    see Details.
  }
  \item{fit.method}{character; only for \code{class(x)=="RFfit"}; a
    vector of slot names for which the fitted covariance or variogram
    model is to be plotted; should be a subset of 
    \code{slotNames(x)} for which the corresponding slots are of class
    \code{"RMmodelFit"}; by default, the maximum likelihood fit will be
    plotted} 
  \item{nmax.phi}{integer; only for \code{class(x)=="RFempVario"};
    the maximal number of bins of angle phi that are to be
    plotted}
  \item{nmax.theta}{integer; only for \code{class(x)=="RFempVario"};
    the maximal number of bins of angle theta that are to be
    plotted} 
  \item{nmax.T}{integer; only for \code{class(x)=="RFempVario"};
    the maximal number of different time bins that are to be
    plotted} 
  \item{plot.nbin}{logical; only for \code{class(x)=="RFempVario"};
    indicates whether the number of pairs per bin are to be plotted}
  \item{plot.sd}{logical; only for \code{class(x)=="RFempVario"};
    indicates whether the calculated standard deviation (\code{x@sd}) is
    to be plotted (in form of arrows of length +-1*sd)}
  \item{variogram}{logical; This argument should currently not be set
    by the user. If \code{TRUE} then the empirical variogram
    is plotted, else an estimate for the covariance function} 
  \item{boundaries}{logical; only for \code{class(x)=="RFempVario"} and
    the anisotropic case where \code{model} is given.
    As the empirical variogram is calculated on a sector of angles,
    no exact variogram curve corresponds to the mean values in this
    sector.
    If \code{boundaries=TRUE} the values of the variogram on the sector
    boundaries are plotted. If \code{FALSE} some kind of mean
    model values are plotted. Neither the boundaries may contain
    the values of emprical variogram nor does the mean values
    need to be close the emprical variogram.  
  } 
  \item{dim}{must equal 1 or 2; only for \code{class(x)=="RMmodel"}; the
    covariance function and the variogram are plotted as a function of
    \eqn{\R^\code{dim}}{R^\code{dim}}.}
  \item{n.points}{integer; only for \code{class(x)=="RMmodel"}; the
    number of points at which the model 
    evaluated (in each dimension); defaults to 200} 
  \item{fct.type}{character; only for \code{class(x)=="RMmodel"}; must
    equal \code{NULL}, \code{"Cov"} or \code{"Variogram"}; controls
    whether the covariance (\code{fct.type="Cov"}) or the
    variogram (\code{fct.type="Variogram"}) is plotted; \code{NULL}
    implies automatic choice, where \code{"Cov"} is chosen whenever the model is
    stationary} 
  \item{model}{object of class
    \command{\link[=RMmodel-class]{RMmodel}}; only for
    \code{class(x)=="RFempVario"} or \code{class(x)=="RFfit"}; a list of
    covarianve or variogram models that are to be plotted into the same
    plot as the empirical variogram (and the fitted models)}
  \item{plot.variance}{logical, whether variances should be plotted if
    available}
  \item{select.variables}{vector of integers or list of vectors.
    The argument is only of interest for multivariate models.
    Here, \code{length(select.variables)} gives the number of pictures shown
    (excuding the plot for the variances, if applicable).
    If \code{select.variables} is a vector of integers then exactly these
    components are shown. If \code{select.variables} is a list, each element
    should be a vector up to length \eqn{l\le 3}{l <= 3}:
    \itemize{
      \item \eqn{l=1}\cr the component is shown in the usual way
      \item \eqn{l=2} \cr the two components are interpreted as vector and
      arrows are plotted
      \item \eqn{l=3} \cr the first component is shown as single component;
      the remaining two component are interpreted as a vector and
      plotted into the picture of the first component 
    }
  }
  \item{legend}{logical, whether a legend should be plotted}
  \item{zlim}{vector of length 2 with the usual meaning.
    In case of multivariate random fields, \code{zlim} can also be a character
    wih the value \sQuote{joint} indicating that all plotted compoments shall
    have the same zlim OR a matrix with two rows, where the i-th column
    gives the zlim of the i-th variable OR a list with entries named
    \code{data} and \code{var} if a separate \code{zlim}
    for the Kriging variance is to
    be used.}
  \item{maxchar}{maximum number of characters used in the legend
    (for multiplots only)
  }
}

\details{
  Internally, \code{...} are passed to \code{image} and
  \code{plot.default}, respectively; if, by default, multiple colors,
  xlabs or ylabs are used, also vectors of suitable length can be
  passed as \code{col}, \code{xlab} and \code{ylab}, respectively.

  One exception is the use of \code{...} in \command{plot} for class
  \code{'RMmodel'}.
  Here, further models might be passed. All models must have names
  starting with \code{model}. If \code{'.'} is following in the name,
  the part of the name after the dot is shown in the legend. Otherwise
  the name is ignored and a standardized name derived from the model
  definition is shown in the legend. Note that for the first argument
  a name cannot be specified.
}

\author{Alexander Malinowski,
  Martin Schlather, \email{schlather@math.uni-mannheim.de} \url{http://ms.math.uni-mannheim.de/de/publications/software}
}


\keyword{methods}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


## define the model:
model <- RMtrend(mean=0.5) + # mean
         RMstable(alpha=1, var=4, scale=10) + # see help("RMstable")
                                        ## for additional arguments
         RMnugget(var=1) # nugget


#############################################################
## Plot of covariance structure

plot(model)
plot(model, xlim=c(0, 30))
plot(model, xlim=c(0, 30), fct.type="Variogram")
plot(model, xlim=c(-10, 20), fct.type="Variogram", dim=2)

#############################################################
## Plot of simulation results

## define the locations:
from <- 0
step <- .1 
len <- if (interactive()) 50 else 3 # nicer, but also time consuming if len=100
x1D <- GridTopology(from, step, len)
x2D <- GridTopology(rep(from, 2), rep(step, 2), rep(len, 2))
x3D <- GridTopology(rep(from, 3), rep(step, 3), rep(len, 3))

## 1-dimensional
sim1D <- RFsimulate(model = model, x=x1D, n=6) 
plot(sim1D, nmax=4)

## 2-dimensional
sim2D <- RFsimulate(model = model, x=x2D, n=6) 
plot(sim2D, nmax=4)
plot(sim2D, nmax=4, col=terrain.colors(64),
main="My simulation", xlab="my_xlab")

## 3-dimensional
model <- RMmatern(nu=1.5, var=4, scale=2)
sim3D <- RFsimulate(model = model, x=x3D) 
plot(sim3D, MARGIN=c(2,3), MARGIN.slices=1, n.slices=4)

 
#############################################################
## empirical variogram plots

x <- seq(0, 10, 0.05)
bin <- seq(from=0, by=.2, to=3)

model <- RMexp()
X <- RFsimulate(x=cbind(x), model=model)
ev1 <- RFempiricalvariogram(data=X, bin=bin)
plot(ev1)

model <- RMexp(Aniso = cbind(c(10,0), c(0,1)))
X <- RFsimulate(x=cbind(x,x), model=model)
ev2 <- RFempiricalvariogram(data=X, bin=bin, phi=3)
plot(ev2, model=list(exp = model))


#############################################################
## plot Fitting results
x <- seq(0, 2, len=if (interactive()) 21 else 6)
 
model <- RMexp(Aniso = cbind(c(10,0), c(0,1)))
X <- RFsimulate(x=cbind(x,x), model=model)
fit <- RFfit(~RMexp(Aniso=diag(c(NA, NA))), data=X, fit.nphi = 2,
             modus="easygoing")
plot(fit) 


#############################################################
## plot Kriging results 
model <- RMwhittle(nu=1.2, scale=2)
n <- if (interactive()) 200 else 5
x <- runif(n, max=step*len/2) 
y <- runif(n, max=step*len/2) # 200 points in 2 dimensional space
sim <- RFsimulate(model = model, x=x, y=y)

interpolate <- RFinterpolate(model=model, x=x2D, data=sim)
plot(interpolate)
plot(interpolate, sim)


#############################################################
## plotting vector-valued results
model <- RMdivfree(RMgauss(), scale=4)
x <- y <- seq(-10,10, if (interactive()) 0.5 else 10)
simulated <- RFsimulate(model = model, x=x, y=y, n=1)
plot(simulated)
plot(simulated, select.variables=list(1, 1:3, 4))



#############################################################
## options for the zlim argument
model <- RMdelay(RMstable(alpha=1.9, scale=2), s=c(0, 4)) +
         RMdelay(RMstable(alpha=1.9, scale=2), s=c(4, 0))
simu <- RFsimulate(model, x, y)

plot(simu, zlim=list(data=cbind(c(-6,2), c(-2,1)), var=c(5,6)))
plot(simu, zlim=cbind(c(-6,2), c(-2,1)))
plot(simu, zlim=c(-6,2))
plot(simu, zlim="joint")

\dontshow{FinalizeExample()}

}

