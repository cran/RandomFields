\name{plot-method}
\docType{methods}
\alias{plot-method}
\alias{plot,RFgridDataFrame,missing-method}
\alias{plot,RFpointsDataFrame,missing-method}
\alias{plot,RFspatialGridDataFrame,missing-method}
\alias{plot,RFspatialPointsDataFrame,missing-method}

\alias{plot,RFgridDataFrame,matrix-method}
\alias{plot,RFpointsDataFrame,matrix-method}
\alias{plot,RFspatialGridDataFrame,matrix-method}
\alias{plot,RFspatialPointsDataFrame,matrix-method}

\alias{plot,RFgridDataFrame,data.frame-method}
\alias{plot,RFpointsDataFrame,data.frame-method}
\alias{plot,RFspatialGridDataFrame,data.frame-method}
\alias{plot,RFspatialPointsDataFrame,data.frame-method}

\alias{plot,RFgridDataFrame,RFgridDataFrame-method}
\alias{plot,RFgridDataFrame,RFpointsDataFrame-method}
\alias{plot,RFpointsDataFrame,RFgridDataFrame-method}
\alias{plot,RFpointsDataFrame,RFpointsDataFrame-method}
\alias{plot,RFspatialGridDataFrame,RFspatialGridDataFrame-method}
\alias{plot,RFspatialGridDataFrame,RFspatialPointsDataFrame-method}
\alias{plot,RFspatialPointsDataFrame,RFspatialGridDataFrame-method}
\alias{plot,RFspatialPointsDataFrame,RFspatialPointsDataFrame-method}

\alias{persp,RFspatialGridDataFrame-method}
\alias{contour.RFspatialGridDataFrame}

\alias{RFplotSimulation}
\alias{RFplotSimulation1D}

\title{Methods for function \code{plot} in package \pkg{RandomFields}}

\description{
  Plot methods are implemented for simulated random fields (objects of
  class \code{\link[=RFsp-class]{RFsp}}), explicit covariance models
  (objects of class \code{\link[=RMmodel-class]{RMmodel}}),  
  empirical variograms (objects of class
  \code{\link[=RFempVariog-class]{RFempVariog}}) and fitted models
  (objects of class \code{\link[=RFfit-class]{RFfit}}).
  
  The plot methods not described here are described together with the
  class itself, for instance,
  %\code{\link[=RFempVariog-class]{RFempVariog}})
  \code{\link[=RFfit-class]{RFfit}},
  \code{\link[=RFempVariog-class]{RFempVariog}}
  \code{\link[=RMmodel-class]{RMmodel}}.
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
RFplotSimulation(x, y, MARGIN=c(1,2), MARGIN.slices=NULL,
 n.slices = if (is.null(MARGIN.slices)) 1 else 10, nmax=6, 
 plot.variance = !is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance,
 select.variables, zlim, legend=TRUE,
 MARGIN.movie = NULL, file=NULL, speed = 0.3,
 height.pixel=300, width.pixel=height.pixel,
 ..., plotmethod="image")

RFplotSimulation1D(x, y, nmax=6,
  plot.variance=!is.null(x@.RFparams$has.variance) && x@.RFparams$has.variance,
  legend=TRUE, ...)

\S4method{plot}{RFgridDataFrame,missing}(x, y, ...)
\S4method{plot}{RFpointsDataFrame,missing}(x, y, ...)
\S4method{plot}{RFspatialGridDataFrame,missing}(x, y, ...)
\S4method{plot}{RFspatialPointsDataFrame,missing}(x, y, ...)

\S4method{plot}{RFgridDataFrame,matrix}(x, y, ...)
\S4method{plot}{RFpointsDataFrame,matrix}(x, y, ...)
\S4method{plot}{RFspatialGridDataFrame,matrix}(x, y, ...)
\S4method{plot}{RFspatialPointsDataFrame,matrix}(x, y, ...)

\S4method{plot}{RFgridDataFrame,data.frame}(x, y, ...)
\S4method{plot}{RFpointsDataFrame,data.frame}(x, y, ...)
\S4method{plot}{RFspatialGridDataFrame,data.frame}(x, y, ...)
\S4method{plot}{RFspatialPointsDataFrame,data.frame}(x, y, ...)

\S4method{plot}{RFgridDataFrame,RFgridDataFrame}(x, y, ...)
\S4method{plot}{RFgridDataFrame,RFpointsDataFrame}(x, y, ...)
\S4method{plot}{RFpointsDataFrame,RFgridDataFrame}(x, y, ...)
\S4method{plot}{RFpointsDataFrame,RFpointsDataFrame}(x, y, ...)
\S4method{plot}{RFspatialGridDataFrame,RFspatialGridDataFrame}(x, y, ...)
\S4method{plot}{RFspatialGridDataFrame,RFspatialPointsDataFrame}(x, y, ...)
\S4method{plot}{RFspatialPointsDataFrame,RFspatialGridDataFrame}(x, y, ...)
\S4method{plot}{RFspatialPointsDataFrame,RFspatialPointsDataFrame}(x, y, ...)

\S4method{persp}{RFspatialGridDataFrame}(x, ..., zlab="")
\S3method{contour}{RFspatialGridDataFrame}(x, ...)

}


\arguments{
  \item{x}{object of class \code{\link[=RFsp-class]{RFsp}} or
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
 % \item{fixed.MARGIN}{only for \code{class(x)=="RMmodel"} and if
 %   \code{dim > 2}; a vector of length \code{dim}-2 with distance values
 %   for the coordinates that are not displayed} 
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
  \item{file, speed, height.pixel, width.pixel}{
    In case \code{MARGIN.movie} and \code{file} is given an 'avi' movie
    is stored using the \command{mencoder} command with speed argument
    \code{speed}. As temporary files \code{file__###.png}  
    of size \code{width.pixel} x \code{height.pixel} are created.
  }
  \item{...}{arguments to be passed to methods; mainly graphical
    arguments, or further models in case of class \code{'RMmodel'},
    see Details.
  }
%  \item{f i tmethod}{character; only for \code{class(x)=="RFfit"}; a
%    vector of slot names for which the fitted covariance or variogram
%    model is to be plotted; should be a subset of 
%    \code{slotNames(x)} for which the corresponding slots are of class
%    \code{"RMmodelFit"}; by default, the maximum likelihood fit will be
%    plotted} 
 %  \item{dim}{must equal 1 or 2; only for \code{class(x)=="RMmodel"}; the
%    covariance function and the variogram are plotted as a function of
%    \eqn{\R^\code{dim}}{R^\code{dim}}.}
%  \item{n.points}{integer; only for \code{class(x)=="RMmodel"}; the
 %   number of points at which the model 
 %   evaluated (in each dimension); defaults to 200} 
%  \item{fct.type}{character; only for \code{class(x)=="RMmodel"}; must
%    equal \code{NULL}, \code{"Cov"} or \code{"Variogram"}; controls
%    whether the covariance (\code{fct.type="Cov"}) or the
%    variogram (\code{fct.type="Variogram"}) is plotted; \code{NULL}
%  implies automatic choice, where \code{"Cov"} is chosen whenever the model is
%    stationary} 
%  \item{model}{object of class
%    \command{\link[=RMmodel-class]{RMmodel}}; only for
%    \code{class(x)=="RFempVario"} or \code{class(x)=="RFfit"}; a list of
%    covarianve or variogram models that are to be plotted into the same
%    plot as the empirical variogram (and the fitted models)}
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
  \item{plotmethod}{string or function. Internal.
  }
  \item{zlab}{character. See \command{\link[graphics]{persp}}}
%  \item{maxchar}{maximum number of characters used in the legend
%    (for multiplots only)  }
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


\seealso{
   \command{\link{RFpar}}.
}

\keyword{methods}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}


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
image(model, xlim=c(-10, 20), fct.type="Variogram")
persp(model, xlim=c(-10, 20), fct.type="Variogram")

#############################################################
## Plot of simulation results

## define the locations:
from <- 0
step <- .1 
len <- 50  # nicer if len=100 %ok
\dontshow{if(RFoptions()$internal$examples_red){warning("reduced 'len'"); len<-3}} 
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
x <- seq(0, 1, len=21)
 
model <- RMexp(Aniso = cbind(c(10,0), c(0,1)))
X <- RFsimulate(x=cbind(x,x), model=model)
fit <- RFfit(~RMexp(Aniso=diag(c(NA, NA))), data=X, fit.nphi = 2,
             modus="easygoing")
plot(fit) 


#############################################################
## plot Kriging results 
model <- RMwhittle(nu=1.2, scale=2)
n <- 200
x <- runif(n, max=step*len/2) 
y <- runif(n, max=step*len/2) # 200 points in 2 dimensional space
sim <- RFsimulate(model = model, x=x, y=y)

interpolate <- RFinterpolate(model=model, x=x2D, data=sim)
plot(interpolate)
plot(interpolate, sim)


#############################################################
## plotting vector-valued results
model <- RMdivfree(RMgauss(), scale=4)
x <- y <- seq(-10,10, 0.5)
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

