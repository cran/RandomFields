\name{GaussRF}
\alias{GaussRF}
\alias{InitGaussRF}
\title{Gaussian Random Fields}
\description{
  These functions simulate stationary spatial and spatio-temporal
  Gaussian random fields using turning bands/layers, circulant embedding,
  direct methods, and the random coin method.
}
\usage{
GaussRF(x, y=NULL, z=NULL, T=NULL, grid, model, param, trend,
        method=NULL, n=1, register=0, gridtriple=FALSE)

InitGaussRF(x, y=NULL, z=NULL, T=NULL, grid, model, param, trend,
            method=NULL, register=0, gridtriple=FALSE)
}
\arguments{
  \item{x}{matrix of coordinates, or vector of x coordinates}
  \item{y}{vector of y coordinates}
  \item{z}{vector of z coordinates}
  \item{T}{vector of time coordinates, may only be given if
    the random field is defined as an anisotropic random field,
    i.e. if \code{model=list(list(model=,var=,k=,aniso=),...)}.
    \code{T} must always be given in the \code{gridtriple} format,
    independently how the spatial part is defined.
  }
  \item{grid}{logical; determines whether the vectors \code{x},
    \code{y}, and \code{z} should be
    interpreted as a grid definition, see Details.  \code{grid}
    does not apply for \code{T}.}
  \item{model}{string or list; covariance or variogram model,
    see \code{\link{CovarianceFct}}, or
    type \code{\link{PrintModelList}()} to get the list of all implemented
    models; see Details.}
  \item{param}{vector or matrix of parameters or missing, see Details
    and \code{\link{CovarianceFct}}; 
    The simplest form is that \code{param} is vector of the form
    \code{param=c(NA,variance,nugget,scale,...)}, in this order;\cr
    The dots \code{...} stand for additional parameters of the
    model.}
  \item{trend}{Not programmed yet.
    trend surface: number (mean) or a vector of length
    \eqn{d+1} (linear trend \eqn{a_0+a_1 x_1 + \ldots + a_d x_d}{
      a_0 +a_1 x_1 + ... + a_d x_d}), or function(x)}
  \item{method}{\code{NULL} or string; method used for simulating,
    see \code{\link{RFMethods}}, or
    type \code{\link{PrintMethodList}()} to get all options.
    If \code{model} is given as list then \code{method} may not be
    set if \code{model[[i]]$method}, \eqn{i=1,3,..} is given, and vice
    versa.}
  \item{n}{number of realisations to generate}
  \item{register}{0:9; place where intermediate calculations are stored;
    the numbers are aliases for 10 internal registers}
  \item{gridtriple}{logical. Only relevant if \code{grid==TRUE}.
    If \code{gridtriple==TRUE}
    then \code{x}, \code{y}, and \code{z} are of the
    form \code{c(start,end,step)}; if
    \code{gridtriple==FALSE} then \code{x}, \code{y}, and \code{z}
    must be vectors of ascending values
  }
}
\details{
  \code{GaussRF} can use different methods for the simulation,
  i.e., circulant embedding, turning bands, direct methods, and random
  coin method. 
  If \code{method==NULL} then \code{GaussRF} searches for a
  valid method.  \code{GaussRF} may not find the fastest method neither the
  most precise one.  It just finds any method among the available
  methods. (However it guesses what is a good choice.) 
  Note that some of the methods do not work for all covariance 
  or variogram models.

  \itemize{
    \item
    An isotropic random field is created by 
    \code{GaussRF} where
    \code{model} is the covariance or variogram model and the
    parameter is \code{param=c(mean,variance,nugget,scale, ...)}.
    Alternatively the \code{trend} can be given (not programmed yet); then
    \code{param=c(variance,nugget,scale, ...)}.

    \item
    Nested models can be defined in the same way as a nested
    \code{\link{CovarianceFct}}.  If the \code{trend} is not given
    it is set to 0.

    \item
    An anisotropic random field (i.e. zonal anisotropy, geometrical
    anisotropy, separable models, non-separable space-time models)
    and a random field based on multiplicative or nested models
    is defined as in the case of
    an anisotropic \code{\link{CovarianceFct}}. 
    If the \code{trend} is not given it is set to 0.
    The \code{method} may be specified by the global \code{method}
    or for each model separately, as additional parameter
    \code{method} for each entry of the list;
    note that methods can not be mixed within a multiplicative part.

    If \code{model=list(list(model=,var=,k=,aniso=),...)} then a time
    component might be given.  In case of \code{model="nugget"},
    \code{aniso} must still be given as a matrix.  Namely if
    \code{aniso} is a singular matrix then a zonal nugget effect
    is obtained.
   }  
 
  \code{GaussRF} calls initially \code{InitGaussRF},
  which does some basic checks on the validity of the parameters.  Then,
  \code{InitGaussRF} performs some first calculations, like the first
  Fourier transform in the circulant embedding method or the matrix
  decomposition for the direct methods.  Random numbers are not involved. 
  \code{GaussRF} then calls \code{\link{DoSimulateRF}} which uses the
  intermediate results and random numbers to create a simulation.

  When \code{InitGaussRF} checks the validity of the parameters, it
  also checks whether the previous simulation has had the same
  specification of the random field.  If so (and if
  \code{\link{RFparameters}()$STORING==TRUE}), the stored intermediate
  results are used instead of being recalculated. 
  
  Comments on specific parameters:
  \itemize{
  \item \code{grid==FALSE} : the vectors \code{x}, \code{y},
  and \code{z} are interpreted as vectors of coordinates
  \item \code{(grid==TRUE) && (gridtriple==FALSE)} : the vectors
  \code{x}, \code{y}, and \code{z}
  are increasing sequences with identical lags for each sequence. 
  A corresponding
  grid is created (as given by \code{expand.grid}). 
  \item \code{(grid==TRUE) && (gridtriple==FALSE)} : the vectors
  \code{x}, \code{y}, and \code{z}
  are triples of the form (start,end,step) defining a grid
  (as given by \code{expand.grid(seq(x$start,x$end,x$step),
                          seq(y$start,y$end,y$step),
                          seq(z$start,z$end,z$step))})
			
  %\item \code{model="nugget"} is one possibility to create independent
  %Gaussian random variables.  Without loss of efficiency, any
  %covariance function with parameter vector
  %\code{c(mean, 0, nugget, scale, ...)}
  %can also be used.  If \code{model="nugget"} is used
  %the second component of \code{param} must be zero.
  %% this has to be changed in later versions; 

  %\item The sum of the components variance and nugget in the argument
  %\code{param} equals the sill of the variogram. 

  \item \code{register} is a parameter which may never be
  used by most of the users (please let me know if you use it!).  In
  other words,
  the package will work fine if you ignore this parameter. 
  The parameter \code{register} is of interest in the following
  situation.  Assume you wish to create sequentially
  several realisations of two random fields \eqn{Z_1}{Z1} and
  \eqn{Z_2}{Z2} that have different
  specifications of the covariance/variogram models, i.e.
  \eqn{Z_1}{Z1}, \eqn{Z_2}{Z2}, \eqn{Z_1}{Z1}, \eqn{Z_2}{Z2},...
  Then, without using different registers, the algorithm
  will not be able to profit from already calculated intermediate
  results, as the specifications of the covariance/variogram model
  change every time. 
  However, using different registers allows for profiting from
  up to 10 stored intermediate results. 

  \item The strings for \code{model} and \code{method} may
  be abbreviated as long as the abbreviations match only one
  option.  See also \code{\link{PrintModelList}()} and
  \code{\link{PrintMethodList}()}
  
  \item Further control parameters for the simulation are set by means of
  \code{\link{RFparameters}(...)}.
  
  }
}
\note{
  The algorithms for all the simulation methods are controlled by
  additional parameters, see \code{\link{RFparameters}()}.  These
  parameters have an influence on the speed of the algorithm
  and the precision of the result. 
  The default parameters are chosen such that
  the simulations are fine for many models and their parameters. 
  If in doubt modify the example in \code{\link{EmpiricalVariogram}()}
  to check the precision.
}
\value{

  \code{InitGaussRF} returns 0 if no error has occurred and a positive value
  if failed.\cr

  \code{GaussRF} and \code{\link{DoSimulateRF}} return \code{NULL}
  if an error has occurred; otherwise the returned object
  depends on the parameters \code{n} and \code{grid}:\cr
    \code{n==1}:\cr
    * \code{grid==FALSE}.  A vector of simulated values is
    returned (independent of the dimension of the random field)\cr
    * \code{grid==TRUE}.  An array of the dimension of the
    random field is returned.\cr
    
    \code{n>1}:\cr
    * \code{grid==FALSE}.  A matrix is returned.  The columns
    contain the repetitions.\cr
    * \code{grid==TRUE}.  An array of dimension
    \eqn{d+1}{d+1}, where \eqn{d}{d} is the dimension of
    the random field, is returned.  The last
    dimension contains the repetitions.
}
\references{
  See \link{RFMethods} for the references.
}
\author{Martin Schlather, \email{martin.schlather@cu.lu}
  \url{http://www.cu.lu/~schlathe}}
\seealso{
  \code{\link{CovarianceFct}},
  \code{\link{DeleteRegister}},
  \code{\link{DoSimulateRF}},
  \code{\link{GetPracticalRange}},
  \code{\link{EmpiricalVariogram}},
  \code{\link{mleRF}},
  \code{\link{MaxStableRF}},
  \code{\link{RFMethods}},
  \code{\link{RandomFields}},
  \code{\link{RFparameters}},
  \code{\link{ShowModels}},
%  \code{\link{winddata}}.
}

\examples{
 #############################################################
 ##                                                         ##
 ## Examples using the symmetric stable model, also called  ##
 ## "powered exponential model"                             ## 
 ##                                                         ##
 #############################################################
 PrintModelList()    ## the complete list of implemented models
 model <- "stable"   
 mean <- 0
 variance <- 4
 nugget <- 1
 scale <- 10
 alpha <- 1   ## see help("CovarianceFct") for additional
              ## parameters of the covariance functions
 x <- seq(0, 20, 0.1) 
 y <- seq(0, 20, 0.1)     
 f <- GaussRF(x=x, y=y, model=model, grid=TRUE,
              param=c(mean, variance, nugget, scale, alpha))
 image(x, y, f)


 #############################################################
 ## ... using gridtriple
 x <- c(0, 20, 0.1)  ## note: vectors of three values, not a 
 y <- c(0, 20, 0.1)  ##       sequence
 f <- GaussRF(grid=TRUE, gridtriple=TRUE,
               x=x ,y=y, model=model,  
               param=c(mean, variance, nugget, scale, alpha))
 image(seq(x[1],x[2],x[3]), seq(y[1],y[2],y[3]), f)


 #############################################################
 ## arbitrary points
 x <- runif(100, max=20) 
 y <- runif(100, max=20)
 z <- runif(100, max=20) # 100 points in 3 dimensional space
(f <- GaussRF(grid=FALSE,
              x=x, y=y, z=z, model=model, 
              param=c(mean, variance, nugget, scale, alpha)))


 #############################################################
 ## usage of a specific method
 ## -- the complete list can be obtained by PrintMethodList()
 x <- runif(100, max=20) # arbitrary points
 y <- runif(100, max=20)
 (f <- GaussRF(method="dir",  # direct matrix decomposition
              x=x, y=y, model=model, grid=FALSE, 
              param=c(mean, variance, nugget, scale, alpha)))


 #############################################################
 ## simulating several random fields at once
 x <- seq(0, 20, 0.1)  # grid
 y <- seq(0, 20, 0.1)
 f <- GaussRF(n=3,  # three simulations at once
              x=x, y=y, model=model, grid=TRUE,  
              param=c(mean, variance, nugget, scale, alpha))
 image(x, y, f[,,1])
 image(x, y, f[,,2])
 image(x, y, f[,,3])
        

 #############################################################
 ##                                                         ##
 ##      Examples using the extended definition forms       ##
 ##                                                         ##
 ##                                                         ##  
 #############################################################

## note that the output seems plausible but not checked!!!!

## tbm may also be used for multiplicate models (if they have
## *exactly* the same anisotropy parameters)
x <- (0:100)/10
m <- matrix(c(1,2,3,4),ncol=2)/5
z <- GaussRF(x=x, y=x, grid=TRUE,
              model=list(
                list(m="power",v=1,k=2,a=m),
                "*", list(m="sph", v=1, a=m)
                ),
              me="TBM3", reg=0,n=1)
print(c(mean(as.double(z)),var(as.double(z))))
image(z,zlim=c(-3,3))


## non-separable space-time model applied for two space dimensions
## note that tbm method does not work nicely, but at least
## in some special cases.
x <- y <- (1:32)/2     ## grid definition, but as a sequence
T <- c(1,32,1)*10      ## note necessarily gridtriple definition
aniso <- diag(c(0.5,8,1)) 
k <- c(1,phi=1,1,1,psi=1,dim=2)
model <- list(list(m="nsst", v=1, k=k, a=aniso))
z <- GaussRF(x=x, y=y, T=T, grid=TRUE, model=model)
for (i in 1:dim(z)[3]) { image(z[,,i]); readline();}
for (i in 1:dim(z)[2]) { image(z[,i,]); readline();}
for (i in 1:dim(z)[1]) { image(z[i,,]); readline();}


 #############################################################
 ##                                                         ##
 ##                    Brownian motion                      ##
 ##                  using Stein's method                   ##
 ##                                                         ##  
 #############################################################
 # 2d
 x <- (0:100)/10
 kappa <- 1   # in [0,2)
 z <- GaussRF(x=x, y=x, grid=TRUE, model="2d", param=c(0,1,0,1,kappa))
 image(z,zlim=c(-3,3))

 # 3d
 x <- (0:30)/10
 kappa <- 1   # in [0,2)
 z <- GaussRF(x=x, y=x, z=x, grid=TRUE, model="3d",
              param=c(0,1,0,1,kappa))
 for (i in 1:dim(z)[1]) { image(z[i,,]); readline();}



 #############################################################
 ## This example shows the benefits from stored,            ##
 ## intermediate results: in case of the circulant          ##
 ## embedding method, the speed is doubled in the second    ##
 ## simulation.                                             ##  
 #############################################################
 DeleteAllRegisters()
 RFparameters(Storing=TRUE,PrintLevel=1)
 y <- x <- seq(0, 50, 0.2)
 (p <- c(runif(3), runif(1)+1))
 ut <- unix.time(f <- GaussRF(x=x,y=y,grid=TRUE,model="exponen",
                              method="circ", param=p))
 image(x, y, f)
 hist(f)
 c( mean(as.vector(f)), var(as.vector(f)) )
 cat("unix time (first call)", format(ut,dig=3),"\n")

 # second call with the *same* parameters is much faster:
 ut <- unix.time(f <- GaussRF(x=x,y=y,grid=TRUE,model="exponen",
                              method="circ",param=p)) 
 image(x, y, f)
 hist(f)
 c( mean(as.vector(f)), var(as.vector(f)) )
 cat("unix time (second call)", format(ut,dig=3),"\n")



%##############################################################
%##       for a further, more complicated example on         ##
%##        space-time modelling, see data(winddata)          ##
%##############################################################


}
\keyword{spatial}


