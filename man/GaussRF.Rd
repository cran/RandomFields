\name{GaussRF}
\alias{InitGaussRF}
\title{Gaussian Random Fields}
\description{
  These functions simulate isotropic Gaussian random fields
  using turning bands, circulant embedding, direct methods,
  and the random coin method.
}
\usage{
GaussRF(x, y=NULL, z=NULL, grid, model, param, method=NULL,
          n=1, register=0, gridtriple=FALSE)

InitGaussRF(x, y=NULL, z=NULL, grid, model, param,
               method=NULL, register=0, gridtriple=FALSE)
}
\arguments{
\item{x}{matrix of coordinates, or vector of x coordinates}
\item{y}{vector of y coordinates}
\item{z}{vector of z coordinates}
\item{grid}{logical; determines whether the vectors \code{x},
  \code{y}, and \code{z} should be
  interpreted as a grid definition, see Details.}
\item{model}{string; covariance or variogram model,
  see \code{\link{CovarianceFct}}, or
  type \code{\link{PrintModelList}()} to get all options}
\item{param}{parameter vector:
  \code{param=c(mean, variance, nugget, scale,...)};
  the parameters must be given
  in this order; further parameters are to be added in case of a
  parametrised class of models, see \code{\link{CovarianceFct}}}
\item{method}{\code{NULL} or string; Method used for simulating,
  see \code{\link{RFMethods}}, or
  type \code{\link{PrintMethodList}()} to get all options}
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
}\details{
  \code{GaussRF} creates an isotropic Gaussian random field with
  \code{model} as covariance function/variogram model and parameters
  \code{param=c(mean,variance,nugget,scale,...)}. 
  The sill of the variogram equals \code{variance + nugget}.
  
  \code{GaussRF} can use different methods for the simulation,
  i.e., circulant embedding, turning bands, direct methods, and random
  coin method. 
  If \code{method==NULL} then \code{GaussRF} searches for a
  valid method.  \code{GaussRF} may not find the fastest method neither the
  most precise one.  It just finds any method among the available
  methods. (However it guesses what is a good choice.) 
  Note that some of the methods do not work for all covariance 
  or variogram models.
  
  \code{GaussRF} is split up in an initial \code{InitGaussRF},
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
  
  Using \code{InitGaussRF} and \code{\link{DoSimulateRF}} in sequence might
  be slightly faster than \code{GaussRF} (but less convenient). 
  
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
  \item \code{model="nugget"} is one possibility to create independent
  Gaussian random variables.  Without loss of efficiency, any
  covariance function with parameter vector
  \code{c(mean, 0, nugget, scale, ...)}
  can also be used.  If \code{model="nugget"} is used
  the second component of \code{param} must be zero.
  %% this has to be changed in later versions; 

  \item The sum of the components variance and nugget in the argument
  \code{param} equals the sill of the variogram. 

  \item \code{register} is a parameter which may never be
  used by most users (please let me know if you use it!).  In
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
  \code{InitGaussRF} returns 0 if no error has occured and a positive value
  if failed.\cr

  \code{GaussRF} and \code{\link{DoSimulateRF}} return \code{NULL}
  if an error has occured; otherwise the returned object
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
  Gneiting, T. and Schlather, M. (2001)
  Statistical modeling with covariance functions.
  \emph{In preparation}.
  
  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields.}
  Technical report ST 99-10, Dept. of Maths and Statistics,
  Lancaster University. 
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
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
  \code{\link{ShowModels}}.
}

\examples{
 #############################################################
 ## Examples using the symmetric stable model, also called  ##
 ## "powered exponential model"                             ## 
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
 f <- GaussRF(grid=FALSE,
              x=x, y=y, z=z, model=model, 
              param=c(mean, variance, nugget, scale, alpha))
 f


 #############################################################
 ## usage of a specific method
 ## -- the complete list can be obtained by PrintMethodList()
 x <- runif(100, max=20) # arbitrary points
 y <- runif(100, max=20)
 f <- GaussRF(method="dir",  # direct matrix decomposition
              x=x, y=y, model=model, grid=FALSE, 
              param=c(mean, variance, nugget, scale, alpha))
 f


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
 ## This example shows the benefits from stored,            ##
 ## intermediate results: in case of the circulant          ##
 ## embedding method, the speed is doubled in the second    ##
 ## simulation.                                             ##  
 #############################################################
 DeleteAllRegisters()
 RFparameters(Storing=TRUE,PrintLevel=1)$null
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
}
\keyword{spatial}


