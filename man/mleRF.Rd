\name{fitvario}
\alias{mleRF}% obsolete
\alias{fitvario}
\alias{fitvario.default}
\title{LSQ and Maximum Likelihood Estimation
  of Random Field Parameters}
\description{
  This function estimates arbitrary parameters of
  a random field specification.
  }
\usage{
fitvario(x, y=NULL, z=NULL, T=NULL, data, model, param,
         lower=NULL, upper=NULL, sill=NA, ...)

fitvario.default(x, y=NULL, z=NULL, T=NULL, data, model, param,
         lower=NULL, upper=NULL, sill=NA, trend,
         use.naturalscaling=TRUE, PrintLevel=RFparameters()$Print,
         trace.optim=0, bins=20, nphi=1, ntheta=1, ntime=20,
         distance.factor=0.5,
         upperbound.scale.factor=10, lowerbound.scale.factor=20,
         lowerbound.scale.LS.factor=5,
         upperbound.var.factor=10, lowerbound.var.factor=100,
         lowerbound.sill=1E-10, scale.max.relative.factor=1000,
         minbounddistance=0.001, minboundreldist=0.02,
         approximate.functioncalls=50, pch="*", var.name="X",
         time.name="T", transform=NULL, standard.style=NULL)
}
\arguments{
  \item{x}{(\eqn{n\times2}{n x 2})-matrix of coordinates, or vector of
    x-coordinates}
  \item{y}{vector of y coordinates}
  \item{z}{vector of z coordinates}
  \item{T}{vector of T coordinates; these coordinates are given in
    triple notation, see \code{\link{GaussRF}}}
  \item{data}{vector or matrix of values measured at \code{coord};
    If also a time component is given, the in the data the indices for
    the spatial components run the fastest.
  }
  \item{model}{string or list;
    covariance model, see \code{\link{CovarianceFct}}, or
    type \code{\link{PrintModelList}()} to get all options.
    See also t
    
    If \code{model} is a list, then the parameters with value \code{NA}
    are estimated. Parameters that have value \code{NaN} should be
    explicitely be defined by the function \code{transform}.
    An alternative to define \code{NaN} values and
    the function \code{transform}, is to replace the \code{NaN}
    by a real-valued function with solely parameter a list defining
    a covariance model. In case of the anisotropy matrix, the matrix
    must be replaced by a list if functions are introduced.
    Only the list elements
    variance, scale or anisotropy, and kappas can be used, and
    not the mean or the trend.
    Further, the mean or the trend cannot be set by such a function.
    See also \code{transform} below.
  }
  \item{param}{
    vector or matrix or NULL.
    If vector then
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order. Further parameters are to be added in case of a
    parametrised class of covariance functions, see
    \code{\link{CovarianceFct}}.
    \emph{Any} components set to \code{NA} are estimated; the others
    are kept fix.
    See also \code{model} above. 
  }
  \item{lower}{list or vector. Lower bounds for the parameters.
    If \code{lower} and \code{param} are vectors and \code{length(lower) <
      length(param)} 
    then \code{lower} must match the number of additional parameters
    \eqn{\kappa}{a,b,c,...}.
    
    If \code{param} is matrix the length of \code{lower} must match the
    number columns of param or being 2 elements smaller (then
    \code{lower} is filled with \code{NA} from the left.
    The bounds are equally applied to all rows.

    If \code{lower} is a list, then elements that are not given are
    considered as \code{NA}.

    If \code{lower} is not given, or \code{lower} contains \code{NA}
    then the missing bounds are generated automatically.
  }
  \item{upper}{list or vector. Upper bounds for the parameters.
    See also lower.
  }
  \item{sill}{If not \code{NA} the sill is kept fix. Only used if the
    standard format for the covariance model is given. See Details.}
  \item{trend}{
    Not programmed yet.
    May only be set if \code{missing(param)};
    linear formula : uses \code{X1}, \code{X2},... and \code{T} as internal
    parameters for the coordinates; all parameters are estimated\cr
    matrix   : must have the same number of rows as \code{x}\cr
    \emph{fixed} mean + matrix or linear formula : not possible within
    this function (just subtract the mean from your data before calling
    this function)    
  }
  \item{...}{arguments as given in \code{mleRF.default} and listed in the
    following.}
  \item{use.naturalscaling}{
    logical.  Only used if model is given in standard (simple) way.
    If \code{TRUE} then \emph{internally}, rescaled
    covariance functions will be used for which
    cov(1)\eqn{\approx}{~=}0.05.
    \code{use.naturalscaling} has the advantage that \code{scale}
    and the form parameters of the model get \sQuote{orthogonal},
    but \code{use.naturalscaling} does not work for all models.
    See Details.}
  \item{PrintLevel}{level to which messages are shown.  See Details.}
  \item{trace.optim}{tracing of the function \code{\link{optim}}}
  \item{bins}{number of bins of the empirical variogram.  See Details.}
  \item{nphi}{scalar or vector of 2 components.
    If it is a vector then the first component gives the first angle
    of the xy plane
    and the second one gives the number of directions on the half circle.
    If scalar then the first angle is assumed to be zero}
  \item{ntheta}{scalar or vector of 2 components.
    If it is a vector then the first component gives the first angle
    in the third direction
    and the second one gives the number of directions on the half circle.
    If scalar then the first angle is assumed to be zero.}
  \item{ntime}{scalar or vector of 2 components.
    if \code{ntimes} is a vector, then the first component are the
    maximum time distance (in units of the grid length \code{T[3]}) and the
    second component gives the step size (in units of the grid length
    \code{T[3]}). If scalar then the step size is assumed to 1 (in units
    of the grid length \code{T[3]}.
  }
  \item{distance.factor}{relative right bound for the bins.  See Details.}
  \item{upperbound.scale.factor}{relative upper bound for scale in LSQ
    and MLE. 
    See Details.}
  \item{lowerbound.scale.factor}{relative lower bound for scale in MLE. 
    See Details.}
  \item{lowerbound.scale.LS.factor}{relative lower bound for scale
    in LSQ.  See Details.}
  \item{upperbound.var.factor}{relative upper bound for variance
    and nugget.  See Details.}
  \item{lowerbound.var.factor}{relative lower bound for variance.
    See Details.}
  \item{lowerbound.sill}{absolute lower bound for variance
    and nugget.  See Details.}
  \item{scale.max.relative.factor}{relative lower bound
    for scale below
    which an additional nugget effect is detected.  See Details.}
  \item{minbounddistance}{absolute distance to the bounds
    below which a part of the algorithm is considered as
    having failed.  See Details.}
  \item{minboundreldist}{relative distance to the bounds
    below which a part of the algorithm is considered as
    having failed.  See Details.}
  \item{approximate.functioncalls}{approximate evaluations
    of the ML target function on a grid.  See Details.}
  \item{pch}{character shown before each step of calculation;
    depending on the specification there are two to five steps.
    Default: "*".}
  \item{var.name}{basic name for the coordinates in the formula of the
    trend.  Default: \sQuote{X}}
  \item{time.name}{basic name for the time component in the formula of the
    trend.  Default: \sQuote{X}}
  \item{transform}{vector of strings.
    Essentially, \code{transform} allows for the definition of a parameter as a
    function of other estimated or fixed parameters.
    All the parameters are supposed to be in a vector called \sQuote{param}
    where the positions are given by \code{\link{parampositions}}.
    An example of \code{transform} is
    \code{function(param) \{param[3] <- 5 - param[1]; param\}}.
    Any parameter that is set by \code{transform}, should be \code{NaN}
    in the model definition. If it is \code{NA} a warning is given.
    Note that the mean and the trend of the model can be neither set nor used in
    \code{transform}.  See also \code{standard.style}.
    
    Instead of giving \code{transform}, in the model definition,
    all NaN values are replaced by functions whose only parameter
    is a bare model list, i.e., only the list elements
    variance, scale or anisotropy, and kappas can be used, and
    not the mean or the trend.
    Further, the mean or the trend cannot be set by such a function.
    
    Default: \code{NULL}}
  \item{standard.style}{logical or \code{NULL}.  This variable should only be
    set by the advanced user.  If \code{NULL} \code{standard.style} will be
    \code{TRUE} if the covariance model allows for a \sQuote{standard}
    definition (see \command{\link{convert.to.readable}} and
    \command{\link{CovarianceFct}}) and \code{transform} is \code{NULL}.
    Otherwise \code{standard.style} will be \code{FALSE}.
    If a \\sQuote{standard} definition is given and both the variance and the
    nugget are either not estimated or do not appear on the right hand
    side of the \code{transform}, then \code{standard.style} might be set
    to \code{TRUE} by the user. This accelerates the MLE algorithm.
    The responsibility is completely left to the user, then.

    Currently \command{mleRF} is only implemented for the
    \\sQuote{standard} definition of the covariance model.
    Hence \code{standard.style} must always be \code{TRUE} and consequently,
    neither the variance nor the nugget may appear on either side
    of the \code{transform}
  }
}

\details{
  The maximisation is performed using \code{\link{optim}}.  Since
  \code{\link{optim}}
  needs as input parameter an initial vector of parameters, \code{mleRF}
  takes the initial parameter from the LSQ estimation.

  If the best parameter vector of the MLE found so far is too close
  to some given bounds, see the specific parameters below, it is
  assumed that
  \code{\link{optim}} ran into a local minimum because of a bad starting
  value.
  In this case the MLE target function is calculated on a grid, the
  best parameter vector is taken, and the optimisation is restarted with
  this parameter vector.
  
  Comments on specific parameters:
  \itemize{
    \item \code{lower}\cr
    The lower bounds are technical bounds that should not
    restrict really restrict the domaine of the value.
    However, if these values are too small the optimisation
    algorithm will frequently run into local minima or get
    stuck close the border of the parameter domain.
    It is advised to limit seriously the domain of the additional
    parameters of the covariance model and/or the total number of
    parameters to be estimated, if \dQuote{many} parameters
    of the covariance model are  estimated.
    
    If the model is given in standard form, the user may supply
    the lower bounds for the whole parameter vector, or only for
    the additional form parameters of the model.
    The lower bound for the mean will be ignored.
    \code{lower} may contain \code{NA}s, then these values
    are generated by the 

    If a nested model is given, the bounds may again be supplied
    for all parameters or only for the additional form parameters
    of the model. The bounds given apply uniformely to all
    submodels of the nested model.

    If the \code{model} is given in list format, then
    \code{lower} is a list, where components may be missing
    or \code{NA}. These are generated by the algorithm, then.

    If \code{lower} is \code{NULL} all lower bounds are generated
    automatically. 
     
    \item \code{upper.kappa}\cr
    See \code{lower.kappa}.
    
    \item \code{sill}\cr
    Additionally to estimating \code{nugget} and \code{variance}
    separately, they may also be estimated together under the
    condition that \code{nugget} + \code{variance} = \code{sill}.
    For the latter a finite value for \code{sill} has to be supplied,
    and  \code{nugget} and  \code{variance} are set to \code{NA}.

    \code{sill} is only used for the standard model. 
    
   }\itemize{
   
    \item \code{use.naturalscaling}\cr
    logical. If \code{TRUE} then internally, rescaled
    covariance functions will be used for which
    cov(1)\eqn{\approx}{~=}0.05. However this parameter
    does not influence
    the output of \code{mleRF}: the parameter vector
    returned by \code{mleRF} refers
    \emph{always} to the standard covariance model as given in
    \code{\link{CovarianceFct}}. (In contrast to \code{PracticalRange}
    in \code{\link{RFparameters}}.)\cr
    Advantages if \code{use.naturalscaling==TRUE}:
    \itemize{
      \item \code{scale} and the shape parameter of a parameterised
      covariance model can be estimated better if they are estimated
      simultaneously.
      \item The estimated bounds calculated  by means of
      \code{upperbound.scale.factor} and \code{lowerbound.scale.factor},
      etc. might be more realistic.
      \item in case of anisotropic models, the inverse of the elements
      of the anisotropy matrix should be in the above bounds.
    }
    Disadvantages if \code{use.naturalscaling==TRUE}:
    \itemize{
      \item For some covariance models with additional parameters, the
      rescaling factor has to be determined numerically.
      Then, more time is needed to perform \code{mleRF}.
    }
    Default: \code{TRUE}.
    
    \item \code{PrintLevel}\cr
    0 : no message\cr
    1 : error messages\cr
    2 : warnings\cr
    3 : minimum debugging information\cr
    5 : extended debugging information, including graphics\cr
    Default: \code{0}.
    
   }\itemize{
   \item \code{trace.optim}\cr
    see control parameter \code{trace} of
    \code{\link{optim}}.  Default: \code{0}.

    \item \code{bins}\cr
    vector of explicit boundaries for the bins or the
    number of bins for the empirical variogram (used in the
    LSQ target function, which is described at the beginning
    of the Details).  Default: \code{20}.
    
    \item \code{distance.factor}\cr
    right boundary of the last
    bin is calculated as \code{distance.factor} * (maximum distance
    between all pairs of points). Only used if \code{bins} is a scalar.
    Default: \code{0.5}.
    
    \item \code{upperbound.scale.factor}\cr
    The upper bound for the scale is determined
    as \code{upperbound.scale.factor} * (maximum distance
    between all pairs of points). 
    Default: \code{10}.
    
    \item \code{lowerbound.scale.factor}\cr
    The lower bound for the scale is determined
    as

    \eqn{\hbox{(minimum distance
      between different pairs of points)} /
      \code{lowerbound.scale.factor}}{(minimum distance
      between different pairs of points) /
      \code{lowerbound.scale.factor}}.
    
    Default: \code{20}.\cr
    
  }\itemize{
    \item \code{lowerbound.scale.LS.factor}\cr
    For the LSQ target function a different lower bound
    for the scale is used. It is determined as
    
    \eqn{\hbox{(minimum distance
    between different pairs of points)} /
    \code{lowerbound.scale.LS.factor}}{(minimum distance
    between different pairs of points) /
    \code{lowerbound.scale.LS.factor}}.
    
    Default: \code{5}.\cr
    
    \item \code{upperbound.var.factor}\cr
    The upper bound for the variance and the nugget is determined
    as  
    \deqn{\code{upperbound.var.factor} *
      {\hbox{var}}(\code{data}).}{\code{upperbound.var.factor} *
      var(\code{data}).}
    Default: \code{10}.
    
  }\itemize{
    \item \code{lowerbound.var.factor}\cr
    The lower bound for the variance and the nugget is determined
    as
    \deqn{\hbox{var}(\code{data}) /
      \code{lowerbound.var.factor}.}{var(\code{data}) /
      \code{lowerbound.var.factor}.} 
    If a standard model definition is given and
    either the nugget or the variance is fixed,
    the parameter to be estimated
    must also be greater than \code{lowerbound.sill}.
    If a non-standard model definition is given
    then \code{lowerbound.var.factor} is only used
    for the first model; the other lower bounds for the
    variance are zero.
    Default: \code{100}.
    
    \item \code{lowerbound.sill}\cr
    See \code{lowerbound.var.factor}. 
    Default: \code{1E-10}.
    
    \item \code{scale.max.relative.factor}\cr
    If the initial scale value for the ML estimation
    obtained by the LSQ target function is
    less than

    \eqn{\hbox{(minimum distance
    between different pairs of points)} /
    \code{scale.max.relative.factor}}{ (minimum distance
    between different pairs of points) /
    \code{scale.max.relative.factor}}
    
    it is assumed that a nugget effect
    is present.  In case the user set \code{nugget=0},
    the ML estimation is automatically performed
    for \code{nugget=NA} instead of \code{nugget=0}. 
    Note: if \code{scale.max.relative.factor} is greater
    than \code{lowerbound.scale.LS.factor} then
    \code{nugget} is never set to \code{NA} as
    the scale has the lower bound (minimum distance
    between different pairs of points) /
    \code{lowerbound.scale.LS.factor}. 
    Default: \code{1000}.
    
  }\itemize{
    \item \code{minbounddistance}\cr
    If any value of the parameter vector
    returned from the ML estimation
    is closer than \code{minbounddistance}
    to any of the bounds or if any value
    has a relative distance smaller than
    \code{minboundreldist}, then it is assumed that
    the MLE algorithm has dropped into a local minimum,
    and it will be continued with evaluating the
    ML target function on a grid, cf. the beginning paragraphs
    of the Details.
    Default: \code{0.001}.
    
    \item \code{minboundreldist}\cr
    See \code{minbounddistance}. 
    Default: \code{0.02}.
    
    \item \code{approximate.functioncalls}\cr
    In case the parameter vector is too close to the given
    bounds, the ML target function is evaluated on a grid
    to get a new initial value for the ML estimation. 
    The number of points of the grid is approximately
    \code{approximate.functioncalls}. 
    Default: \code{50}.
  }
  Another maximum likelihood estimator for random fields
  exists as part of the package \code{geoR} whose homepage
  is at \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html},
  with a different philosophy behind.
}
\value{
  the function returns a list with the following elements
  \item{mle.value}: the log-likelihood for the fitted parameters
  \item{trend.coff}{parameters for linear trend (optional)}
  \item{mle}{fitted model}
  \item{ev}{list returned by \code{\link{EmpiricalVariogram}}}
  \item{lsq}{model fitted by least squares; trends are never taken into account}
  \item{nlsq}{weighted lsq.  Weight is the square root of the number
    of points in the bin}
  \item{slsq}{weighted lsq.  Weight is the inverse the standard deviation of the
    variogram cloud within the bin}
  \item{flsq}{weighted lsq.  Weights are the values of the fitted variogram
    to the power of -2}
  \code{mle.lower}{lower bounds for the parameters used in the
    optimisation algorithm} 
  \code{mle.upper}{upper bounds for the parameters used in the
    optimisation algorithm}
}
\references{
  Ribeiro, P. and Diggle, P. (2001) Software for geostatistical analysis
  using R and S-PLUS: geoR and geoS, version 0.6.15.
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}
\author{Martin Schlather, \email{martin.schlather@cu.lu}
  \url{http://www.cu.lu/~schlathe}}

\note{This function does not depend on the value of
  \code{\link{RFparameters}()$PracticalRange}. 
  The function \code{mleRF} always uses the standard specification
  of the covariance model as given in \code{\link{CovarianceFct}}.
  
%  The function is still in an experimental stage (and will probably
%  never leave it).  The function works well if
%  only a few components of the parameter vector are to be estimated.
}
\section{Acknowledgement}{
  Thanks to Paulo Ribeiro for hints and comparing \code{mleRF} to
  \code{\link[geoR]{likfit}} of the package \code{geoR} whose homepage is at 
  \url{http://www.est.ufpr.br/geoR/}.
}
\seealso{
  \code{\link{CovarianceFct}},
  \code{\link{GetPracticalRange}},
  \command{\link{parampositions}}
  \code{\link{RandomFields}},
  \code{\link{winddata}}.
}
\examples{
 % library(RandomFields, lib="~/TMP"); source("~/article/R/NEW.RF/RandomFields/R/MLE.R");

 model <-"gencauchy"
 param <- c(0, 1, 0, 1, 1, 2)
 estparam <- c(0, NA, 0, NA, NA, 2) ## NA means: "to be estimated"
 ## sequence in `estparam' is
 ## mean, variance, nugget, scale, (+ further model parameters)
 ## So, mean, variance, and scale will be estimated here.
 ## Nugget is fixed and equals zero.
 points <- 100
 x <- runif(points,0,3)
 y <- runif(points,0,3) ## 100 random points in square [0, 3]^2
 d <- GaussRF(x=x, y=y, grid=FALSE, model=model, param=param, n=10)
 str(fitvario(x=cbind(x,y), data=d, model=model, param=estparam,
              lower=c(0.1, 0.1), upper=c(1.9, 5)))


## The next two estimations give about the same result.
## For the first the sill is fixed to 1.5. For the second the sill
## is reached if the estimated variance is smaller than 1.5
estparam <-  c(0, NA, NA, NA, NA, NA) 
str(fitvario(x=cbind(x,y), data=d, model=model, param=estparam, sill=1.5))

estparam <-  c(0, NA, NaN, NA, NA, NA) 
parampositions(model=model, param=estparam)
f <- function(param) {
   param[5] <- max(0, 1.5 - param[1])
   return(param)
}
str(fitvario(x=cbind(x,y), data=d, model=model, param=estparam, 
    sill=1, transform=f))

## the next call gives a warning, since the user may programme
## strange things in this setup, and the program cannot check it.
estparam <- c(0, NA, NA, NA, NA, NaN) 
parampositions(model=model, param=estparam)
f <- function(param) {param[3] <- param[2]; param}
unix.time(str(fitvario(x=cbind(x,y), data=d, model=model,
      param=estparam, transform=f, standard.style=TRUE), vec.len=6))

## much better programmed, but also much slower:
estmodel <- list(list(model="gencauchy", var=NA, scale=NA,
                      kappa=list(NA, function(m) m[[1]]$kappa[1])))
unix.time(str(fitvario(x=cbind(x,y), data=d, model=estmodel),
              vec.len=6))


%#################################################################
%##                                                             ##
%## for a further, more complicated example, see data(winddata) ##
%##                                                             ##
%#################################################################
%
% check space time extension in MLE\n add examples\n add test

}
\keyword{spatial}







%  LocalWords:  fitvario mleRF LSQ param NA naturalscaling PrintLevel optim pch
%  LocalWords:  RFparameters nphi ntheta ntime lowerbound minbounddistance var
%  LocalWords:  minboundreldist functioncalls eqn coord CovarianceFct emph cr
%  LocalWords:  PrintModelList cov variogram uote parampositions NaN kappas MLE
