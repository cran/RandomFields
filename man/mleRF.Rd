\name{mleRF}
\alias{mleRF}
\title{Maximum Likelihood Estimation of Random Field Parameters}
\description{
  This function estimates arbitrary parts of the parameter vector of
  a random field specification by maximising the likelihood.
}
\usage{
mleRF(coord, data, model, param, lower.kappa=NULL, upper.kappa=NULL,
      sill=NA, ...) }
\arguments{
  \item{coord}{(\eqn{n\times2}{n x 2})-matrix of coordinates.}
  \item{data}{vector of values measured at \code{coord}.}
  \item{model}{string; covariance model, see \code{\link{CovarianceFct}}, or
  type \code{\link{PrintModelList}()} to get all options.}
  \item{param}{vector of parameters for a random field specification:
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order. Further parameters are to be added in case of a
    parametrised class of covariance functions, see
    \code{\link{CovarianceFct}}.
    \emph{Any} components set to \code{NA} are estimated; the others
    are kept fix.
  }
  \item{lower.kappa}{lower bounds for additional parameters of a
    parametrised class.}
  \item{upper.kappa}{upper bounds for additional parameters of a
    parametrised class.}
  \item{sill}{If not \code{NA} the sill is kept fix.  See Details.}
  \item{use.naturalscaling}{
    logical.  If \code{TRUE} then \emph{internally}, rescaled
    covariance functions will be used for which
    cov(1)\eqn{\approx}{~=}0.05.  See Details.}
  \item{PrintLevel}{level to which messages are shown.  See Details.}
  \item{trace.optim}{tracing of the function \code{optim}}
  \item{bins}{number of bins of the empirical variogram.  See Details.}
  \item{distance.factor}{relative right bound for the bins.  See Details.}
  \item{upperbound.scale.factor}{relative upper bound for scale. 
    See Details.}
  \item{lowerbound.scale.factor}{relative lower bound for scale. 
    See Details.}
  \item{lowerbound.scale.LS.factor}{relative lower bound for scale
    in an auxiliary function call.  See Details.}
  \item{upperbound.var.factor}{relative upper bound for variance
    and nugget.  See Details.}
 \item{lowerbound.var.factor}{relative lower bound for variance
    and nugget.  See Details.}
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
}
\details{
  The maximisation is performed using \code{optim}.  Since \code{optim}
  needs as input parameter an initial vector of parameters, \code{mleRF}
  estimates this initial parameter first using an auxiliary target
  function, namely weighted least squares for the empirical variogram.
  
  If the best parameter vector of the MLE found so far is too close
  to some given bounds, see the specific parameters below, it is
  assumed that
  \code{optim} ran into a local minimum because of a bad starting
  value by the auxiliary target function.
  In this case the MLE target function is calculated on a grid, the
  best parameter vector is taken, and the optimisation is restarted with
  this parameter vector.
  
  Comments on specific parameters:
  \itemize{
    \item \code{lower.kappa}\cr
    If a covariance function possesses additional parameters, see
    \code{stable} in \code{\link{CovarianceFct}} for example,
    one may supply (vectors of) lower and upper limits by
    \code{lower.kappa} and \code{upper.kappa}. 
    Note that the values of these two parameters are \emph{not}
    recycled if
    the number of parameters of a covariance function
    is a multiple of the length of \code{lower.kappa} or
    \code{upper.kappa}. 
    In case a covariance model has more than one parameter, and
    if the second parameter is to be estimated, say, then limits
    for the first \emph{and} the second parameters should be supplied. 
    (The function assigns the given limits to the parameters by
    always starting with the first parameter, whether it is estimated
    or not.) 
    It is advised to limit seriously the domain of the additional
    parameters of the covariance model and/or the total number of
    parameters to be estimated, if additional parameters
    of the covariance model are also estimated.
    
    \item \code{upper.kappa}\cr
    See \code{lower.kappa}.
    
    \item \code{sill}\cr
    Additionally to estimating \code{nugget} and \code{variance}
    separately, they may also be estimated together under the
    condition that \code{nugget} + \code{variance} = \code{sill}.
    For the latter a finite value for \code{sill} has to be supplied,
    and  \code{nugget} and  \code{variance} are set to \code{NA}.
    
    
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
      \code{upperbound.scale.factor}, \code{lowerbound.scale.factor},
      etc. might be more realistic.
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
    
    \item \code{trace.optim}\cr
    see control parameter \code{trace} of
    \code{\link{optim}}.  Default: \code{0}.

    \item \code{bins}\cr
    number of bins for the empirical variogram (used in the
    auxiliary target function, which is described at the beginning
    of the Details).  Default: \code{20}.
    
    \item \code{distance.factor}\cr
    right boundary of the last
    bin is calculated as \code{distance.factor} * (maximum distance
    between all pairs of points). 
    Default: \code{0.5}.
    
    \item \code{upperbound.scale.factor}\cr
    The upper bound for the scale is determined
    as \code{upperbound.scale.factor} * (maximum distance
    between all pairs of points). 
    Default: \code{10}.
    
    \item \code{lowerbound.scale.factor}\cr
    The lower bound for the scale is determined
    as (minimum distance
    between different pairs of points) /
    \code{lowerbound.scale.factor}. 
    Default: \code{20}.
    
    \item \code{lowerbound.scale.LS.factor}\cr
    For the auxiliary target function a different lower bound
    for the scale is used. It is determined as (minimum distance
    between different pairs of points) /
    \code{lowerbound.scale.LS.factor}. 
    Default: \code{5}.\cr
    
    \item \code{upperbound.var.factor}\cr
    The upper bound for the variance and the nugget is determined
    as \code{upperbound.var.factor} * var(data). 
    Default: \code{10}.
    
    \item \code{lowerbound.var.factor}\cr
    The lower bound for the variance and the nugget is determined
    as var(data) / \code{lowerbound.var.factor}. 
    If either the nugget or the variance is fixed,
    the parameter to be estimated
    must also be greater than \code{lowerbound.sill}. 
    Default: \code{100}.
    
    \item \code{lowerbound.sill}\cr
    See \code{lowerbound.var.factor}. 
    Default: \code{1E-10}.
    
    \item \code{scale.max.relative.factor}\cr
    If the initial scale value for the ML estimation
    obtained by the auxiliary target function is
    less than (minimum distance
    between different pairs of points) /
    \code{scale.max.relative.factor}
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
  is at \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}. 
  However, the philosophies of \code{mleRF} and \code{likfit} of the
  package \code{geoR}
  differ.  The function \code{likfit} offers more possibilities and higher
  flexibility concerning co-variates and value transformations,
  for example.  By way of contrast, \code{mleRF} is still restricted to
  estimating parts of the parameter vector of a random field (as a
  consistent part of the package \code{RandomFields}).  However, \code{mleRF}
  allows for estimating any combination of the components of this
  parameter vector.  In case both functions can be used,
  the parameters estimated
  by \code{likfit} and \code{mleRF} agree in about 95
  percent of the cases.  Concerning the remaining cases, \code{likfit}
  beats \code{mleRF} more frequently than vice versa; however
  \code{mleRF} does not need initial values for the optimisation,
  and the code of \code{mleRF} is in general faster.
}
\value{
  A vector of the form
  \code{c(mean,variance,nugget,sill,...)} is returned, i.e.
  the parameter vector of the ML estimation including the fixed
  components. 
  Here, the dots \code{...}
  stand for additional parameters of the covariance model, see
  \code{\link{CovarianceFct}}.
}
\references{
  Ribeiro, P. and Diggle, P. (2001) Software for geostatistical analysis
  using R and S-PLUS: geoR and geoS, version 0.6.15.
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
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
  \code{likfit} of the package \code{geoR} whose homepage is at 
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}
\seealso{
  \code{\link{CovarianceFct}},
  \code{\link{GetPracticalRange}},
  \code{\link{RandomFields}}.
}
\examples{
 model <-"expon"
 param <- c(0,1,0,1)
 estparam <- c(NA,NA,0,NA) ## NA means here: "to be estimated"
 ## sequence in `estparam' is
 ## mean, variance, nugget, scale, (+ further model parameters)
 ## So, mean, variance, and scale will be estimated here.
 ## Nugget is fixed and equals zero.
 points <- 100
 x <- runif(points,0,3)
 y <- runif(points,0,3) ## 100 random points in square [0, 3]^2
 f <- GaussRF(x=x, y=y, grid=FALSE, model=model, param=param)
 mleRF(coord=cbind(x,y), data=f, model=model, param=estparam)
}
\keyword{spatial}







