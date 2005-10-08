\name{RFparameters}
\alias{RFparameters}
\title{Control Parameters}
\description{
  \code{RFparameters} sets and returns control parameters for the simulation
  of random fields
}
\usage{
   RFparameters(..., no.readonly=FALSE)
}
\arguments{
  \item{...}{arguments in \code{tag = value} form, or a list of tagged
    values.}
  \item{no.readonly}{If \command{RFparameters} is called without
    parameter then all parameters are returned in a list. If
    \code{no.readonly=TRUE} then only rewritable parameters are returned.
  }
}
\details{
  The possible parameters are
  %
  \cr\cr\bold{General options}\cr
  \describe{
    \item{\code{PracticalRange}}{logical or 2, 3, 11, 12, 13.
      If not \code{FALSE} the range of the covariance functions is
      adjusted so that cov(1) is about 0.05 (for \code{scale==1}).
      \itemize{
	\item \code{FALSE} : the practical range ajustment is not used.
	\item \code{TRUE} : \code{PracticalRange} is applicable only if
	the value is known exactly.
	\item \code{2} : \code{PracticalRange} is applicable if the value is
	known pretty well
	\item \code{3} : \code{PracticalRange} is applicable if the value is
	roughly known
	\item \code{11} : if the practical range is not known exactly it
	is approximated numerically.
	\item \code{12} : if the practical range is not known pretty well it
	is approximated numerically.
	\item \code{13} : if the practical range is even not known
	approximately it is approximated numerically.
      }
      Note that values beyond \code{FALSE}, \code{TRUE}, and \code{11},
      are only used for specialists' purposes.
      
      Default: \code{FALSE} [init].}
    
    \item{\code{PrintLevel}}{If \code{PrintLevel}\eqn{\le0}{<=0}
      there is not any output on the screen.  The
      higher the number the more tracing information is given. 
      Default: 1 [init, do].\cr
      1 : error messages\cr
      2 : messages about partial failures of the algorithm\cr
      >2 : additional informations
      
      Note that \code{PrintLevel} is also used in other packages
      as a default, for example in \link[SoPhy]{SoPhy}
      (\command{\link[SoPhy]{risk.index}} and
      \command{\link[SoPhy]{create.roots}}). The changing of
      \code{PrintLevel} here may cause some unexpected effects in these
      functions. See the documentation there.
    }
  }
  %
  \bold{General options for simulating}\cr
  \describe{
    \item{\code{pch}}{Character or empty string.
      The character is printed after each
      performed simulation if more than one simulation is performed at
      once. If \code{pch='!'} then a counter is shown instead of the character.
      Note that also '\eqn{\mbox{\textasciicircum}}{^}H's are printed if
      the counter (\code{pch='!'}) is shown, 
      which may have undesirable interactions with some few other R
      functions, e.g.  \command{\link[utils]{Sweave}}.
      Default: \code{'*'} [do]. 
    }
   \item{\code{Storing}}{Logical.
      If \code{TRUE} then intermediate results are kept
      after each simulation; if several simulation are performed with the same
      model parameters then
      \code{Storing=TRUE} accelerates the simulations, but needs additional
      memory. \cr
      Note that in subsequent calls of \command{\link{GaussRF}}
      intermediate changes of RFparameters with flag "[init]"
      do not have any influence on the algorithm.
      Hence, for studying the effects for
      divers values of technical parameters like
      \code{CE.force}, \code{CE.mmin}, etc. the parameter \code{Storing} must be
      \code{FALSE}.  See also last paragraphs in the Details.
      Default: \code{FALSE} [init, do].
    }
    \item{\code{stationary.only}}{Logical or NA. Used for the automatic
      choice of methods. See also \link{RFMethods}.
      \itemize{
	\item \code{TRUE}: the simulation of non-stationary random
	fields is refused. In particular, the intrinsic
	embedding method is excluded and
	the simulation of Brownian motion is rejected.
	\item \code{FALSE}: intrinsic embedding is always allowed,
	actually it's the first one considered in the automatic
	selection algorithm.
	\item \code{NA}: the simulation of the Brownian motion allowed,
	but intrinsic embedding is not used for stationary random fields.
      }
      Default: \code{NA} [init].}
    \item{\code{exactness}}{logical or NA.
      \itemize{
	\item \code{TRUE}: add.MPP, hyperplanes and
	all turning bands methods are excluded.
	If the circulant embedding method is considered as badly
	behaved, then the matrix decomposition methods are preferred.
	\item \code{FALSE}: if the circulant embedding method is
	considered as badly behaved or the number of points to be
	simulated is large, the turning bands methods are
	rather preferred.
	\item \code{NA}: approximative non-exact methods are excluded,
	i.e. TBM2 if the Abel transform of the covariance function
	cannot be given explicitely.
      }
      Default: \code{NA} [init].
    }
  }
  %
  
  \bold{Options for simulating with the circulant embedding method}\cr
  \describe{
    \item{\code{CE.force}}{Logical.  Circulant embedding does not work if a
      certain circulant
      matrix has negative eigenvalues.  Sometimes it is convenient
      to replace all the negative eigenvalues by zero
      (\code{CE.force==TRUE}) after \code{CE.trials} number of trials. 
      Default: \code{FALSE} [init].
    }
    \item{\code{CE.mmin}}{Scalar or vector, integer if positive.
      \code{CE.mmin} determines the initial size of the circulant
      matrix. If \code{CE.mmin=0} the minimal starting size is
      determined automatically according to the
      dimensions of the grid.
      If \code{CE.mmin>0} then the absolute starting size is given.
      If \code{CE.mmin<0} then the automatically determined
      matrix size is multiplied
      by \eqn{|\code{CE.mmin}|}; here \code{CE.mmin} must be smaller than -1; the
      value -1 takes over the minimal starting size.\cr
      Note: in any cases, the initial size might be increased according
      to \code{CE.useprimes}.\cr
      Default: \code{0} [init].} 
    \item{\code{CE.strategy}}{0 : if the circulant
      matrix has negative eigenvalues then the
      size in each direction is doubled; \cr 1 : the size is enhanced only in
      one direction, namely that one where the covariance function has the
      largest value at the end point of the grid --- note that
      the default value of \code{CE.trials} is probably too small
      in that case. \cr Default: \code{0}
      [init].}
    \item{\code{CE.maxmem}}{
      maximal total size of the circulant matrix. 
      The total amount of memory needed for the internal calculations
      is about 16 (=2 * sizeof(double))
      times as large as \code{CE.maxmem}
      if \command{\link{RFparameters}}\code{()$Storing=FALSE}
      and 32 (=4 * sizeof(double)) time as large if \code{Storing=TRUE}.

      Note that \code{CE.maxmem} can be used to control the automatic
      choice of the simulation algorithm. Namely, in case of huge
      circulant matrices, other simulation
      methods (TBM) are faster and might be preferred be the user.
      
      \cr Default: \code{20000000} [init].
    }
    \item{\code{CE.tolIm}}{
      If the modulus of the imaginary part is less than
      \code{CE.tolIm} then the eigenvalue is considered as real. 
      Default: \code{1E-3} [init].}
    \item{\code{CE.tolRe}}{
      Eigenvalues between \code{CE.tolRe} and 0 are considered as
      0 and set 0.  Default: \code{-1E-7} [init].}
    \item{\code{CE.trials}}{
      A larger circulant matrix is likely to make more eigenvalues
      non-negative. If at least one of the thresholds \code{CE.tolRe} and
      \code{CE.tolIm} are missed then the matrix size is doubled
      according to \code{CE.strategy},
      and the matrix is checked again.  This procedure is repeated
      up to \code{CE.trials-1} times.  If there are still negative
      eigenvalues, the simulation method fails if \code{CE.force=FALSE}. 
      Default: \code{3} [init].
    }
    \item{\code{CE.useprimes}}{Logical.  If \code{FALSE}
      the columns of the circulant matrix
      have length \eqn{2^k} for some \eqn{k}.  Otherwise the algorithm
      tries to find a nicely factorizable number close to the size of the
      given matrix.  Default: \code{TRUE}  [init].}
    \item{\code{CE.dependent}}{Logical.  If \code{FALSE}
      then independent random fields are created. If \code{TRUE}
      then at least 4 non-overlapping rectangles are taken out of the
      the expanded grid defined by the circulant matrix.
      These simulations are dependent.
      See below for an example.
      See \code{CE.trials} for some
      more information on the circulant matrix.
      Default: \code{FALSE} [init].}
  }
  
  
  \bold{Options for simulating with the local ce
    methods (cutoff, intrinsic)}\cr
  \describe{
    \item{\code{local.force}}{see CE.force above.
      Default: \code{FALSE} [init].
    }
    \item{\code{local.mmin}}{see CE.mmin above.
      
      Difference: if \code{local.mmin=0} the automatic determination
      of the initial size of the circulant matrix takes into
      account an expansion factor. This expansion factor is
      intended to make the circulant matrix positive definite
      and is either theoretically or numerically known, or guessed.

      If the usual strategy of circulant embedding
      (doubling the grid sizes) should be taken over then
      \code{local.mmin} must
      be set to \code{-1}.
      
      Default: \code{0} [init].}    
    \item{\code{local.strategy}}{see CE.strategy above.
      Default: \code{0} [init].}
    \item{\code{local.maxmem}}{see CE.maxmem above.
      Default: \code{20000000} [init].}
    \item{\code{local.tolIm}}{see CE.tolIm above.
      Default: \code{1E-7} [init].}
    \item{\code{local.tolRe}}{see CE.tolRe above.
      Default: \code{-1E-9} [init].}
    \item{\code{local.trials}}{see CE.trials above.  
      Default: \code{1} [init].}
    \item{\code{local.useprimes}}{see CE.useprimes above.
      Default: \code{TRUE} [init].}
    \item{\code{local.dependent}}{see CE.dependent above.
      Default: \code{FALSE} [init].}
  }
  %

  
  \bold{Options for simulating by simple matrix decomposition}\cr
  \describe{
    \item{\code{direct.checkprecision}}{Gaussian
      random vectors can be generated
      by means of the square root of the covariance matrix. 
      By default Cholesky decomposition is used.  If
      \code{direct.checkprecision==TRUE} then the precision is checked. 
      Default: \code{FALSE} [init].}
    \item{\code{direct.maxvariables}}{Decomposition of the covariance matrix.
      If the number of variables to generate is
      greater than \code{direct.maxvariables}, then any matrix decomposition
      method is rejected.  It is important that this option is set
      conveniently to avoid great losses of time during the automatic
      search of a simulation method (\code{method=NULL} in
      \command{\link{GaussRF}}).  
      Default: \code{1800} [init]}
    \item{\code{direct.method}}{Decomposition of the covariance matrix.
      If \code{direct.method==1}, Cholesky
      decomposition will not be attempted, but singular value
      decomposition
      used instead.
      Default: \code{0} [init].}
    \item{\code{direct.requiredprecision}}{Decomposition of the covariance
      matrix. 
      If \code{direct.checkprecision==TRUE} and
      the \code{direct.requiredprecision} is not reached then Cholesky
      decomposition fails, and singular value decomposition is used. 
      Default: \code{1e-11} [init].
    }
  }
  %

  
  \bold{Options for simulating nugget effects}\cr
  Simulating a nugget effect seems trivial. It gets complicated
  and best methods (including \code{direct} and \code{circulant
    embedding}!) fail if zonal anisotropies are considered,
  where sets of points have to be identified that belong to the
  same subspace of eigenvalue 0 of the anisotropy matrix.
  \describe{
    \item{nugget.tol}{
      points at a distance less than or equal to \code{nugget.tol}
      are considered as being identical. This strategy applies to
      the simulation method and the covariance function itself.
      Hence, the covariance function is only positive definite
      if \code{nugget.tol=0.0}. The value of \code{nugget.tol}
      should be of order \eqn{10^{-15}}{1e-15}.
      Default: \code{0.0} [init].
    }
  }


  \bold{Options for simulating with a turning bands method}\cr
  Currently, there are 3 variants of the turning bands method
  implemented:
  \describe{
    \item{\code{spectral}}{The spectral turning bands method
      is implemented for 2 (and 1)
      dimensions only.}
    \item{\code{TBM2}}{It is based on the two dimensional turning bands operator
      and is applicable for 1 and 2 dimensions. As an additional dimension
      the time dimension can be added.
      }
    \item{\code{TBM3}}{It is based on the three dimensional turning bands operator
      and is applicable for 1,2,3 dimensions. As an additional dimension
      the time dimension can be added.
      }
    }
   The following parameters are used. 
  \describe{
    \item{\code{spectral.grid}}{Logical.    The angle of the lines is random if
      \code{spectral.grid==FALSE}, 
      and \eqn{k\pi/}{k*pi/}\code{spectral.lines}
      for \eqn{k}{k} in \code{1:spectral.lines},
      otherwise.  Default: \code{TRUE} [do].}
    \item{\code{spectral.lines}}{Spectral turning bands.
      Number of lines used (in total for all additive components of the
      covariance function).  Default: \code{500} [do].}
    %\item{\code{TBM.method}}{Set at init time; setting ignored and stored
    %setting used if other parameters are identical to former parameters!
    %-- use DeleteKey, to make sure that the current setting is used.
    % [init]}
    \item{\code{TBM.method}}{character.
      The preferred method to simulate on the line for \code{TBM2} and
      \code{TBM3};
      currently either \code{'circulant embedding'} or \code{'direct'}.
      If \code{'direct'}
      then the method is overwritten if the number of points on the grid
      is larger than \code{direct.maxvariables}.
      If the circulant embedding method is used, then the \code{TBMCE}
      parameters below determine the behaviour of the circulant
      embedding algorithm.
      
      Default: \code{"circulant embedding"} [init].
    }
    \item{\code{TBM.center}}{Scalar or vector.
      If not \code{NA}, the \code{TBM.center} is used as the center of
      the turning bands for \code{TBM2} and \code{TBM3}.
      Otherwise the center is determined
      automatically such that the line length is minimal.
      See also \code{TBM.points} and the examples below.
      Default: \code{NA}  [init].
    }
    \item{\code{TBM.points}}{integer. If greater than 0,
      \code{TBM.points} gives the number of points simulated on the TBM
      line, hence 
      must be greater than the minimal number of points given by
      the size of the simulated field and the two paramters
      \code{TBMx.linesimufactor} and \code{TBMx.linesimustep}.
      If \code{TBM.points} is not positive the number of points is
      determined automatically.
      The use of \code{TBM.center} and \code{TBM.points} is highlighted
      in an example below.
      Default: \code{0}  [init].
      }
    \item{\code{TBM2.every}}{If \code{TBM2.every>0} then every
      \code{TBM2.every}th iteration is announced.
      Default: \code{0} [do].}
    \item{\code{TBM2.lines}}{
      Number of lines used.
      Default: \code{60} [do].}
    \item{\code{TBM2.linesimufactor}}{If \code{TBM2.linesimufactor} and
      \code{TBM2.linesimustep} must be non-negative; if
      \code{TBM2.linesimustep}
      is positive then \code{TBM2.linesimufactor} is ignored.
      If both
      parameters are naught then \code{TBM.points} is used (and must be
      positive).
      The grid on the line is \code{TBM2.linesimufactor}-times
      finer than the smallest distance. 
      See also \code{TBM2.linesimustep}. 
      Default: \code{2.0} [init].}
    \item{\code{TBM2.linesimustep}}{
      If \code{TBM2.linesimustep} is positive the grid on the line has lag
      \code{TBM2.linesimustep}. 
      See also \code{TBM2.linesimufactor}. 
      Default: \code{0.0} [init].}
    \item{\code{TBM2.num}}{
      Logical. If \code{TRUE} then the covariance function on the line
      is approximated numerically.
      If \code{FALSE} only those models are allowed
      that have an analytic representation on
      the line.
      Default: \code{TRUE} [init].}
    \item{\code{TBM3.every}}{If \code{TBM3.every>0} then every
      \code{TBM3.every}th iteration is announced.
      Default: \code{0} [do].}
    \item{\code{TBM3.lines}}{ 
      Number of lines used. 
      Default: \code{500} [do].}
    \item{\code{TBM3.linesimufactor}}{See \code{TBM2.linesimufactor} for
      the meaning.
      Default: \code{2.0} [init].}
    \item{\code{TBM3.linesimustep}}{See 
      \code{TBM2.linesimustep} for the meaning.   
      Default: \code{0.0} [init].}
    \item{\code{TBMCE.force}}{see \code{TBM.method} and \code{CE.force}
      Default: \code{FALSE} [init].}
    \item{\code{TBMCE.mmin}}{see \code{TBM.method} and 
      \code{CE.mmin}.  Default: \code{0} [init].}
    \item{\code{TBMCE.strategy}}{see \code{TBM.method} and
      \code{CE.strategy}.  Default: \code{0} [init].} 
    \item{\code{TBMCE.maxmem}}{see \code{TBM.method} and
      \code{CE.maxmem}.  Default: \code{10000000} [init].} 
    \item{\code{TBMCE.tolIm}}{see \code{TBM.method} and
      \code{CE.tolIm}.  Default: \code{1E-3} [init].} 
    \item{\code{TBMCE.tolRe}}{see \code{TBM.method} and
      \code{CE.tolRe}. Default: \code{-1E-7} [init].} 
    \item{\code{TBMCE.trials}}{see \code{TBM.method} and
      \code{CE.trials}.  Default: \code{3} [init].} 
    \item{\code{TBMCE.useprimes}}{see \code{TBM.method} and
      \code{CE.useprimes}.  Default: \code{TRUE} [init].} 
    \item{\code{TBMCE.dependent}}{see \code{TBM.method} and \code{CE.dependent}.
      Default: \code{FALSE} [init].}
  }
  %

  \bold{Options for simulating with Poisson point processes}\cr
  \describe{
    \item{\code{add.MPP.realisations}}{
      Number of superposed
      realisations (to approximate the normal distribution; total number
      for all (additive) components with same anisotropy);
      Default: \code{100} [do].}
    \item{\code{MPP.approxzero}}{Functions that
      do not have 
      compact support are set to zero outside the ball outside for which the
      function has absolute values less than \code{MPP.approxzero}. 
      Default: \code{0.001} [init].}
    \item{\code{MPP.radius}}{
      In order avoid edge effects, the simulation area is enlarged by
      a constant \eqn{r}{r} so that all marks have their
      (supposed) support in the ball with radius \eqn{r}{r} centred at
      the origin; see also \code{MPP.approxzero}.
      If \code{MPP.radius>0} the true radius \eqn{r}{r} is replaced by
      \code{MPP.radius}.
      Default: \code{0.0} [init].}
  }
%

  \bold{Options for simulating hyperplane tessellations}\cr
  \describe{
    \item{\code{hyper.superpos}}{integer.
      number of superposed hyperplane tessellations.
      Default: \code{300} [do].
    }
    \item{\code{hyper.maxlines}}{integer.
      Maximum number of allowed lines.
      Default: \code{1000} [init].
    }
    \item{\code{hyper.mar.distr}}{integer.
      code for the marginal distribution used in the
      simulation:
      \describe{
	\item{\code{0}}{uniform distribution}
	\item{\code{1}}{Frechet distribution with form parameter
	  \code{hyper.mar.param}}
	\item{\code{2}}{Bernoulli distribution (Binomial with \eqn{n=1}) with
	  parameter \code{hyper.mar.param}}
      }
      The parameter should not be changed yet.
      Default: \code{0} [do].
    }
    \item{\code{hyper.mar.param}}{Parameter used for the marginal
      distribution. The parameter should not be changed yet.
      Default: \code{0} [do].
    }
  }
  
  \bold{Options specific to simulating max-stable random fields}
  \describe{
    \item{\code{maxstable.maxGauss}}{Max-stable random fields.
      The simulation of the max-stable process based on random fields uses
      a stopping rule that necessarily needs a finite upper endpoint
      of the marginal distribution of the random field.
      In the case of extremal Gaussian random fields,
      see \command{\link{MaxStableRF}}, the upper endpoint is
      approximated by \code{maxstable.maxGauss}.
      Default: \code{3.0} [init].
    }
  }

  The following refers to the simulation of Gaussian random fields
  (\command{\link{InitGaussRF}}, \command{\link{GaussRF}}), but most
  parts also apply
  for the simulation of max-stable random fields
  (\command{\link{InitMaxStableRF}}, \command{\link{MaxStableRF}}).
  
  Some of the global parameters determine the basic settings of a
  simulation, e.g. \code{direct.method} (which chooses a square
  root of a positive definite matrix).  The values of
  such parameters are read by
  \command{\link{InitGaussRF}} and stored in an internal register.
  Changing
  such a parameter between calling \command{\link{InitGaussRF}} and calling
  \command{\link{DoSimulateRF}} or between subsequent calls of
  \command{\link{GaussRF}} will not have any effect.  These parameters have
  the flag "[init]".
  
  Parameters like \code{TBM2.lines} (which determines the number of
  i.i.d. processes to be simulated on the line)
  are only relevant when generating
  random numbers.  These parameters are read by \code{DoSimulateRF}
  (or by the second part of \command{\link{GaussRF}}), and
  are marked by "[do]".
  
  \code{Storing} has an influence on both, \command{\link{InitGaussRF}} and
  \command{\link{DoSimulateRF}}.  \command{\link{InitGaussRF}} may reserve
  more memory if \code{Storing==TRUE}.  \command{\link{DoSimulateRF}} will
  free the register
  if \code{Storing==FALSE}, whatever the value of \code{Storing} was
  when \command{\link{InitGaussRF}} was called. 

  The distinction between [init] and [do] is also relevant if 
  \command{\link{GaussRF}} is used and called a second time
  with the same parameters for the random field and if
  \code{RFparameters()$Storing==TRUE}.  
  Then \command{\link{GaussRF}} realises that the second call has the
  same random field parameters, and
  takes over the stored intermediate results (that have been calculated
  with the \code{RFparameters()} at that time).  To prevent the use of
  stored intermediate results or to take into account intermediate
  changes of RFparameters
  set \code{RFparameters(Storing==FALSE)} or use
  \command{\link{DeleteRegister}()} between calls of \code{GaussRF}.

  A programme that checks whether the parameters are well
  adapted to a specific simulation problem is given as an example of
  \command{\link{EmpiricalVariogram}()}.

  For further details on the implemented methods, see \link{RFMethods}.
}
\value{
  If any parameter has been given
  \code{RFparameters} returns an invisible list of
  the given parameters in full name.

  Otherwise the complete list of parameters is returned. Further the
  values of the following internal readonly variables are returned:
  \cr
  * \code{covmaxchar}:  max. name length for variogram/covariance
  models
  \cr
  * \code{covnr}: number of currently implemented
  variogram/covariance models
  \cr
  * \code{distrmaxchar}: max. name length for a distribution
  \cr
 * \code{distrnr}: number of currently implemented
    distributions
  \cr
 * \code{maxdim}: maximum number of dimensions for a
    random field
  \cr
  * \code{maxmodels}: maximum number of elementary models in
    definition of a complex covariance model
  \cr
 * \code{methodmaxchar}: max. name length for methods
  \cr
 * \code{methodnr}: number of currently implemented simulation methods
}
\references{
  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields.}
  Technical report ST 99-10, Dept. of Maths and Statistics,
  Lancaster University. 
}
\author{Martin Schlather, \email{schlath@hsu-hh.de}
  \url{http://www.unibw-hamburg.de/WWEB/math/schlath/schlather.html}}

\seealso{\command{\link{GaussRF}},
  \command{\link{GetPracticalRange}}, 
  \command{\link{MaxStableRF}},
  \code{\link{RandomFields}},
  and \command{\link{RFMethods}}.}

\examples{
 RFparameters(Storing=TRUE)
 str(RFparameters())


############################################################
##                                                        ##
##       use of TBM.points and TBM.center                 ##
##                                                        ##
############################################################
## The following example shows that the same realisation      
## can be obtained on different grid geometries (or point     
## configurations) using TBM3 (or TBM2)
                                             
% library(RandomFields, lib="~/TMP")

x1 <- seq(-150,150,1)
y1 <- seq(-15, 15, 1)
x2 <- seq(-50, 50, 1)
model <- "exponential"
param <- c(0, 1, 0, 10)

meth <- "TBM3"
###### simulation of a random field on long thing stripe
runif(1)
rs <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
DeleteAllRegisters()
z1 <- GaussRF(x1, y1, model=model, param=param, grid=TRUE, register=1,
              method=meth, TBM.center=0, Storing=TRUE)
% str(GetRegisterInfo(1))
get(getOption("device"))(height=1.55, width=12)
par(mar=c(2.2, 2.2, 0.1, 0.1))
image(x1, y1, z1, col=rainbow(100))
polygon(range(x2)[c(1,2,2,1)], range(y1)[c(1,1,2,2)], border="red", lwd=3)


###### definition of a random field on a square of shorter diagonal
assign(".Random.seed", rs, envir=.GlobalEnv)
z2 <- GaussRF(x2, x2, model=model, param=param, grid=TRUE, register=2,
              method=meth, TBM.center=0,
              TBM.points=length(GetRegisterInfo(1)$method[[1]]$mem$l))
% str(GetRegisterInfo(2))
get(getOption("device"))(height=4.3, width=4.3)
par(mar=c(2.2, 2.2, 0.1, 0.1))
image(x2, x2, z2, zlim=range(z1), col=rainbow(100))
polygon(range(x2)[c(1,2,2,1)], range(y1)[c(1,1,2,2)], border="red", lwd=3)


############################################################
##                                                        ##
##                   use of exactness                     ##
##                                                        ##  
############################################################
x <- seq(0, 1, 1/24)
model <- list(list(model="stable", var=1, scale=1, kappa=1.0),
              "+",
              list(model="gencauchy", var=1, scale=1, kappa=c(1, 2)),
              )

for (exactness in c(NA, FALSE, TRUE)) { 
  readline(paste("\n\nexactness: ", exactness, "; press return"))
  DeleteRegister()
  z <- GaussRF(x, x,  grid=TRUE, gridtriple=FALSE,
                model=model, exactness=exactness,
                stationary.only=NA, Print=4, n=1,
                TBM2.linesimustep=1, Storing=TRUE)
   str(GetRegisterInfo()$method)
 }


 #############################################################
 ## The following gives a tiny example on the advantage of  ##
 ## local.dependent=TRUE (and CE.dependent=TRUE) if in a    ##
 ## study most of the time is spent with simulating the     ##
 ## Gaussian random fields. Here, the covariance at a pair  ##
 ## of points is estimated.                                 ##
 #############################################################

# In the example below, local.dependent speeds up the simulation
# by about factor 27 at the price of an increased variance of
# factor 1.5

x <- seq(0, 1, len=10)
y <- seq(0, 1, len=10)
grid.size <- c(length(x), length(y))
model <- list(list(model="exp", var=1.1, aniso=c(2,1,0.5,1)))
CovarianceFct(matrix(c(1, -1), ncol=2), model=model) ## true value

# determine number of non-overlapping realisations on the torus
DeleteRegister()
InitGaussRF(x, y, model=model, grid=TRUE, method="cu")
blocks <- GetRegisterInfo()$method[[1]]$mem$new$method[[1]]$mem$simupos
(n <- prod(blocks))

m <- if (interactive()) 1000 else 2

# using local.dependent=TRUE...
c1 <- numeric(m)
DeleteRegister()
unix.time(
for (i in 1:m) {
  cat("\n", i)
  z <- GaussRF(x, y, model=model, grid=TRUE, method="cu", n=n,
               local.dependent=TRUE, pch="")
  c1[i] <- cov(z[1,length(y),], z[length(x), 1 , ])
}) # about 35 0.3 35 0  0
var(c1) # about 0.013

# using local.dependent=FALSE...
c2 <- numeric(m)
DeleteRegister()
unix.time(
for (i in 1:m) {
  cat("\n", i)
  z <- GaussRF(x, y, model=model, grid=TRUE, method="cu", n=n,
               local.dependent=FALSE, pch="")
  c2[i] <- cov(z[1,length(y),], z[length(x), 1 , ])
})  # about 950  3 950 0 0
var(c2) # about 0.0087

}

\keyword{spatial}
