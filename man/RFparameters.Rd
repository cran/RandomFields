\name{RFparameters}
\alias{RFparameters}
\title{Control Parameters}
\description{
  \code{RFparameters} sets and returns control parameters for the simulation
  of random fields
}
\usage{
   RFparameters(...)
}
\arguments{
  \item{...}{arguments in \code{tag = value} form, or a list of tagged
    values.}
}
\details{
  The possible parameters are
  \describe{
    \item{\code{Storing}}{Logical.
      If \code{TRUE} then intermediate results are kept
    after each simulation; if several simulation are performed with the same
    model parameters then
    \code{Storing=TRUE} accelerates the simulations, but needs additional
    memory. \cr
    Note that in subsequent calls of \code{\link{GaussRF}}
    intermediate changes of RFparameters with flag "[init]"
    do not have any influence on the algorithm.
    Hence, for studying the effects for
    divers values of technical parameters like
    \code{CE.force}, \code{CE.mmin}, etc. the parameter \code{Storing} must be
    \code{FALSE}.  See also last paragraphs in the Details.
    Default: \code{TRUE} [init, do].}
  \item{\code{PrintLevel}}{If \code{PrintLevel}\eqn{\le0}{<=0}
    there is not any output on the screen.  The
    higher the number the more tracing information is given. 
    Default: 1 [init, do].\cr
    1 : error messages\cr
    2 : messages about partial failures of the algorithm
  }
  \item{\code{PracticalRange}}{logical or 2, 3, 11, 12, 13.
    If not \code{FALSE} the range of the covariance functions is
    adjusted so that cov(1) is about 0.05 (for \code{scale==1}).
    \itemize{
      \item \code{TRUE} : \code{PracticalRange} is applicable only if the value
      is known exactly.
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
  
   \item{\code{CE.force}}{Logical.  Circulant embedding does not work if a
     certain matrix has negative eigenvalues.  Sometimes it is convenient
     to replace all the negative eigenvalues by zero
     (\code{CE.force==TRUE}) after \code{CE.trials} number of trials. 
     Default: \code{FALSE} [init].
   }
   \item{\code{CE.mmin}}{Scalar or vector. 
     Circulant embedding usually uses the smallest matrix
     possible; by \code{CE.mmin} the minimum number of rows and columns
     of the matrix are given.  If \code{CE.mmin>=0} then the minimum
     absolute size is given.\cr
     If \code{CE.mmin<0} then the follow holds.
     If \code{CE.userfft==FALSE} then the matrix
     has \code{-CE.mmin} times the size of the size of the smallest
     matrix.  If \code{CE.userfft==TRUE} then the matrix
     has \code{max(-CE.mmin/2,1)} times the size of the size of the smallest
     matrix.  Default: \code{0} [init].}    
   \item{\code{CE.tolRe}}{Circulant embedding.
     Threshold above which eigenvalues are considered as
    non-negative.  Default: \code{-1E-5} [init].}
  \item{\code{CE.tolIm}}{Circulant embedding.
    If the modulus of the imaginary part is less than
    \code{CE.tolIm} then the eigenvalue is considered as real. 
    Default: \code{1E-3} [init].}
  \item{\code{CE.trials}}{Circulant embedding.
    A larger embedding matrix is likely to make more eigenvalues
    non-negative. If at least one of the thresholds \code{CE.tolRe} and
    \code{CE.tolIm} are missed then the matrix size is doubled,
    and the matrix is checked again.  This procedure is repeated
    up to \code{CE.trials-1} times.  If there are still negative
    eigenvalues, the simulation method fails if \code{CE.force==FALSE}. 
    Default: \code{3} [init].
  }
  \item{\code{CE.userfft}}{Logical.  If \code{FALSE}
    the columns of the circulant matrix
    have length \eqn{2^k} for some \eqn{k}.  Otherwise the algorithm
    tries to find a nicely factorizable number close to the size of the
    given matrix.  Default: \code{FALSE}  [init].}
  \item{\code{CE.strategy}}{0 : if matrix has negative eigenvalues then the
    size in each direction is doubled; \cr 1 : the size is enhanced only in
    one direction, namely that where the covariance function has the
    largest value at the end point of the grid - note that
    the default value of \code{CE.trials} is probably too small
    in that case. \cr Default: \code{0}
    [init].}
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
    conveniently if \code{method==NULL} in  \link{GaussRF}. 
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
  \item{\code{spectral.lines}}{Spectral turning bands.
    Number of lines used (in total for all additive components of the
    covariance function).  Default: \code{500} [do].}
  \item{\code{spectral.grid}}{Logical.  Spectral turning bands is implemented
    for 2 dimensions only.  The angle of the lines is random if
    \code{spectral.grid==FALSE}, 
    and \eqn{k\pi/}{k*pi/}\code{spectral.lines}
    for \eqn{k}{k} in \code{1:spectral.lines},
    otherwise.  Default: \code{TRUE} [do].}
  %\item{\code{TBM.method}}{Set at init time; setting ignored and stored
  %setting used if other parameters are identical to former parameters!
  %-- use DeleteKey, to make sure that the current setting is used.
  % [init]}
  \item{\code{TBMCE.force}}{Ordinary TBM methods.  At the moment only the
    circulant embedding method on the line is implemented; this
    parameter corresponds to \code{CE.force}. 
    Default: \code{FALSE} [init].}
  \item{\code{TBMCE.mmin}}{Ordinary TBM methods.  This parameter corresponds to
    \code{CE.mmin}.  Default: \code{0} [init].}
  \item{\code{TBMCE.tolRe}}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.tolRe}. Default: \code{-1E-5} [init].}
  \item{\code{TBMCE.tolIm}}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.tolIm}.  Default: \code{1E-3} [init].}
  \item{\code{TBMCE.trials}}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.trials}.  Default: \code{3} [init].}
  \item{\code{TBMCE.userfft}}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.userfft}.  Default: \code{true} [init].}
  \item{\code{TBMCE.strategy}}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.strategy}.  Default: \code{0} [init].}
  \item{\code{TBM2.lines}}{Ordinary 2-dimensional turning bands method. 
    Number of lines used. 
    Default: \code{60} [do].}
  \item{\code{TBM2.linesimufactor}}{Either \code{TBM2.linesimufactor} or
    \code{TBM2.linesimustep} must be greater than zero.  The parameter
    that is zero is ignored.  The grid on the line is
    \code{TBM2.linesimufactor}-times
    finer than the smallest distance. 
    See also \code{TBM2.linesimustep}. 
    Default: \code{2.0} [init].}
  \item{\code{TBM2.linesimustep}}{
    If \code{TBM2.linesimustep} is positive the grid on the line has lag
    \code{TBM2.linesimustep}. 
    See also \code{TBM2.linesimufactor}. 
    Default: \code{0.0} [init].}
  \item{\code{TBM2.every}}{If \code{TBM2.every>0} then every
    \code{TBM2.every}th iteration is announced.
    Default: \code{0} [do].}
  \item{\code{TBM3D2.lines}}{Ordinary 3-dimensional turning bands method,
     simulation of a \emph{2-dimensional} field. 
     Number of lines used. 
    Default: \code{500} [do].}
  \item{\code{TBM3D2.linesimufactor}}{Either \code{TBM3D2.linesimufactor} or
    \code{TBM2.linesimustep} must be greater than zero.  The parameter
    that is zero is ignored.  The grid on the line is
    \code{TBM3D2.linesimufactor}-times
    smaller than the smallest distance. See also \code{TBM3D2.linesimustep}. 
    Default: \code{2.0} [init].}
  \item{\code{TBM3D2.linesimustep}}{The grid on the line has lag
    \code{TBM3D2.linesimustep}.  See also \code{TBM3D2.linesimufactor}. 
    Default: \code{0.0} [init].}
  \item{\code{TBM3D2.every}}{If \code{TBM3D2.every>0} then every
    \code{TBM3D2.every}th iteration is announced.
    Default: \code{0} [do].}
  \item{\code{TBM3D3.lines}}{Ordinary 3-dimensional turning bands method,
     simulation of a \emph{3-dimensional field}. 
     Number of lines used. 
     Default: \code{500} [do].}
  \item{\code{TBM3D3.linesimufactor}}{Either \code{TBM3D3.linesimufactor} or
    \code{TBM2.linesimustep} must be greater than zero.  The parameter
    that is zero is ignored.  The grid on the line is
    \code{TBM3D3.linesimufactor}-times smaller than the smallest
    distance.  See also \code{TBM3D3.linesimustep}.
    Default: \code{2.0} [init].}
  \item{\code{TBM3D3.linesimustep}}{The grid on the line has lag
    \code{TBM3D3.linesimustep}.  See also \code{TBM3D3.linesimufactor}. 
    Default: \code{0.0} [init].}
  \item{\code{TBM3D3.every}}{If \code{TBM3D3.every>0} then every
    \code{TBM3D3.every}th iteration is announced.
    Default: \code{0} [do].}
  \item{\code{MPP.approxzero}}{Marked point processes. Functions that do not have
    compact support are set to zero outside the ball outside which the
    function has absolute values less than \code{MPP.approxzero}. 
    Default: \code{0.001} [init].}
  \item{\code{add.MPP.realisations}}{Random coins.
    Number of superposed
    realisations (to approximate the normal distribution; total number
    for all (additive) components with same anisotropy);
    Default: \code{100} [do].}
  \item{\code{MPP.radius}}{Marked point processes.
    In order avoid edge effects, the simulation area is enlarged by
    a constant \eqn{r}{r} so that all marks have their
    (supposed) support in the ball with radius \eqn{r}{r} centred at
    the origin; see also \code{MPP.approxzero}.
    If \code{MPP.radius>0} the true radius \eqn{r}{r} is replaced by
    \code{MPP.radius}.
    Default: \code{0.0} [init].}
  \item{\code{maxstable.maxGauss}}{Max-stable random fields.
    The simulation of the max-stable process based on random fields uses
    a stopping rule that necessarily needs a finite upper endpoint
    of the marginal distribution of the random field.
    In the case of extremal Gaussian random fields,
    see \code{\link{MaxStableRF}}, the upper endpoint is
    approximated by \code{maxstable.maxGauss}.
    Default: \code{3.0} [init].
  }
  \item{\code{pch}}{Character. The character is printed after each
    performed simulation if more than one simulation is performed at
    once.  Default: \code{'#'} [do]. 
  }
}



  The following refers to the simulation of Gaussian random fields
  (\code{\link{InitGaussRF}}, \code{\link{GaussRF}}), but most
  parts also apply
  for the simulation of max-stable random fields
  (\code{\link{InitMaxStableRF}}, \code{\link{MaxStableRF}}).
  
  Some of the global parameters determine the basic settings of a
  simulation, e.g. \code{direct.method} (which chooses a square
  root of a positive definite matrix).  The values of
  such parameters are read by
  \code{\link{InitGaussRF}} and stored in an internal register.
  Changing
  such a parameter between calling \code{\link{InitGaussRF}} and calling
  \code{\link{DoSimulateRF}} or between subsequent calls of
  \code{\link{GaussRF}} will not have any effect.  These parameters have
  the flag "[init]".
  
  Parameters like \code{TBM2.lines} (which determines the number of
  i.i.d. processes to be simulated on the line)
  are only relevant when generating
  random numbers.  These parameters are read by \code{DoSimulateRF}
  (or by the second part of \code{\link{GaussRF}}), and
  are marked by "[do]".
  
  \code{Storing} has an influence on both, \code{\link{InitGaussRF}} and
  \code{\link{DoSimulateRF}}.  \code{\link{InitGaussRF}} may reserve
  more memory if \code{Storing==TRUE}.  \code{\link{DoSimulateRF}} will
  free the register
  if \code{Storing==FALSE}, whatever the value of \code{Storing} was
  when \code{\link{InitGaussRF}} was called. 

  The distinction between [init] and [do] is also relevant if 
  \code{\link{GaussRF}} is used and called a second time
  with the same parameters for the random field and if
  \code{RFparameters()$Storing==TRUE}.  
  Then \code{\link{GaussRF}} realises that the second call has the
  same random field parameters, and
  takes over the stored intermediate results (that have been calculated
  with the \code{RFparameters()} at that time).  To prevent the use of
  stored intermediate results or to take into account intermediate
  changes of RFparameters
  set \code{RFparameters(Storing==FALSE)} or use
  \code{\link{DeleteRegister}()} between calls of \code{GaussRF}.

  A programme that checks whether the parameters are well
  adapted to a specific simulation problem is given as an example of
  \code{\link{EmpiricalVariogram}()}.

  For further details on the implemented methods, see \link{RFMethods}.
 }
\value{
  If any parameter has been given
  \code{RFparameters} returns an invisible list of
  the given parameters in full name.

  Otherwise the complete list of parameters is returned. Further the
  values of the following internal readonly variables are returned
  
  \item{\code{covmaxchar}}{max. name length for variogram/covariance
  models}
  \item{\code{methodmaxchar}}{max. name length for methods}
  \item{\code{distrmaxchar}}{max. name length for a distribution}
  \item{\code{covnr}}{number of currently implemented
    variogram/covariance models}
  \item{\code{methodnr}}{number of currently implemented
    variogram/covariance models}
  \item{\code{distrnr}}{number of currently implemented
    distributions} 
  \item{\code{maxdim}}{maximum number of dimensions for a
    random field}
  \item{\code{maxmodels}}{maximum number of models, i.e. the possible
      register numbers in \command{GaussRF} for example, are
      \eqn{1,\ldots,\code{maxmodels}}{1,...,\code{maxmodels}}.
    }
}
\references{
  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields.}
  Technical report ST 99-10, Dept. of Maths and Statistics,
  Lancaster University. 
}
\author{Martin Schlather, \email{martin.schlather@cu.lu}
  \url{http://www.cu.lu/~schlathe}}

\seealso{\code{\link{GaussRF}},
  \code{\link{GetPracticalRange}}, 
  \code{\link{MaxStableRF}},
  \code{\link{RandomFields}},
  and \code{\link{RFMethods}}.}

\examples{
 RFparameters(Storing=TRUE)
 str(RFparameters())
}
\keyword{spatial}
