\name{RFparameters}
\alias{RFparameters}
\alias{RFparameters.default}
\title{Control Parameters}
\description{
  \code{RFparameters} sets and returns control parameters for the simulation
  of random fields
}
\usage{
  RFparameters(...)

  RFparameters.default(Storing=storing, PrintLevel=printlevel,
                       PracticalRange=practicalrange, CE.force=ce.force,
                       CE.mmin=ce.mmin, CE.tolRe=ce.tolRe,
                       CE.tolIm=ce.tolIm, CE.trials=ce.trials,
                       direct.checkprecision=directcheckprecision,
                       direct.maxvariables=directmaxvariables,
                       direct.method=directmethod,
                       direct.requiredprecision=directrequiredprecision,
                       spectral.lines=spectrallines,
                       spectral.grid=spectralgrid,
                       TBMCE.force=tbmceforce, TBMCE.mmin=tbmcemmin,
                       TBMCE.tolRe=tbmcetolre, TBMCE.tolIm=tbmcetolim,
                       TBMCE.trials=tbmcetrials, TBM2.lines=tbm2lines,
                       TBM2.linesimufactor=tbm2linesimufactor,
                       TBM2.linesimustep=tbm2linesimustep,
                       TBM3D2.lines=tbm3D2lines,
                       TBM3D2.linesimufactor=tbm3D2linesimufactor,
                       TBM3D2.linesimustep=tbm3D2linesimustep,
                       TBM3D3.lines=tbm3D3lines,
                       TBM3D3.linesimufactor=tbm3D3linesimufactor,
                       TBM3D3.linesimustep=tbm3D3linesimustep,
                       MPP.approxzero=mppapproxzero,
                       add.MPP.realisations=addmpprealisations,
                       MPP.radius=mppradius,
                       maxstable.maxGauss=maxstablemaxGauss, pch=pchx)
}
\arguments{
  \item{...}{arguments as given in \code{RFparameters.default} and
    listed in the following.}
  \item{Storing}{logical.  If \code{TRUE} then intermediate results are kept
    after each simulation; if several simulation are made with the same
    parameters (e.g., by \code{n}\eqn{>1}{>1} in \link{GaussRF} or
    several calls of \link{GaussRF}) then
    \code{Storing=TRUE} accelerates the simulations, but needs additional
    memory.  Default: \code{TRUE} [init, do].}
  \item{PrintLevel}{If \code{PrintLevel}\eqn{\le0}{<=0}
    there is not any output on the screen.  The
    higher the number the more tracing information. 
    Default: 1 [init, do].\cr
    1 : messages about errors occurred\cr
    2 : messages about partial failures of the algorithm
  }
  \item{PracticalRange}{The range of the covariance functions can be
    adjusted so that cov(1) is about 0.05 (for \code{scale==1}).
    Default: \code{FALSE} [init].}
  
   \item{CE.force}{logical.  Circulant embedding does not work if a
     certain matrix has negative eigenvalues.  Sometimes it is convenient
     to replace all the negative eigenvalues by zero
     (\code{CE.force==TRUE}) after \code{CE.trials} number of trials. 
     Default: \code{FALSE} [init].
   }
   \item{CE.mmin}{Circulant embedding usually uses the smallest matrix
     possible; by \code{CE.mmin} the minimum number of rows and columns
     of the matrix are given.  Default:  \code{0} [init].}
   \item{CE.tolRe}{Circulant embedding.
     Threshold above which eigenvalues are considered as
    non-negative.  Default: \code{-1E-5} [init].}
  \item{CE.tolIm}{Circulant embedding.
    If the modulus of the imaginary part is less than
    \code{CE.tolIm} then the eigenvalue is considered as real. 
    Default: \code{1E-3} [init].}
  \item{CE.trials}{Circulant embedding.
    A larger embedding matrix is likely to make more eigenvalues
    non-negative. If at least one of the thresholds \code{CE.tolRe} and
    \code{CE.tolIm} are missed then the matrix size is doubled,
    and the matrix is checked again.  This procedure is repeated
    up to \code{CE.trials-1} times.  If there are still negative
    eigenvalues, the simulation method fails if \code{CE.force==FALSE}. 
    Default: \code{3} [init].
    }
  
  \item{direct.checkprecision}{Gaussian
    random vectors can be generated
    by means of the square root of the covariance matrix. 
    By default Cholesky decomposition is used.  If
    \code{direct.checkprecision==TRUE} then the precision is checked. 
    Default: \code{FALSE} [init].}
  \item{direct.maxvariables}{Decomposition of the covariance matrix.
    If the number of variables to generate is
    greater than \code{direct.maxvariables}, then any matrix decomposition
    method is rejected.  It is important that this option is set
    conveniently if \code{method==NULL} in  \link{GaussRF}. 
    Default: \code{1000} [init]}
  \item{direct.method}{Decomposition of the covariance matrix.
    If \code{direct.method==1}, Cholesky
    decomposition will not be attempted, but singular value
    decomposition
    used instead.
    Default: \code{0} [init].}
  \item{direct.requiredprecision}{Decomposition of the covariance
    matrix. 
    If \code{direct.checkprecision==TRUE} and
    the \code{direct.requiredprecision} is not reached then Cholesky
    decomposition fails, and singular value decomposition is used. 
    Default: \code{1e-11} [init].
  }
  \item{spectral.lines}{Spectral turning bands.
    Number of lines used.  Default: \code{500} [do].}
  \item{spectral.grid}{Logical.  Spectral turning bands is implemented
    for 2 dimensions only.  The angle of the lines is random if
    \code{spectral.grid==FALSE}, 
    and \eqn{k\pi/}{k*pi/}\code{spectral.lines}
    for \eqn{k}{k} in \code{1:spectral.lines},
    otherwise.  Default: \code{TRUE} [do].}

  %\item{TBM.method}{Set at init time; setting ignored and stored
  %setting used if other parameters are identical to former parameters!
  %-- use DeleteKey, to make sure that the current setting is used.
  % [init]}
  \item{TBMCE.force}{Ordinary TBM methods.  At the moment only the
    circulant embedding method on the line is implemented; this
    parameter corresponds to \code{CE.force}. 
    Default: \code{FALSE} [init].}
  \item{TBMCE.mmin}{Ordinary TBM methods.  This parameter corresponds to
    \code{CE.mmin}.  Default: \code{0} [init].}
  \item{TBMCE.tolRe}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.tolRe}. Default: \code{-1E-5} [init].}
  \item{TBMCE.tolIm}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.tolIm}.  Default: \code{1E-3} [init].}
  \item{TBMCE.trials}{Ordinary TBM methods.  This parameter corresponds
    to \code{CE.trials}.  Default: \code{3} [init].}
  
  \item{TBM2.lines}{Ordinary 2-dimensional turning bands method. 
    Number of lines used. 
    Default: \code{60} [do].}
  \item{TBM2.linesimufactor}{Either \code{TBM2.linesimufactor} or
    \code{TBM2.linesimustep} must be greater than zero.  The parameter
    that is zero is ignored.  The grid on the line is
    \code{TBM2.linesimufactor}-times
    smaller than the smallest distance. 
    See also \code{TBM2.lines}. 
    Default: \code{2.0} [init].}
  \item{TBM2.linesimustep}{The grid on the line has lag
    \code{TBM2.linesimustep}.  See also \code{TBM2.linesimufactor}. 
    Default: \code{0.0} [init].}
  \item{TBM3D2.lines}{Ordinary 3-dimensional turning bands method,
     simulation of a \emph{2-dimensional} field. 
     Number of lines used. 
    Default: \code{500} [do].}
  \item{TBM3D2.linesimufactor}{Either \code{TBM3D2.linesimufactor} or
    \code{TBM2.linesimustep} must be greater than zero.  The parameter
    that is zero is ignored.  The grid on the line is
    \code{TBM3D2.linesimufactor}-times
    smaller than the smallest distance. See also \code{TBM3D2.lines}. 
    Default: \code{2.0} [init].}
  \item{TBM3D2.linesimustep}{The grid on the line has lag
    \code{TBM3D2.linesimustep}.  See also \code{TBM3D2.linesimufactor}. 
    Default: \code{0.0} [init].}
  \item{TBM3D3.lines}{Ordinary 3-dimensional turning bands method,
     simulation of a \emph{3-dimensional field}. 
     Number of lines used. 
     Default: \code{500} [do].}
  \item{TBM3D3.linesimufactor}{Either \code{TBM3D3.linesimufactor} or
    \code{TBM2.linesimustep} must be greater than zero.  The parameter
    that is zero is ignored.  The grid on the line is
    \code{TBM3D3.linesimufactor}-times smaller than the smallest
    distance.  See also \code{TBM3D3.lines}.
    Default: \code{2.0} [init].}
  \item{TBM3D3.linesimustep}{The grid on the line has lag
    \code{TBM3D3.linesimustep}.  See also \code{TBM3D3.linesimufactor}. 
    Default: \code{0.0} [init].}
  \item{MPP.approxzero}{Marked point processes. Functions that do not have
    compact support are set to zero outside the ball outside which the
    function has absolute values less than \code{MPP.approxzero}. 
    Default: \code{0.001} [init].}
  \item{add.MPP.realisations}{Random coins.
    Number of superposed
    realisations (to approximate the normal distribution).
    Default: \code{100} [do].}
  \item{MPP.radius}{Marked point processes.
    In order avoid edge effects, the simulation area is enlarged by
    a constant \eqn{r}{r} so that all marks have their
    (supposed) support in the ball with radius \eqn{r}{r} centred at
    the origin; see also \code{MPP.approxzero}.
    If \code{MPP.radius>0} the true radius \eqn{r}{r} is replaced by
    \code{MPP.radius}.
    Default: \code{0.0} [init].}
  \item{maxstable.maxGauss}{Max-stable random fields.
    The simulation of the max-stable process based on random fields uses
    a stopping rule that necessarily needs a finite upper endpoint
    of the marginal distribution of the random field.
    In the case of extremal Gaussian random fields,
    see \code{\link{MaxStableRF}}, the upper endpoint is
    approximated by \code{maxstable.maxGauss}.
    Default: \code{3.0} [init].
  }
  \item{pch}{character. The character is printed after each
    performed simulation if more than one simulation is performed at
    once.  Default: \code{"#"} [do]. 
  }
}
\details{
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
  \code{\link{DoSimulateRF}} will not have any effect.  These parameters have
  the flag "[init]".
  
  Parameters like \code{TBM2.lines} (which determines the number of
  i.i.d. proceses to be simulated on the line)
  are only relevant when generating
  random numbers.  These parameters are read by \code{DoSimulateRF}, and
  are marked by "[do]".
  
  \code{Storing} has an influence on both, \code{\link{InitGaussRF}} and
  \code{\link{DoSimulateRF}}.  \code{\link{InitGaussRF}} may reserve
  more memory if \code{Storing==TRUE}.  \code{\link{DoSimulateRF}} will
  free the register
  if \code{Storing==FALSE}, whatever the value of \code{Storing} was
  when \code{\link{InitGaussRF}} was called. 

  The distinction between [init] and [do] is relevant even if 
  \code{\link{GaussRF}} is used, but called a second time
  with the same parameters for the random field and if
  \code{RFparameters()$Storing==TRUE}.  
  Then \code{\link{GaussRF}} realises that the second call has the
  same parameters, and
  takes over the stored intermediate results (that have been calculated
  with the \code{RFparameters()} at that time).  To prevent this
  put \code{RFparameters(Storing==FALSE)} or use
  \code{\link{DeleteRegister}()}.

  A programme that checks whether the parameters are well
  adapted to a specific simulation problem is given as an example of
  \code{\link{EmpiricalVariogram}()}.

  For further details on the implemented methods, see \link{RFMethods}.
 }
\value{
  returns \code{NULL} if any parameter has been given,
  and the list of all parameter values otherwise.
}
\references{
  Schlather, M. (1999) \emph{An introduction to positive definite
    functions and to unconditional simulation of random fields.}
  Technical report ST 99-10, Dept. of Maths and Statistics,
  Lancaster University. 
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}

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
