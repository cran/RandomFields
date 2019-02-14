% M-x auto-fill-mode

\newcommand{\Examples}{\examples{\dontshow{StartExample()}RFoptions(seed=0)#1\dontshow{FinalizeExample()}}}

\newcommand{\me}{\author{\martin}}

\newcommand{\argModel}{ object of class  \code{\link[=RMmodel-class]{RMmodel}}, \command{\link{RFformula}} or  \command{\link[stats]{formula}}; best is to consider the examples below, first. \cr The argument \code{params} is a list that specifies free parameters in a formula description, see \link{RMformula}.}

\newcommand{\argX}{vector of x coordinates, or object of class \code{\link[sp:GridTopology-class]{GridTopology}} or \code{\link[raster]{raster}}; for more options see \link{RFsimulateAdvanced}.}

\newcommand{\argYz}{optional vectors of y (z) coordinates, which should not be given if \code{x} is a matrix.}

\newcommand{\argT}{optional vector of time coordinates, \code{T} must always be an equidistant vector. Instead of \code{T=seq(from=From, by=By, len=Len)}, one may also write \code{T=c(From, By, Len)}.}

\newcommand{\argGrid}{logical; the function finds itself the correct value in nearly all cases, so that usually \code{grid} need not be given. See also \link{RFsimulateAdvanced}.}

\newcommand{\argData}{matrix, data.frame or object of class \command{\link[=RFsp-class]{RFsp}}; \cr If a matrix is given the ordering of the colums is the following: space, time, multivariate, repetitions, i.e. the index for the space runs the fastest and that for repetitions the slowest.}

\newcommand{\argBin}{a vector giving the borders of the bins; If not specified an array describing the empirical (pseudo-)(cross-) covariance function  in every direction is returned.}

\newcommand{\argPhi}{an integer defining the number of sectors one half of the X/Y plane shall be divided into. If not specified, either an array is returned (if bin missing) or isotropy is assumed (if bin specified).}

\newcommand{\argTheta}{an integer defining the number of sectors one half of the X/Z plane shall be divided into. Use only for dimension \eqn{d=3}{d=3} if phi is already specified.}

\newcommand{\argDeltaT}{vector of length 2, specifying the temporal bins. The internal bin vector becomes \code{seq(from=0, to=deltaT[1], by=deltaT[2])}}

\newcommand{\argDistances}{another alternative for the argument \code{x} to pass the (relative) coordinates, see \link{RFsimulateAdvanced}.}

\newcommand{\argVdim}{the number of variables of a multivariate data set. If not given and \code{data} is an \code{RFsp} object created by \pkg{RandomFields}, the information there is  taken from there. Otherwise \code{vdim} is assumed to be one. NOTE: still the argument \code{vdim} is an experimental stage.} % to do

\newcommand{\argDots}{for advanced use: further options and control arguments for the simulation that are passed to and processed by \command{\link{RFoptions}}. If \code{params} is given, then \code{...} may include also the variables used in \code{params}.}

\newcommand{\argLower}{list or vector. Lower bounds for the parameters. If \code{lower} is a vector, \code{lower} has to be a vector as well and its length must equal the number of parameters to be estimated. The order of \code{lower} has to be maintained. A component being \code{NA} means that no manual lower bound for the corresponding parameter is set. \cr  If \code{lower} is a list, \code{lower} has to be of (exactly) the same structure of the model.}

\newcommand{\argUpper}{list or vector. Upper bounds for the parameters. See \code{lower}.}

\newcommand{\argUsersguess}{ User's guess of the parameters. All the parameters must be given using the same rules as for \code{lower} (except that no NA's should be contained).}

\newcommand{\argTransform}{obsolete for users; use \code{param} instead. \code{transform=list()} will return structural information to set up the correct function.}

\newcommand{\argOptimcontrol}{control list for \command{\link[stats]{optim}}, which uses \sQuote{L-BFGS-B}. However \code{parscale} may not be given.}

\newcommand{\argGiven}{optional, matrix or list.  If \code{given} matrix then the coordinates  can be given separately, namely by \code{given} where, in each row, a single location is given.\cr  If \code{given} is a  list, it may consist of \code{x}, \code{y}, \code{z}, \code{T},  \code{grid}.\cr  If \code{given} is provided, \code{data} must be a matrix or an array containing the data only. }

\newcommand{\argDataGiven}{ If \code{given} is not given and \code{data} is a matrix or \code{data} is  a data.frame, \pkg{RandomFields} tries to identify  where the data and the coordinates are, e.g. by names in formulae or by fixed names,  see \link{Coordinate systems}. See also \link{RFsimulateAdvanced}.  If all fails, the first columns are interpreted as coordinate vectors, and the last column(s) as (multiple) measurement(s) of the field. Notes that also lists of data can be passed.}

\newcommand{\argErrmodel}{In case of (assumed) error-free measurements (which is mostly the case in geostatistics) the argument \code{err.model} is not given. In case of measurement errors we have  \code{err.model=RMnugget(var=var)}.\cr \code{err.param} plays the same role as \code{params} for \code{model}.}

\newcommand{\argFitmethods}{ Main methods to be used for estimating. If several methods are given, estimation will be performed with each method and the results reported.}

\newcommand{\argFitsubmethods}{variants of the least squares fit of the variogram. variants of the maximum likelihood fit of the covariance function.}


\newcommand{\argOrigin}{character; one of \code{"original"}, \code{"MLE conform"}, \code{"all"}. This argument determines the parameters that are returned.}

\newcommand{\argRaw}{ logical. If \code{FALSE} then the data are interpolated. This approach is always save, but might be slow. If \code{TRUE} then the data may be accessed when covariance  matrices are calculated. No rescaling or anisotropy definition  is allowed in combination with the model. The use is dangerous, but fast.\cr  Default: \code{FALSE}.
}

\newcommand{\argNorm}{optional model that gives the norm between locations}

\newcommand{\GEV}{The argument \code{xi} is always a number, i.e. \eqn{\xi} is constant in space. In contrast, \eqn{\mu} and \eqn{s} might be constant numerical values or (in future!) be given by an \code{\link{RMmodel}}, in particular by an \code{\link{RMtrend}} model.  \cr For \eqn{xi=0}, the default values of \eqn{mu} and \eqn{s} are \eqn{0} and \eqn{1}, respectively. For \eqn{xi\not=0}, the default values of \eqn{mu} and \eqn{s} are \eqn{1} and \eqn{|\xi|}, respectively, so that it defaults to the standard Frechet case if \eqn{\xi > 0}.}
 


%  \newcommand{\arg}{ }
