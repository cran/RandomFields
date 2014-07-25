\name{RFgui}
\alias{RFgui}
\title{Graphical User Interface For Fitting Covariance Models And Variograms}
\description{
  This is a nice instructive graphical tool useful in particular for
  teaching classes 
}
\usage{
RFgui(data, x, y, same.algorithm = TRUE, ev, xcov, ycov,
      sim_only1dim=FALSE, wait = 0, ...)
}
\arguments{
 \item{data}{see \command{\link{RFempiricalvariogram}}. If \code{data}
   is given, the empirical variogram is shown}
 \item{x}{a \code{\link[base]{seq}}uence of  the
   locations of the
   simulated process; if not given,
   \code{x} is determined by \code{data} and if \code{data} is not
   given by default values}
 \item{y}{a \code{\link[base]{seq}}uence of numbers if a simulation
   on \eqn{R^d} is performed. Default is \code{y = x}; see \code{x} for details.}
 \item{same.algorithm}{Force the picture being simulated with the same
 algorithm so that the pictures are always directly comparable. The
 disadvantage is that some models are simulated only (very) approximatively.}
 \item{ev}{instead of the data, the emprirical variogram itself might be
 passed}
 \item{xcov}{ \code{\link[base]{seq}}uence of the
   locations where the covariance function is plotted}
 \item{ycov}{Only for anisotropic models. \code{\link[base]{seq}}uence of the
   locations where the covariance function is also plotted}
 \item{sim_only1dim}{Logical. The argument determines whether a process
   should be simulated on the line or on the plane}
 \item{wait}{integer. See details.
 }
 \item{...}{further options and control arguments for the simulation
 that are passed to and processed by \command{\link{RFoptions}}.}
}
\details{
 If \code{wait} is negative the xterm does not wait for the
 tkltk-window to be finished. Further the variable \code{RFgui.model}
 is created in the environment \code{.GlobalEnv} and contains the
 currently chosen variable in the gui. \command{\link{RFgui}} always return
 \code{NULL}. 

 If \code{wait} is non-negative the xterm waits for the
 tkltk-window to be finished. \command{\link{RFgui}} returns invisibly
 the last chosen model (or \code{NULL} if no model has been chosen).
 \command{\link{RFgui}} idles a lot when \code{wait=0}. It idles less for
 higher values by sleeping about \code{wait} microseconds. Of course
 the handling in the tkltk window get slower also. Reasonable values
 for \code{wait} are within \code{[0,1000]}.

 \cite{same.alg = TRUE} is equivalent to setting
 \code{circulant.trials=1},
 \code{circulant.simu_method = "RPcirculant"},
 \code{circulant.force=TRUE},
 \code{circulant.mmin=-2}.
 
}
 
\value{
 If \code{wait < 0} the the function returns \code{NULL}
 else it returns the last chosen \command{\link{RMmodel}}.

 If \code{wait < 0},
 a side effect effect of \command{\link{RFgui}}
 is the creation of the variable \code{RFgui.model}
 on \code{.GlobalEnv}.
}

\author{
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}

 \subsection{Author(s) of the code}{
 Daphne Boecker \email{d.boecker@gmx.de}

 Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
 }
}

\seealso{
 \code{\link{soil}} for a further example
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
RFgui()
\dontshow{FinalizeExample()}
}

\keyword{spatial}





