\name{ShowModels}
\alias{ShowModels}
\title{Interactive Choice of Models and Parameters}
\description{
 \code{ShowModels} is an interactive plot of 
 models and one- or two-dimensional simulations;
 it allows for fitting variogram models by eye.
}
\usage{
ShowModels(covx=seq(0,4,0.02),x=NULL,y=NULL,legend=TRUE,
           register=0,fixed.rs=FALSE,method=NULL,empirical=NULL,...)
}
\arguments{
  \item{covx}{distances at which the covariance functions or the
    variograms are to be evaluated and plotted}
  \item{x}{if \code{NULL} no simulation; otherwise it gives the
    \eqn{x}{x} coordinates
   }
   \item{y}{if \code{NULL} at most one-dimensional simulations
     are given (depending on the value of \code{x}); otherwise it gives
     the \eqn{y}{y} coordinates of a two-dimensional grid}
%  \item{legend}{if \code{TRUE} then a legend is added to the
%    two-dimensional plot}
  \item{register}{register where intermediate results of the simulations
    are stored, see also \code{\link{GaussRF}}}
  \item{fixed.rs}{if \code{TRUE} then the same random seed is used for
    all simulations}
  \item{method}{simulation method, see \link{RFMethods}; if \code{NULL}
    then a suitable simulation method is chosen automatically}
  \item{empirical}{empirical variogram; list as returned by
  \code{\link{EmpiricalVariogram}}} 
  \item{...}{additional graphical options for the plot of the one- or
    two-dimensional simulations}
}
\details{
  The interactive plot consists of 4 parts:
  \itemize{
    \item top left: the graph of the covariance function or the
    variogram. In case \code{empirical} is given the empirical
    variogram is also plotted.  
    \item bottom left: a one- or two-dimensional simulation
    \item top right: list of implemented models; a specific
    model is chosen by the left mouse button.
    \item  bottom right: parameter list for the chosen model.
    This part is active if the name of the model is surrounded by
    asterixes. It is deactivated by clicking on the upper parts
    of the window. If the bottom right part is active then the
    numerical parameters of the model can be increased by clicking
    on the right hand side of the grey vertical line, and decreased
    by clicking on the left hand side. The farther away from the
    vertical line the larger the added or subtracted values.
    \cr
    The bottom lines give \emph{global} parameters, which determine
    \itemize{
      \item the \code{PracticalRange}, see
      \code{\link{RFparameters} for details}
      
      \item the way the model is shown in the
      upper left part, namely as covariance function or as variogram.
    }
  }
  The interactive plot is left by clicking any mouse button different
  from the left when the top right part is active.
}
\value{
  returns a list of the last parameters.
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
\seealso{\code{\link{CovarianceFct}}, \code{\link{GaussRF}},
  \code{\link{RFMethods}}, \code{\link{RandomFields}}.}
\examples{
# first example: one-dimensional simulations
close.screen(close.screen())
RFparameters(Print=1)$null
x <-  seq(1,10,0.1);
ShowModels(x=x)

# second example: two-dimensional simulations and
#                 empirical variogram
close.screen(close.screen())
RFparameters(Print=1)$null
dx <- runif(300,0,8)
dy <- runif(300,0,8)
dz <- GaussRF(x=dx, y=dy, grid=FALSE, model="gaus",
      param=c(1,2,1,2))
ev <- EmpiricalVariogram(x=dx, y=dy, data=dz, grid=FALSE,
      bin=(-1:20)/4)
x <-  seq(1,5,0.1);
ShowModels(x=x,y=x,empirical=ev)
}
\keyword{spatial}



