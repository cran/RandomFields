\name{ShowModels}
\alias{ShowModels}
\title{Interactive Choice of Models and Parameters}
\description{
 \code{ShowModels} is an interactive plot of 
 models and one- or two-dimensional simulations;
 it allows for fitting variogram models by eye.
}
\usage{
ShowModels(covx=ifelse(is.null(empirical),4,max(empirical$c)),
           x=NULL, y=NULL, fixed.rs=FALSE, method=NULL, empirical=NULL,
           model=NULL, param=NULL, all.param=NULL, PracticalRange=FALSE, 
           legend=TRUE, register=0,...)
}
\arguments{
  \item{covx}{if a sinlge value is given, it is the largest distance
    for which the covariance functions or the
    variograms are plotted; otherwise the models are plotted for the given
    values, and the origin.}
  \item{x}{if \code{NULL} simulations are not performed; otherwise it gives the
    \eqn{x}{x} coordinates of a grid.
   }
   \item{y}{if \code{NULL} at most one-dimensional simulations
     are performed (depending on the value of \code{x}); otherwise
     \code{y} gives
     the \eqn{y}{y} coordinates of a two-dimensional grid.}
  \item{fixed.rs}{if \code{TRUE} then the same random seed is used for
    all simulations.}
  \item{method}{simulation method, see \link{RFMethods}; if \code{NULL}
    then a suitable simulation method is chosen automatically.}
  \item{empirical}{empirical variogram; a list as returned by
    \code{\link{EmpiricalVariogram}}.}
  \item{model}{string; covariance model, see \code{\link{CovarianceFct}}, or
    type \code{\link{PrintModelList}()} to get all options. If given
    this model is shown at the beginning.
  }
  \item{param}{parameter vector:
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order; see \code{\link{CovarianceFct}} for more details.
    Only considered if \code{model} is given.
    If given \code{model} is initialized by \code{param}.}
  \item{all.param}{\code{all.param=c(mean, variance, nugget, scale)};
    the parameters must be given in this order;
    If \code{all.param} is given then the parameters of all covariance
    functions are set to the given values.  The values are overwritten
    for a specific model if \code{model} and \code{param} are given.
    Note that it is not possible to set the values of additional parameters
    of a parametrised class by means of \code{all.param}.}
  \item{PracticalRange}{logical.  If \code{TRUE} the range of the
    covariance functions \eqn{C(h)}{C(h)}
    is adjusted so that the covariance function is about 0.05 for
    \code{scale==1} and distance \eqn{h=1}{h=1}.}
  \item{legend}{if \code{TRUE} then a legend is added to the
    two-dimensional plot; note that \code{\link{legend}}
    has a bug up to R-1.2.3.}
  \item{register}{register where intermediate results of the simulations
    are stored, see also \code{\link{GaussRF}}.}
  \item{...}{additional graphics options for the plot of the one- or
    two-dimensional simulations, see \code{\link{plot}} and
    \code{\link{image}}.}
}
\details{
  The interactive plot consists of 4 parts:
  \itemize{
    \item top left: graph of the covariance function or the
    variogram. In case \code{empirical} is given the empirical
    variogram is also plotted.  
    \item bottom left: one- or two-dimensional simulation
    \item top right: list of implemented models; a specific
    model is chosen by the left mouse button.
    \item  bottom right: parameter list for the chosen model.
    This part is active if the name of the model is surrounded by
    asterixes. It is deactivated by clicking on the upper parts
    of the window. If the bottom right part is active then the
    numerical parameters of the model can be increased by clicking
    on the right hand side of the grey vertical line, and decreased
    by clicking on the left hand side. The farer away from the
    vertical line the larger the added or subtracted values.
    \cr
    The bottom lines give \emph{global} parameters, which determine
    \itemize{
      \item the \code{PracticalRange}, see
      \code{\link{RFparameters}} for details
      
      \item the way the model is shown in the
      upper left part, namely as covariance function or as variogram.
    }
  }
  The interactive plot is left by clicking any mouse button different
  from the left when the top right part is active.
}
\value{
   list of the last model and its parameters.
}
\author{Martin Schlather, \email{Martin.Schlather@uni-bayreuth.de}
  \url{http://www.geo.uni-bayreuth.de/~martin}}
\seealso{\code{\link{CovarianceFct}}, \code{\link{GaussRF}},
  \code{\link{RFMethods}}, \code{\link{RandomFields}}.}
\examples{
# first example: one-dimensional simulations
close.screen(close.screen())
RFparameters(Print=1)
x <-  seq(1,10,0.1);
ShowModels(x=x)

# second example: two-dimensional simulations and
#                 empirical variogram
close.screen(close.screen())
RFparameters(Print=1)
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



