\name{ShowModels}
\alias{ShowModels}
\title{Interactive Choice of Models and Parameters}
\description{
  \code{ShowModels} is an interactive plot for
  the selection of 
 models and their one- or two-dimensional simulations;
 it also allows for the fitting of variogram models by eye.

 CURRENTLY NOT AVAIABLE IN VERSION 2 OF RANDOMFIELDS.
 PLEASE USE VERSION 1.3.47 INSTEAD.
}
\usage{
ShowModels(x, y=NULL,
           covx=ifelse(is.null(empirical), diff(range(x))/5,
                       max(empirical$c)), 
           fixed.rs=TRUE, method=NULL, empirical=NULL,
           model=NULL, param=NULL, anisotropy = FALSE, all.param=NULL,  
           legends = TRUE, 
           register=0, Mean=NULL, erase=TRUE,
           x.fraction=0.60, cex.names=1, covx.default = 100,
           link.fct=NULL, Zlim=NULL,
           Col.rect="red", Col.bg="blue", Col.sep="grey",
           Col.line="red", Col.txt="black", Col.flash="red",
           Col.vario="blue", Col.main="black",
           Col.model=c("red", "black"),  vario.lty=c(1,2), 
           cex.leg = 0.7 * cex.names,  cex.eval = 0.8 * cex.names,
           update=TRUE, screen.new=TRUE, use.outer.RFparameters=FALSE,
           debug=FALSE, ...)
}
\arguments{
  \item{x}{if \code{NULL} simulations are not performed; otherwise it gives the
    \eqn{x}{x} coordinates of a grid as a sequence of increasing numbers
   }
   \item{y}{if \code{NULL} at most one-dimensional simulations
     are performed (depending on the value of \code{x}); otherwise
     \code{y} gives
     the \eqn{y}{y} coordinates of a two-dimensional grid (as a sequence
     of increasing numbers).}
  \item{covx}{if a single value is given, it is the largest distance
    for which the covariance functions or the
    variograms are plotted; otherwise the models are plotted for the given
    values, and the origin.}
  \item{fixed.rs}{if \code{TRUE} then the same random seed is used for
    all simulations until the user clicks on the formula, the title
    or the subtitles.}
  \item{method}{simulation method, see \link{RFMethods}; if \code{NULL}
    then a suitable simulation method is chosen automatically.}
  \item{empirical}{empirical variogram; a list as returned by
    \command{\link{EmpiricalVariogram}}. Also empirical variograms with
    a pair number of anisotropy directions may be passed. Then the first
    and the middle one are taken.
  }
  \item{model}{covariance model, see \command{\link{CovarianceFct}}, or
    type \command{\link{PrintModelList}()} to get all options. If given,
    this model is shown at the beginning.
    Additive or multiplicative models are not allowed.

    However, model can also be given by a simple list definition,
    see \command{\link{CovarianceFct}}. Then \code{param} must not be
    given.
    In this case also anisotropic models can be defined.
  }
  \item{param}{parameter vector:
    \code{param=c(mean, variance, nugget, scale,...)};
    the parameters must be given
    in this order; see \command{\link{CovarianceFct}} for more details.
    Only considered if \code{model} is given.
    If given, \code{model} is initialised by \code{param}.}
  \item{anisotropy}{logical. If \code{TRUE} then an isotropic model is
    considered. 
  }
  \item{all.param}{%% ## var, nugg, scale
    \code{all.param=c(mean, variance, nugget, scale)};
    the parameters must be given in this order;
    If \code{all.param} is given then the parameters of all covariance
    functions are set to the given values.  The values are overwritten
    for a specific model if \code{model} and \code{param} are given.
    Note that it is not possible to set the values of additional (form)
    parameters
    of a parametrised class by means of \code{all.param}.

    In case of an anisotropic model the anisotropy matrix is by default
    diagonal with both entries equal to \eqn{1/\code{all.param[4]}}.
  }
  \item{legends}{if \code{TRUE} then a legend is added to the
    two-dimensional plot.}
  \item{register}{register where intermediate results of the simulations
    are stored, see also \command{\link{GaussRF}}.}
  \item{Mean}{mean of the random field}
  \item{erase}{parameter of \command{\link{split.screen}}, which is called
    at the very beginning}
  \item{x.fraction}{the current screen is split into 2 x 2 screens.  
    The parameter \code{x.fraction} gives the size of the left screens in
    the x directions as part of 1. See also the Details.}
  \item{cex.names}{font size for model names}
  \item{covx.default}{if \code{length(cov.x)==1} then \eqn{[0,
      \code{cov.x}]} is covered by \code{covx.default} points of equal
    distance}
  \item{link.fct}{\code{NULL} or \code{function(values)} or "MaxStable".  
    Transformation of the  Gaussian random field.
    \cr\cr
    If \code{link.fct="MaxStable"} then
    max-stable random fields are simulated for the given covariance
    function and the extremal coefficient function is given (up the
    constant -1) instead of the variogram or the covariance function
  }
  \item{Zlim}{Vector of two elements or list of two vectors of two
    elements. Graphical limits for the Gaussian random process (and the
    transformed field).}
  \item{Col.rect}{colour for interactive plot; see
    \command{\link{eval.parameters}}.} 
  \item{Col.bg}{colour for interactive plot; see \command{\link{eval.parameters}}.}
  \item{Col.sep}{colour for interactive plot; see \command{\link{eval.parameters}}.}
  \item{Col.line}{colour for interactive plot; see
    \command{\link{eval.parameters}}.} 
  \item{Col.txt}{colour for interactive plot; see \command{\link{eval.parameters}}.}
  \item{Col.flash}{colour for the previously chosen model}
  \item{Col.vario}{colour for the empirical variogram plot}
  \item{Col.main}{colour for the title of the random field}
  \item{Col.model}{vector of two colours for plotting the variogram of
    the Gaussian random field and the transformed field}
  \item{vario.lty}{vector of two line types for primary and secondary
    axis of the variogram}
  \item{cex.leg}{font size used in the legends}
  \item{cex.eval}{font size used in the menue entries}
 \item{update}{logical.  If \code{TRUE} the plots are updated after
    each interactive change of the values.  Otherwise, the bottom
    'simulate' is added in the menu.}
  \item{screen.new}{logical.  If \code{FALSE} the screen is erased
    before a simulation and completely rebuild;  otherwise the screen is
    updated.  If \code{FALSE} flickering appears during the update of the
    current screen, otherwise it may happen during the reorganisation of
    any window (and may take quite a lot of time).
  }
  \item{use.outer.RFparameters}{logical.
    If \code{FALSE} the following parameters usually set by
    \command{\link{RFparameters}} are internally set
    \itemize{
      \item \code{PracticalRange=FALSE}
      \item \code{PrintLevel=1} if \code{debug=FALSE} and \code{5} otherwise.
      \item \code{maxstable.maxGauss=2}
      \item \code{CE.force=TRUE}
      \item \code{CE.trials=1}
      \item \code{CE.mmin=-4}
      \item \code{CE.useprimes=TRUE}
    }
  }
  \item{debug}{logical. If \code{TRUE} then internally the
    \code{RFparameter()$PrintLevel} is set to 5.
  }
  \item{...}{additional graphics options for the plot of the one- or
    two-dimensional simulations, see \command{\link{plot}} and
    \command{\link{image}}.}
}
\details{
  The interactive plot consists of 3 parts:
  \itemize{
    \item top left: graph of the covariance function or the
    variogram. In case \code{empirical} is given the empirical
    variogram is also plotted. If \code{link.fct} is given, then also
    the covariance function or the variogram is plotted. If the
    correlation model is for a non-stationary random field,
    the variogram for the transformed random field is not estimated in
    a primitive way -- this is indicated with a star in the legend
    \item bottom left: one- or two-dimensional simulation
    \item right:

      -- list of implemented models; a specific
      model is chosen by the left mouse button, or:
      
      -- menu for the parameters of the chosen model.
      The list includes the variance, a nugget effect, the mean and the
      scale or the anisotropy parameters. Further, some global parameters
      can be changed. They are the \code{PracticalRange}  (see
      \command{\link{RFparameters}} for details) and the angle of the
      main variogram direction (or NA, then it follows the angle of the
      anisotropy). Finally, the user can choose between the plot of the
      covariance and the corresponding variogram.
  }
  The interactive plot is left by clicking any mouse button different
  from the left when the top right part is active.
}
\value{
  list of the last model and its parameters.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de}}
\seealso{\command{\link{CovarianceFct}}, \command{\link{eval.parameters}},
  \command{\link{GaussRF}},
  \command{\link{RFMethods}}, \code{\link{RandomFields}}.}
\examples{
\dontrun{

%library(RandomFields, lib="~/TMP"); source("~/R/RF/RandomFields/R/ShowModels.R");source("~/R/RF/RandomFields/R/rf.R");source("~/R/RF/RandomFields/R/auxiliary.R")
# first example: one-dimensional simulations


# library(RandomFields)

options(locatorBell=FALSE)
x <-  seq(1,10,0.1)
ShowModels(x=x, covx=10, cov.def=100, type="l")


x <-  seq(1,10,0.1)
ShowModels(x=x, y=x, covx=10, cov.def=100)


# second example: two-dimensional simulations and
#                 empirical variogram
dx <- runif(300,0,8)
dy <- runif(300,0,8)
dz <- GaussRF(x=dx, y=dy, grid=FALSE, model="gaus",
      param=c(1,2,1,2))
ev <- EmpiricalVariogram(x=dx, y=dy, data=dz, grid=FALSE,
      bin=(-1:20)/4)
x <-  seq(1,5,0.1);
ShowModels(x=x, y=x, empirical=ev)

# third example: two-dimensional anistropic simulations and
#                link function
%  library(RandomFields, lib="~/TMP"); source("~/R/RF/RandomFields/R/ShowModels.R");RFparameters(Print=5,CE.trials=1);ShowModels(x=x, y=x, link=function(x) x^2,model=list(list(model="exp", var=1,aniso=c(1,0,0,5)))); 
x <- seq(1,10,0.1)
ShowModels(x=x, y=x, link=function(x) exp(x),
          model=list(list(model="spheric", var=1, aniso=c(1,0,0,5))))

x <-  seq(1,10,0.1)
ShowModels(x=x, link=function(x) exp(x),
           model=list(list(model="spheric",var=1, scale=1)))


%  library(RandomFields, lib="~/TMP");  source("~/R/RF/RandomFields/R/ShowModels.R");x <-  seq(1,10,0.05);ShowModels(x=x, link="MaxStable", fixed.rs=TRUE,model=list(list(model="gauss",var=1, scale=1)), type="l")
x <-  seq(1,10,0.1)
ShowModels(x=x, link="MaxStable", fixed.rs=TRUE,
           model=list(list(model="gauss",var=1, scale=1)), type="l")
} % dontrun

} % examples
\keyword{spatial}            



