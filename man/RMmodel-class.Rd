\name{RMmodel-class}
\docType{class}
\alias{RMmodel-class}
\alias{show,RMmodel-method}
\alias{print.RMmodel}
\alias{plot,RMmodel,missing-method}
\alias{points.RMmodel}
\alias{lines.RMmodel}
\alias{str.RMmodel}
\alias{[,RMmodel-method}
\alias{[,RMmodel,ANY,ANY,ANY-method}
\alias{[<-,RMmodel-method}
\alias{[<-,RMmodel,ANY,ANY,ANY-method}
\alias{+,RMmodel,RMmodel-method}
\alias{+,numeric,RMmodel-method}
\alias{+,logical,RMmodel-method}
\alias{+,RMmodel,numeric-method}
\alias{+,RMmodel,logical-method}
\alias{*,RMmodel,RMmodel-method}
\alias{*,numeric,RMmodel-method}
\alias{*,logical,RMmodel-method}
\alias{*,RMmodel,logical-method}
\alias{*,RMmodel,numeric-method}
\alias{-,RMmodel,RMmodel-method}
\alias{-,numeric,RMmodel-method}
\alias{-,logical,RMmodel-method}
\alias{-,RMmodel,logical-method}
\alias{-,RMmodel,numeric-method}
\alias{/,RMmodel,RMmodel-method}
\alias{/,numeric,RMmodel-method}
\alias{/,logical,RMmodel-method}
\alias{/,RMmodel,numeric-method}
\alias{/,RMmodel,logical-method}
\alias{^,RMmodel,RMmodel-method}
\alias{^,numeric,RMmodel-method}
\alias{^,logical,RMmodel-method}
\alias{^,RMmodel,numeric-method}
\alias{^,RMmodel,logical-method}
\alias{c,RMmodel-method}
%\alias{plot,RMmodel-method}
%\alias{points.RMmodel}
%\alias{lines.RMmodel}

\title{Class \code{RMmodel}}
\description{ Class for \pkg{RandomField}'s representation of explicit
 covariance models 
}


\usage{
\S4method{plot}{RMmodel,missing}(x, y, dim=1, n.points=200,
 fct.type=NULL, MARGIN, fixed.MARGIN, maxchar=15, ...)

\S3method{points}{RMmodel}(x, y, n.points=200, fct.type=NULL, ...)

\S3method{lines}{RMmodel}(x, y, n.points=200, fct.type=NULL, ...)
}



\arguments{
  \item{x}{object of class \code{\link[=RFsp-class]{RFsp}} or
    \command{\link[=RFempVariog-class]{RFempVario}} or
    \command{\link[=RFfit-class]{RFfit}} or
    \command{\link[=RMmodel-class]{RMmodel}}; in the latter case, \code{x} can
    be any sophisticated model but it must be either stationary or a
    variogram model} 
  \item{y}{ignored in most methods}
  \item{MARGIN}{vector of two; two integer values giving the coordinate
    dimensions 
    w.r.t. which the field or the covariance model
    is to be plotted; in all other directions,
    the first index is taken}
%  \item{MARGIN.slices}{integer value; if \eqn{[space-time-dimension>2]},
%    \code{MARGIN.slices} can specify a third dimension w.r.t. which a
%    sequence of slices is plotted. Currently only works for grids.} 
  \item{fixed.MARGIN}{only for \code{class(x)=="RMmodel"} and if
    \code{dim > 2}; a vector of length \code{dim}-2 with distance values
    for the coordinates that are not displayed}
  \item{maxchar}{integer. Maximum number of characters to print the
    model in the legend.
  }
  \item{...}{arguments to be passed to methods; mainly graphical
    arguments, or further models in case of class \code{'RMmodel'},
    see Details.
  }
  \item{dim}{must equal 1 or 2; only for \code{class(x)=="RMmodel"}; the
    covariance function and the variogram are plotted as a function of
    \eqn{\R^\code{dim}}{R^\code{dim}}.}
  \item{n.points}{integer; only for \code{class(x)=="RMmodel"}; the
    number of points at which the model 
    evaluated (in each dimension); defaults to 200} 
  \item{fct.type}{character; only for \code{class(x)=="RMmodel"}; must
    equal \code{NULL}, \code{"Cov"} or \code{"Variogram"}; controls
    whether the covariance (\code{fct.type="Cov"}) or the
    variogram (\code{fct.type="Variogram"}) is plotted; \code{NULL}
    implies automatic choice, where \code{"Cov"} is chosen whenever the model is
    stationary} 
}


\section{Creating Objects}{
 Objects are created by calling a function of class
 \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
}

\section{Slots}{
 \describe{
 \item{\code{call}:}{language object; the function call by which the
 object was generated} 
 \item{\code{name}:}{character string; nickname of the model, name of
 the function by which the object was generated} 
 \item{\code{submodels}:}{list; contains submodels (if existent)} 
 \item{\code{par.model}:}{list; conatins model specific arguments} 
 \item{\code{par.general}:}{list of 4; contains the four standard
 arguments \code{var}, \code{scale}, \code{Aniso} and \code{proj} that
 can be given for any model; if not specified by the user, the string
 \code{"RFdefault"} is inserted}
 }
}
 
%\section{Extends}{
%}

\section{Methods}{
 \describe{
 \item{+}{\code{signature(x = "RMmodel")}: allows to sum up
 covariance models; internally calls \command{\link{RMplus}}.}
\item{-}{\code{signature(x = "RMmodel")}: allows to substract
 covariance models; internally calls \command{\link{R.minus}}.}
 \item{*}{\code{signature(x = "RMmodel")}: allows to multiply
 covariance models; internally calls \command{\link{R.minus}}.}
  \item{/}{\code{signature(x = "RMmodel")}: allows to devide
 covariance models; internally calls \command{\link{R.div}}.}
 \item{c}{\code{signature(x = "RMmodel")}: concenates
 covariance functions or variogram models}
 \item{plot}{\code{signature(x = "RMmodel")}: gives a plot of the
 covariance function or of the variogram model, for more details see
 \command{\link{plot-method}}.}
 \item{points}{\code{signature(x = "RMmodel")}: adds a covariance plot to
 an existing plot, for more details see
 \command{\link{plot-method}}.}
 \item{lines}{\code{signature(x = "RMmodel")}: adds a covariance plot to
 an existing plot, for more details see
 \command{\link{plot-method}}.}
 \item{str}{\code{signature(x = "RMmodel")}: as the usual
 \code{\link[utils]{str}}-method for S4 objects but where only those
 entries of the 
 'par.general'-slot are shown that contain values different from
 \code{'RFdefault'}}
 \item{show}{\code{signature(x = "RMmodel")}: returns the structure
 of \code{x}}
 \item{print}{\code{signature(x = "RMmodel")}: identical with
 \command{show}-method, additional argument is \code{max.level}}
 \item{[}{\code{signature(x = "RMmodel")}: enables accessing
 the slots via the \code{"["}-operator, e.g. \code{x["par.general"]}}
 \item{[<-}{\code{signature(x = "RMmodel")}: enables replacing
 the slots via the \code{"["}-operator}
 
\item{\code{signature(x = "RMmodel", y = "missing")}}{Generates
 covariance function or variogram function plots in one or two
 dimensions.}

 }
}

%\section{Details}{
%}

\author{Alexander Malinowski \email{Alexander.Malinowski@web.de}}

\seealso{
  \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
  \link{RMmodel}
}

\keyword{classes}
\keyword{print}
\keyword{hplot}


\examples{
# see RMmodel for introductory examples


# Compare:
model <- RMexp(scale=2) + RMnugget(var=3)
str(model)  ## S4 object as default in version 3 of RandomFields

model <- summary(model)
str(model)  ## list style as in version 2 of RandomFields
            ## see also 'spConform' in 'RFoptions' to make this style
            ## the default


\dontshow{FinalizeExample()}
}

