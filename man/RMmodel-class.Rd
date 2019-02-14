\name{RMmodel-class}
\docType{class}
\alias{RMmodel-class}
\alias{show,RMmodel-method}
\alias{print.RMmodel}
\alias{plot,RMmodel,missing-method}
\alias{persp,RMmodel-method}
\alias{image,RMmodel-method}
\alias{points,RMmodel-method}
\alias{lines,RMmodel-method}
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
\alias{RFplotModel}


\title{Class \code{RMmodel}}
\description{ Class for \pkg{RandomFields'} representation of explicit
 covariance models 
}


\usage{
RFplotModel(x, y, dim=1,
           n.points=if (dim==1 || is.contour) 200 else 100,
           fct.type=NULL, MARGIN, fixed.MARGIN, maxchar=15, ...,
           plotmethod=if (dim==1) "matplot" else "contour")

\S4method{plot}{RMmodel,missing}(x, y, ...)
\S4method{points}{RMmodel}(x, ..., type="p")
\S4method{lines}{RMmodel}(x, ..., type="l")
\S4method{image}{RMmodel}(x, ..., dim=2)
\S4method{persp}{RMmodel}(x, ..., dim=2, zlab="")
}



\arguments{
  \item{x}{object of class \code{\link[=RFsp-class]{RFsp}} or
    \command{\link[=RFempVariog-class]{RFempVario}} or
    \command{\link[=RFfit-class]{RFfit}} or
    \command{\link[=RMmodel-class]{RMmodel}}; in the latter case, \code{x} can
    be any sophisticated model but it must be either stationary or a
    variogram model.} 
  \item{y}{ignored in most methods}
  \item{MARGIN}{vector of two; two integer values giving the coordinate
    dimensions 
    w.r.t. whether the field or the covariance model
    is to be plotted; in all other directions,
    the first index is taken.}
%  \item{MARGIN.slices}{integer value; if \eqn{[space-time-dimension>2]},
%    \code{MARGIN.slices} can specify a third dimension w.r.t. which a
%    sequence of slices is plotted. Currently only works for grids.} 
  \item{fixed.MARGIN}{only for \code{class(x)==CLASS_CLIST} and if
    \code{dim > 2}; a vector of length \code{dim}-2 with distance values
    for the coordinates that are not displayed.}
  \item{maxchar}{integer. Maximum number of characters to print the
    model in the legend.
  }
  \item{...}{arguments to be passed to methods; mainly graphical
    arguments, or further models in case of class \code{CLASS_CLIST},
    see Details.
  }
  \item{dim}{must equal 1 or 2; only for \code{class(x)==CLASS_CLIST}; the
    covariance function and the variogram are plotted as a function of
    \eqn{\R^\code{dim}}{R^\code{dim}}.}
  \item{n.points}{integer; only for \code{class(x)==CLASS_CLIST}; the
    number of points at which the model is 
    evaluated (in each dimension); defaults to 200.} 
  \item{fct.type}{character; only for \code{class(x)==CLASS_CLIST}; must
    equal \code{NULL}, \code{"Cov"} or \code{"Variogram"}; controls
    whether the covariance (\code{fct.type="Cov"}) or the
    variogram (\code{fct.type="Variogram"}) is plotted; \code{NULL}
    implies automatic choice, where \code{"Cov"} is chosen whenever the model is
    stationary.}
  \item{plotmethod}{string or function. Internal.
  }
  \item{type}{character. See \command{\link[graphics]{points}}.}
  \item{zlab}{character. See \command{\link[graphics]{persp}}.}
}

\value{
  If \code{RFoptions()$split_screen=TRUE} and
  \code{RFoptions()$close_screen=TRUE} then
  the \code{plot} functions return
  the \code{screen} numbers. Else \code{NULL}.
}

\section{Creating Objects}{
 Objects are created by calling a function of class
 \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}.
}

\section{Slots}{
 \describe{
 \item{\code{call}:}{language object; the function call by which the
 object was generated} 
 \item{\code{name}:}{character string; nickname of the model, name of
 the function by which the object was generated} 
 \item{\code{submodels}:}{list; contains submodels (if existent)} 
 \item{\code{par.model}:}{list; contains model specific arguments} 
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
 \item{+}{\code{signature(x = CLASS_CLIST)}: allows to sum up
 covariance models; internally calls \command{\link{RMplus}}.}
\item{-}{\code{signature(x = CLASS_CLIST)}: allows to substract
 covariance models; internally calls \command{\link{R.minus}}.}
 \item{*}{\code{signature(x = CLASS_CLIST)}: allows to multiply
 covariance models; internally calls \command{\link{R.mult}}.}
  \item{/}{\code{signature(x = CLASS_CLIST)}: allows to divide
 covariance models; internally calls \command{\link{R.div}}.}
 \item{c}{\code{signature(x = CLASS_CLIST)}: concatenates
 covariance functions or variogram models.}
 \item{plot}{\code{signature(x = CLASS_CLIST)}: gives a plot of the
 covariance function or of the variogram model; for more details see
 \command{\link{plot-method}}.}
 \item{points}{\code{signature(x = CLASS_CLIST)}: adds a covariance plot to
 an existing plot; for more details see
 \command{\link{plot-method}}.}
 \item{lines}{\code{signature(x = CLASS_CLIST)}: adds a covariance plot to
 an existing plot; for more details see
 \command{\link{plot-method}}.}
 \item{str}{\code{signature(x = CLASS_CLIST)}: as the usual
 \code{\link[utils]{str}}-method for S4 objects but only those
 entries of the 
 'par.general'-slot are shown that contain values different from
 \code{'RFdefault'}.}
 \item{show}{\code{signature(x = CLASS_CLIST)}: returns the structure
 of \code{x}.}
 \item{print}{\code{signature(x = CLASS_CLIST)}: identical with
 \command{show}-method, additional argument is \code{max.level}.}
 \item{[}{\code{signature(x = CLASS_CLIST)}: enables accessing
 the slots via the \code{"["}-operator, e.g. \code{x["par.general"]}.}
 \item{[<-}{\code{signature(x = CLASS_CLIST)}: enables replacing
 the slots via the \code{"["}-operator.}
 
\item{\code{signature(x = CLASS_CLIST, y = "missing")}}{Generates
 covariance function or variogram function plots in one or two
 dimensions.}

 }
}

\section{Details}{
  All the above arguments apply for all the  S3 and S4
  functions given here as they call \command{RFplotModel} immediately.
}

\author{Alexander Malinowski, \martin}

\seealso{
  \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
  \link{RMmodel}
}

\keyword{classes}
\keyword{print}
\keyword{hplot}


\examples{\dontshow{StartExample()}
# see RMmodel for introductory examples


# Compare:
model <- RMexp(scale=2) + RMnugget(var=3)
str(model)  ## S4 object as default in version 3 of RandomFields

model <- summary(model)
str(model)  ## list style as in version 2 of RandomFields
            ## see also 'spConform' in 'RFoptions' to make this style
            ## the default

\dontshow{FinalizeExample()}}

