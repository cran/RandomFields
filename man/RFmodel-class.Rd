\name{RMmodel-class}
\docType{class}
\alias{RMmodel-class}
\alias{show,RMmodel-method}
\alias{print.RMmodel}
\alias{str.RMmodel}
\alias{[,RMmodel-method}
\alias{[,RMmodel,ANY,ANY,ANY-method}
\alias{[<-,RMmodel-method}
\alias{[<-,RMmodel,ANY,ANY,ANY-method}
\alias{+,RMmodel,RMmodel-method}
%\alias{plot,RMmodel-method}
%\alias{points.RMmodel}
%\alias{lines.RMmodel}

\title{Class \code{RMmodel}}
\description{ Class for \pkg{RandomField}'s representation of explicit
 covariance models 
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
 \item{\code{par.model}:}{list; conatins model specific parameters} 
 \item{\code{par.general}:}{list of 4; contains the four standard
 parameters \code{var}, \code{scale}, \code{Aniso} and \code{proj} that
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
 
 }
}

%\section{Details}{
%}

\author{Alexander Malinowski \email{malinows@math.uni-goettingen.de}}

\seealso{
  \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
  \link{RMmodel}
}

\examples{
# see RMmodel for examples
\dontshow{RFoptions(seed=NA)}
}


\keyword{classes}
