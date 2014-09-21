\name{RFgetModelInfo}
\alias{RFgetModel}
\alias{RFgetModelInfo}
%\alias{GetModelInfo}
%\alias{GetModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal information}
\description{
 The function returns internal information about the simulation
 of a random field and the calculation of covariance functions
}
\usage{
RFgetModel(register, explicite.natscale, show.call=FALSE)

RFgetModelInfo(register, level = 3,
               spConform = RFoptions()$general$spConform,
               which.submodels = c("user", "internal", "both"),
               modelname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{register}{\eqn{0,...,21} or an evaluating function,
   e.g. \command{\link{RFsimulate}}.
   Place where intermediate calculations are stored.
   See also section \code{Registers} in \command{\link{RFoptions}}.
 }
 \item{explicite.natscale}{logical. Advanced option.
   If missing, then the model is returned as stored.
   If \code{FALSE} then any \command{\link{RMnatsc}} is ignored.
   If \code{TRUE} then any  \command{\link{RMnatsc}} is tried to be
   combined with leading \command{\link{RMS}}, or returned as such.
 }
 \item{show.call}{logical. If \code{FALSE} then the model
   is shown as interpreted. If \code{TRUE} then the user's input
   including the calling function is returned. See example below.
 }
 \item{level}{integer [0..5]; level of details, i.e. the higher the
   number the more details are given. 
   If \code{level >= 10} then the leading internal model is also given
   (which is, in general, not of interest by the user).
 }

% \item{max.elements}{integer; since \command{\link{RFgetModelInfo}}
% might be a dump of simulation that needs a huge amount of
% memory and since all entries are copied, the maximal amount
% of available memory might be easily exceeded. Therefore,
% only the size of the critical parts are returned and not the
% vector or matrix itself, if the number of elements exceeds
% \code{max.elements}. 
% }
 \item{spConform}{see \command{\link{RFoptions}}
 }
 \item{which.submodels}{
   Internally, the sub-models are represented in two different ways
   \sQuote{internal} and \sQuote{user}. The latter is very close to
   the model defined by the user.
 }
 \item{modelname}{ string. If \code{modelname} is given then
 it returns the first appearance of the
 covariance model with name \code{modelname}.
 If \code{meth} is given then the model within the method is returned.
 }
}
\details{
 \code{\link{RFgetModelInfo}(register, ignore.active=TRUE)} is
 useful for debugging and specialists' need to control
 the algorithm, see the examples in \command{\link{RFoptions}}
 and \command{\link{RFsimulate}}.
 
 If \command{\link{RFoptions}}\code{()$Storage=FALSE} then values of
 the internal registers are not kept if \command{\link{RFsimulate}}
 has been called.
 Hence \command{\link{RFgetModelInfo}} cannot
 provide any information.

% \code{GetModelInfo} returns the complete information on the internal
% model structure. It allows for \code{register=-1} returning the
% model structure for the last use of \command{\link{RMmodel}},
% \link{sophisticated} models
% or \command{\link{Variogram}} or similar commands.
% \code{register=-2} is for internal use only.

% \code{GetModel} returns only the parts of the internal model structure
% that have been defined by the user.
% The \code{modus} is ignored if \code{PracticalRange=FALSE}
% in \command{\link{RFoptions}}. In case \code{PracticalRange=TRUE}
% the \code{modus} has the following effects
% \enumerate{
% \item{0}{the model is returned without the explicite scale transformation}
% \item{1}{the model is returned the way it is stored, including the
% scale transformation}
% \item{2}{the scale
% transformation of the the practical range is included,
% but a simplified model will be returned
% }
% }
}
\note{
 Put \code{Storing=TRUE}, see \command{\link{RFoptions}}
 if you like to have more
 internal information in case of an expected failure of an initialisation of
 a random field simulation.
}
\value{
 List of internal information is returned.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\seealso{\command{\link{RFsimulate}}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMexp(scale=4, var=2) + RMnugget(var=3) + RMtrend(mean=1)
z <- RFsimulate(model, 1:4, storing=TRUE)
RFgetModelInfo(RFsimulate, level=1)


RFgetModel(RFsimulate, show.call=FALSE)
RFgetModel(RFsimulate, show.call=TRUE)
\dontshow{FinalizeExample()}

}
\keyword{spatial}
