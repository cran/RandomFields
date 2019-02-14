\name{RFgetModel}
\alias{RFgetModel}

%\alias{GetModelInfo}
%\alias{GetModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internally stored model}
\description{
 The function returns the stored model.
}
\usage{
RFgetModel(register, explicite.natscale, show.call=FALSE,
           origin="original")
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
 \item{show.call}{logical or character. If \code{FALSE} then the model
   is shown as interpreted. If \code{TRUE} then the user's input
   including the calling function is returned. See example below.

   If \code{show.call} is a character it behaves as
   \code{\link[=RFgetModelInfo]{which.submodels}}.
 }
 \item{origin}{\argOrigin}

} 
\details{
  Whereas
  \command{\link{RFgetModel}} returns a model that can be re-used
  by the user,
  \command{RFgetModelInfo} can return detailed information.
}
\note{
 Put \code{Storing=TRUE}, see \command{\link{RFoptions}},
 if you like to have (more)
 internal information in case of failure of an initialization of
 a random field simulation.
}
\value{
  The stored model is returned in list format.
}
\me
\seealso{\command{\link{RFgetModelInfo}}, \command{\link{RFsimulate}}.
}
\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMexp(scale=4, var=2) + RMnugget(var=3) + RMtrend(mean=1)
z <- RFsimulate(model, 1:4)
RFgetModel(show.call=FALSE)
RFgetModel(show.call=TRUE)
\dontshow{FinalizeExample()}}
\keyword{spatial}
