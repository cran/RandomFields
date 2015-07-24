\name{RMintern}
\alias{RMintern}

% to do: check the list with the make_internal and internal=TRUE

\alias{RMmixed}
\alias{ROmissing}
\alias{RMnull}
\alias{RMnatsc_intern}

\alias{RPplusproc}
\alias{RP$pow}
\alias{RPaverageIntern}
\alias{RPcutoffIntern}
\alias{RPintrinsIntern}
\alias{RPhyperIntern}
\alias{RPnuggetIntern}
\alias{RPspectralIntern}
\alias{RPtbmIntern}
\alias{RPbrorigIntern}
\alias{RPbrmixedIntern}
\alias{RPbrshiftIntern}

\alias{RO#}
\alias{RO>}

\alias{RMselect}
\alias{RMSpower}
\alias{RMmult_inverse} 
\alias{RMsetparam}
\alias{RMshape.ave}
\alias{RMshape.stp}
\alias{RMr3binner}

\alias{RRarcsqrt}
\alias{RRsetDistr}
\alias{RMptsGivenShape}
\alias{RMstandardShape}
\alias{RMstatShape}

\alias{RPS}
\alias{RPplus}
\alias{RPmult}
\alias{RPmppplus}
\alias{RPtrend}

\title{Internal models}
\description{
  Internal models or model names that may appear in feedbacks from
  'RandomFields'. Those ending by \sQuote{Intern} should appear only
  in very rare cases.
}
\details{
The following and many more internal models exist
\itemize{
  \item RF__Name__ : internal representation of certain functions
  \link[=RF]{RF__name__}

  \item \code{RO#} : model for transforming coordinates within the cartesian
   system

   \item \code{RO>} : model for transforming earth coordinates to cartesian
   coordinates  

   \item \code{ROmissing} : for error messages only

   \item \code{RMmixed} : internal representation of a \link[=formula]{mixed model}

   \item \code{RMselect} : will be obsolete in future

   \item \code{RMsetparam}, \code{RMptsGivenShape}, \code{RMstandardShape}, \code{RMstatShape} :
   for max-stable processes and Poisson processes:
   models that combine shape functions with corresponding point processes

   \item \code{RP__name__Intern} : internal representations of some \link[=RP]{processes}

   \item \code{RPS}, \code{RPplusp}, etc. : specific processes for
   \link{RMS} and  \link{RMplus} etc. (For those covariance models that
   have specific simulation processes programmed.)

   \item \command{\link{RMS}} : internal representation of
   the modifying arguments \code{var}, \code{scale}, \code{Aniso},
   \code{proj}
 }
 }

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\keyword{spatial}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## in the following 'try' the model 'RMtbmIntern' appears
model <- RPtbm(RMexp())
x <- seq(0, 10, 1) 
try(RFsimulate(model, x=x)) ## fails
\dontshow{FinalizeExample()}

}
