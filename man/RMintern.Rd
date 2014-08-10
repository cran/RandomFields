\name{RMintern}
\alias{RMintern}

\alias{RMmixed}
\alias{RMmissing}
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
\alias{RM$pow}
\alias{RMU}
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

\alias{RFCov}
\alias{RFCovMatrix}
\alias{RFDummy}
%\alias{RFget}
\alias{RFFctn}
\alias{RFDistr}
\alias{RFPseudovariogram}
\alias{RFSimulate}
\alias{RFVariogram}
\alias{RFloglikelyhood}

\alias{RP$proc}
\alias{RPplusproc}
\alias{RPmultproc}
\alias{RPmppplusproc}
\alias{RPtrendproc}

\title{Internal models}
\description{
  Internal models that may appear in feedbacks from
  'RandomFields'
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

   \item \code{RMmissing} : for error messages only

   \item \code{RMmixed} : internal representation of a \link[=formula]{mixed model}

   \item \code{RMselect} : will be obsolete in future

   \item \code{RMsetparam}, \code{RMptsGivenShape}, \code{RMstandardShape}, \code{RMpoissonShape} :
   for max-stable processes and Poisson processes:
   models that combine shape functions with corresponding point processes

   \item \code{RP__name__Intern} : internal representations of \link[=RP]{processes}

   \item \code{RP$proc}, \code{RPplusproc}, etc. : specific processes for \link{RMS} and
   \link{RMplus} etc

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