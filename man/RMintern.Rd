\name{RMintern}
\alias{RMintern}
\alias{RMmixed}
\alias{RMmissing}
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
\alias{RM#}
\alias{RMselect}
\alias{RMset}
\alias{RMptsGivenShape}
\alias{RMstandardShape}
\alias{RMpoissonShape}
\alias{RFSimulate}
\alias{RFCov}
\alias{RFCovMatrix}
\alias{RFDummy}
\alias{RFget}
\alias{RFFctn}
\alias{RFDistr}
\alias{RFPseudovariogram}
\alias{RFVariogram}
\alias{RP$proc}
\alias{RPplusproc}
\title{Internal models}
\description{
  Internal models that may appear in feedbacks from
  'RandomFields'
}
\details{
The following internal models exist
\itemize{
  \item RF__Name__ : internal representation of certain functions
  \link[=RF]{RF__name__}
   \item RM# : model only for transforming coordinates
  \item RMmissing : for error messages only
  \item RMmixed : internal representation of a \link[=formula]{mixed model}
   \item RMselect : will be obsolete in future
   \item RMset, RMptsGivenShape, RMstandardShape, RMpoissonShape :
   for max-stable processes and Poisson processes:
   models that combine shape functions with corresponding point processes
   \item RP__name__Intern : internal representations of \link[=RP]{processes}
   \item RP$proc, RPplusproc : internal models for \link{RMS} and
   \link{RMplus}
 }
 }

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}
\examples{
set.seed(0)

## in the following 'try' the model 'RMtbmIntern' appears
model <- RPtbm(RMexp())
x <- seq(0, 10, 1) 
try(RFsimulate(model, x=x)) ## fails

}