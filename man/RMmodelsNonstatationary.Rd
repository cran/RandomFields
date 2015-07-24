\name{RMmodelsNonstationary}
\alias{RMmodelsNonstationary}
\alias{Nonstationary RMmodels}
\alias{non-stationary RMmodels}
\alias{non-stationary}
\title{Nonstationary features of the mdoels}
\description{
 Here, nonstationary covariance models are presented.
}

\details{
\bold{Covariance models}
\tabular{ll}{
 \command{\link{RMnonstwm}} \tab one of Stein's non-stationary
 Wittle-Matern models
 \cr
 \command{\link{RMprod}} \tab scalar product
 \cr
 \code{\link[=RMS]{Aniso}} \tab for space transformation, see
 the example in \link{R.models}.
 \cr
 \code{scale}, cf. \command{\link{RMS}}, can be any non-negative
 function for any scale mixture model, such as
 the \link{whittle-matern}-classes,
 the \link{powered exponential} family, and the
 \command{\link{RMgencauchy}} model.
}

\bold{Trend models}
see \link{RMmodelsTrend}
 
}


\seealso{\command{\link{RFformula}}, \command{\link{RMmodels}},
   \command{\link{RM}},
 \command{\link{RMmodelsAdvanced}}

 \sQuote{\href{../doc/nonstationary_jss.pdf}{nonstationary}}, a vignette
 for nonstationary geostatistics
}

\author{
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}

\examples{
%RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
%##                   RFoptions(seed=NA) to make them all random again


# to do


\dontshow{FinalizeExample()}
}
