\name{RMmodelsAdvanced}
\alias{RMmodelsAdvanced}
\alias{Advanced RMmodels}
\title{Advanced features of the mdoels}
\description{
 Here, further models and advanced comments for \command{\link{RMmodel}}
 are given. See also \command{\link{RFgetModelNames}}.
}

\details{

\bold{Further stationary and isotropic models}

\tabular{ll}{
\command{\link{RMaskey}} \tab Askey model (generalized test or triangle model) \cr
\command{\link{RMbessel}} \tab Bessel family \cr
\command{\link{RMcircular}} \tab circular model \cr
\command{\link{RMcauchy}} \tab modified Cauchy family \cr
\command{\link{RMcubic}} \tab cubic model (see Chiles \& Delfiner) \cr
\command{\link{RMdagum}} \tab Dagum model \cr
\command{\link{RMdampedcos}} \tab exponentially damped cosine \cr
\command{\link{RMqexp}} \tab Variant of the exponential model \cr
\command{\link{RMid}} \tab identity \cr
\command{\link{RMfractdiff}} \tab fractionally differenced process \cr
\command{\link{RMfractgauss}} \tab fractional Gaussian noise \cr
\command{\link{RMgengneiting}} \tab generalized Gneiting model \cr
\command{\link{RMgneitingdiff}} \tab Gneiting model for tapering \cr
\command{\link{RMhyperbolic}} \tab generalised hyperbolic model \cr
\command{\link{RMlgd}} \tab Gneiting's local-global distinguisher\cr
\command{\link{RMma}} \tab one of Ma's model \cr
\command{\link{RMpenta}} \tab penta model (see Chiles \& Delfiner) \cr
\command{\link{RMpower}} \tab Golubov's model \cr
\command{\link{RMwave}} \tab cardinal sine \cr
}

\bold{Variogram models (stationary increments/intrinsically stationary)}

\tabular{ll}{
\command{\link{RMdewijsian}} \tab generalised version of the DeWijsian model \cr
\command{\link{RMgenfbm}} \tab generalized fractal Brownian motion \cr
}

\bold{General composed models (operators)}

Here, composed models are given that can be of any kind, depending on the submodel.

\tabular{ll}{
% \command{\link{RMCauchy}} \tab Cauchy like transform -- TO BE PROGRAMMED (includes \code{ma1})) \cr
\command{\link{RMbernoulli}} \tab Correlation function of a binary field
based on a Gaussian field \cr
 \command{\link{RMexponential}} \tab exponential of a covariance model \cr
 \command{\link{RMintexp}} \tab integrated exponential of a covariance model (INCLUDES \code{ma2})\cr
 \command{\link{RMpower}} \tab powered variograms\cr
 \command{\link{RMqam}} \tab Porcu's quasi-arithmetric-mean model\cr
 \command{\link{RMS}} \tab for details on the optional transformation
 parameters (\code{var}, \code{scale}, \code{Aniso}, \code{proj}).
}

\bold{Stationary and isotropic composed models (operators)}

\tabular{ll}{
 \command{\link{RMcutoff}} \tab Gneiting's modification towards finite range\cr
 \command{\link{RMintrinsic}} \tab Stein's modification towards finite range\cr
 \command{\link{RMnatsc}} \tab natural pratical range operator\cr
 \command{\link{RMstein}} \tab Stein's modification towards finite range\cr 
% \command{\link{RMtbm2}} \tab Turning bands operator in two (spatial)
% dimensions\cr % nicht an user exportiert
 \command{\link{RMtbm}}\tab Turning bands operator in three (spatial) dimensions
}

\bold{Stationary space-time models}

Here, most of the models are composed models (operators).
\tabular{ll}{
 \command{\link{RMave}} \tab space-time moving average model \cr
 \command{\link{RMcoxisham}} \tab Cox-Isham model \cr
 \command{\link{RMcurlfree}} \tab curlfree (spatial) field (stationary and anisotropic)\cr
 \command{\link{RMdivfree}} \tab divergence free (spatial) vector valued field, (stationary and anisotropic)\cr
 \command{\link{RMiaco}} \tab non-separabel space-time model\cr
 % obsolete -- included by Cauchy --- should be given by an example
 \command{\link{RMmastein}} \tab Ma-Stein model\cr
 \command{\link{RMnsst}} \tab Gneiting's non-separable space-time model \cr 
 \command{\link{RMstein}} \tab Stein's non-separabel space-time model\cr
 \command{\link{RMstp}} \tab Single temporal process\cr
 \command{\link{RMtbm}} \tab Turning bands operator}

\bold{Multivariate and vector valued models}
\tabular{ll}{
 \command{\link{RMbiwm}} \tab full bivariate Whittle-Matern model (stationary and isotropic)\cr
 \command{\link{RMbigneiting}} \tab bivariate Gneiting model (stationary and isotropic)\cr
\command{\link{RMcurlfree}} \tab curlfree (spatial) vector-valued field (stationary and anisotropic)\cr
\command{\link{RMdelay}} \tab bivariate delay effect model (stationary)\cr
\command{\link{RMdivfree}} \tab divergence free (spatial) vector valued field, (stationary and anisotropic)\cr
 \command{\link{RMkolmogorov}} \tab Kolmogorov's model of turbulence\cr
\command{\link{RMmatrix}} \tab trivial multivariate model\cr
%\command{\link{RMmqam}} \tab multivariate quasi-arithmetic mean (stationary)\cr
 \command{\link{RMparswm}} \tab multivariate Whittle-Matern model (stationary and isotropic)\cr
\command{\link{RMschur}} \tab element-wise product with a positive definite
matrix\cr 
\command{\link{RMvector}} \tab vector-valued field (combining \command{\link{RMcurlfree}} and \command{\link{RMdivfree}})
}

\bold{Non-stationary models}
\tabular{ll}{
 \command{\link{RMnonstwm}} \tab one of Stein's non-stationary Wittle-Matern model \cr
}



\bold{Models related to max-stable random fields (tail correlation functions)}
\tabular{ll}{
\command{\link{RMaskey}} \tab Askey model (generalized test or triangle
model) with \eqn{\alpha \ge [dim / 2] +1}\cr
\command{\link{RMbessel}} \tab Bessel family \cr
 \command{\link{RMbernoulli}} \tab Correlation function of a binary field
 based on a Gaussian field \cr
 \command{\link{RMbr2bg}} \tab Operator relating a Brown-Resnick process
 to a Bernoulli process\cr
 \command{\link{RMbr2eg}} \tab Operator relating a Brown-Resnick process
 to an extremal Gaussian process\cr
 \command{\link{RMbrownresnick}} \tab tail correlation function
 of Brown-Resnick process\cr
\command{\link{RMgencauchy}} \tab generalized Cauchy family with  \eqn{\alpha\le 1/2}\cr
 \command{\link{RMmatern}} \tab Whittle-Matern model with \eqn{\nu\le 1}\cr
\command{\link{RMschlather}} \tab tail correlation function of the
 extremal Gaussian field \cr
 \command{\link{RMstable}} \tab symmetric stable family or powered
 exponential model with \eqn{\alpha\le 1}\cr
 \command{\link{RMstrokorb}} \tab shapes functions related max-stable
 processes\cr 
 \command{\link{RMwhittle}} \tab Whittle-Matern model, alternative
 parametrization with \eqn{\nu\le 1/2}\cr

 }


\bold{Other covariance models}
\tabular{ll}{
 \command{\link{RMuser}} \tab User defined model \cr
}

\bold{Auxiliary models}
There are models or better function that are not covariance functions,
but can be part of a model definition. See \bold{\link{Auxiliary RMmodels}.}
}

\note{
  \itemize{
    \item
    Note that, instead of the named parameters, a single parameter \code{k}
    can be passed. This is possible if all the parameters
    are scalar. Then \code{k} must have length equal to the number of
    parameters.

    \item
    Note also that a completely different possibility exists to define a
    model, namely by a list. This format allows for easy flexible models
    and modifications. Here, the parameter \code{var}, \code{scale},
    \code{Aniso} and \code{proj} must be passed by the model
    \command{\link{RMS}}. 
    For instance,
    \code{model <- RMexp(scale=2, var=5)} is equivalent to
    \code{model <- list("RMS", scale=2, var=5, list("RMexp"))}.
    and
    \code{model <- RMnsst(phi=RMgauss(var=7), psi=RMfbm(alpha=1.5),
      scale=2, var=5)}
    is equivalent to
    \code{model <- list("RMS", scale=2, var=5,
      list("RMnsst",
      phi=list("RMS", var=7, list("RMgauss")),
      psi=list("RMfbm", alpha=1.5))
    )}.

  \item
  Instead of a fixed parameter, a random parameter might be given.
  These are models starting with \code{RR} or is distribution family,
  e.g. \code{norm}, \code{exp}, or \code{unif}. Note that the behaviour
  of the random parameters varies:

  \itemize{
    \item in Max-stable fields and
    \command{\link{RPpoisson}}, a new realisation of the
    random parameter is drawn for each shape function
    \item in all the other cases: a realisation is only drawn once.
    This effects, in particular, Gaussian fields with parameter
    \code{n>1}, where all the realisations are based on the same
    realisation of the random parameters.
  }
  
  MLE ist not programmed yet.

  Very advanced:
  In case of a distribution family, its parameters might be again given
  by a \link{RMmodel}. Note that the checking of the validity of the
  parameters is not possible anymore, in general.
  If the size of the parameter can be passed by given in the model first
  a deterministic parameter of the given size (with dummy values) and
  then the genuine (random) definition.

  See also
  \link{RMmodelsAuxiliary} and \link{Baysian}.
  }
  
  All models have secondary names that stem from 
  \pkg{RandomFields} versions 2 and earlier and
  that can also be used as strings in the list notation.
  See \code{\link{RFgetModelNames}(internal=FALSE)} for
  the full list. 
}


%\section{Methods}{
% \describe{
% \item{[}{\code{signature(x = "RFgridDataFrame")}: selects
% slot by name}
% \item{[<-}{\code{signature(x = "RFgridDataFrame")}: replaces
% slot by name}
% \item{as}{\code{signature(x = "RFgridDataFrame")}:
% converts into other formats, only implemented for target class
% \command{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}} } 
% \item{cbind}{\code{signature(...)}: if arguments have identical
% topology, combine their attribute values}
% }
%}


\references{
 \itemize{
 \item Chiles, J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 % \item Gneiting, T. and Schlather, M. (2004)
 % Statistical modeling with covariance functions.
 % \emph{In preparation.}
 \item Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University.
 \item Schlather, M. (2011) Construction of covariance functions and
 unconditional simulation of random fields. In Porcu, E., Montero, J.M.
 and Schlather, M., \emph{Space-Time Processes and Challenges Related
 to Environmental Problems.} New York: Springer.
 % \item Schlather, M. (2002) Models for stationary max-stable
 % random fields. \emph{Extremes} \bold{5}, 33-44.
 \item Yaglom, A.M. (1987) \emph{Correlation Theory of Stationary and
 Related Random Functions I, Basic Results.}
 New York: Springer.
 \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 3nd edition.
 }
}

\seealso{\command{\link{RFformula}}, \command{\link{RMmodels}},
 \command{\link{RMmodelsAuxiliary}}
}

\examples{
set.seed(0)
RFgetModelNames(type="positive", group.by=c("domain", "isotropy"))
}

\author{
 Alexander Malinowski, \email{malinowski@math.uni-mannheim.de}
 
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}
