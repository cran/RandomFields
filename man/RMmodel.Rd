\name{RMmodel}
\alias{RM}
\alias{RMmodel}
\alias{RMmodels}
\alias{[,RMmodel,ANY,ANY-method} %]
\alias{[<-,RMmodel,ANY,ANY-method} %]
\title{Covariance and Variogram Models in \pkg{RandomFields} (RM commands)}
\description{
Summary of implemented covariance and variogram models% in \link{RFformula}
}

\details{
 To generate a covariance or variogram model for use within
 \pkg{RandomFields}, calls of the form
 \deqn{RM_name_(..., var, scale, Aniso, proj)}
 can be used,
 where _name_ has to be replaced by a valid model name,
 \itemize{
 \item
 \code{...} can take model specific arguments. %Argument
 %corresponding to specific covariance model
 \item
 \code{var} is the optional variance argument \eqn{v},
 \item
 \code{scale} the optional scale argument \eqn{s},
 \item
 \code{Aniso} an optional anisotropy matrix \eqn{A} or given by \command{\link{RMangle}}, and
 \item
 \code{proj} is the optional projection vector which defines a
 diagonal matrix of zeros and ones and \code{proj} gives the
 positions of the ones (integer values).
 }
 With \eqn{\phi} denoting the original model, the transformed model is
 \eqn{C(h) = v * \phi(A*h/s)}.

 \command{RM_name_} must be a function of class
 \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}.
 The return value of all functions \command{RM_name_} is of class
 \command{\link[=RMmodel-class]{RMmodel}}.\cr
 
 The following models are available
 (cf. \command{\link{RFgetModelNames}}).
 %Choose from the following covariance models

 \bold{Basic stationary and isotropic models}
 \tabular{ll}{
 \command{\link{RMcauchy}} \tab Cauchy family \cr
 \command{\link{RMexp}} \tab exponential model \cr
 \command{\link{RMgencauchy}} \tab generalized Cauchy family \cr
 \command{\link{RMgauss}} \tab Gaussian model \cr
 \command{\link{RMgneiting}} \tab differentiable model with compact support \cr
 \command{\link{RMmatern}} \tab Whittle-Matern model \cr
 \command{\link{RMnugget}} \tab nugget effect model \cr
 \command{\link{RMspheric}} \tab spherical model \cr
 \command{\link{RMstable}} \tab symmetric stable family or powered exponential model \cr
 \command{\link{RMwhittle}} \tab Whittle-Matern model, alternative
 parametrization\cr 
 }

 \bold{Variogram models (stationary increments/intrinsically stationary)}

 \tabular{ll}{
 \command{\link{RMfbm}} \tab fractal Brownian motion\cr
 }


 \bold{Basic Operations}

 \tabular{ll}{
 \command{\link{RMmult}}, \code{*} \tab product of covariance models \cr
 \command{\link{RMplus}}, \code{+} \tab sum of covariance models or variograms\cr
 }



 \bold{Basic models for mixed effect modelling}

 \tabular{ll}{
 \command{\link{RMfixcov}} \tab constant pre-defined covariance \cr
 \command{\link{RMfixed}} \tab fixed or trend effects;
 Caution: \link{RMfixed} is not
 a function and can be used only in \link[=RFformula]{formula notation}\cr
 %\command{\link{RMmixed}} \tab Mixture of fixed, mixed, and random effect
 %model\cr % gibts nur noch intern
 }

 \bold{Others}
 \tabular{ll}{
   \command{\link{RMtrend}} \tab trend \cr
   \command{\link{RMangle}} \tab defines a 2x2 anisotropy matrix by
   rotation and stretch arguments.
 }

  
 \bold{See \link{RMmodelsAdvanced} for many more, advanced models.\cr
 \bold{See \link{spherical models} for models valid on spherical
   coordinate systems.}
 }
}

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
 \item Yaglom, A.M. (1987) \emph{Correlation Theory of Stationary and
 Related Random Functions I, Basic Results.}
 New York: Springer.
 \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 3nd edition.
 }
}

\author{
 Alexander Malinowski, \email{malinowski@math.uni-mannheim.de}
 
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}


\seealso{
  \link{RC}, \link{RF},  \link{RP}, \link{RR}, \link{R.},
  \command{\link{RFcov}},
  \command{\link{RFformula}},
  \command{\link{RMmodelsAdvanced}},
  \command{\link{RMmodelsAuxiliary}}
}

\keyword{spatial}
\keyword{models}


\examples{
% library(RandomFields, lib="~/TMP")
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
RFgetModelNames(type="positive definite", domain="single variable",
                isotropy="isotropic", operator=!FALSE) ## RMmodel.Rd

## an example of a simple model
model <- RMexp(var=1.6, scale=0.5) + RMnugget(var=0) #exponential + nugget
plot(model)

\dontshow{FinalizeExample()}
}