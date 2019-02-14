\name{RMmodelsAdvanced}
\alias{RMmodelsAdvanced}
\alias{Advanced RMmodels}
\title{Advanced features of the models}
\description{
 Here, further models and advanced comments for \command{\link{RMmodel}}
 are given. See also \command{\link{RFgetModelNames}}.
}

\details{

\bold{Further stationary and isotropic models}

\tabular{ll}{
\command{\link{RMaskey}} \tab Askey model (generalized test or triangle model) \cr
\command{\link{RMbcw}} \tab bridging model between
                       \command{\link{RMcauchy}} and \command{\link{RMgenfbm}} \cr
\command{\link{RMbessel}} \tab Bessel family \cr
\command{\link{RMcircular}} \tab circular model \cr
\command{\link{RMconstant}} \tab spatially constant model \cr
\command{\link{RMcubic}} \tab cubic model (see Chiles and Delfiner) \cr
\command{\link{RMdagum}} \tab Dagum model \cr
\command{\link{RMdampedcos}} \tab exponentially damped cosine \cr
\command{\link{RMqexp}} \tab variant of the exponential model \cr
\command{\link{RMfractdiff}} \tab fractionally differenced process \cr
\command{\link{RMfractgauss}} \tab fractional Gaussian noise \cr
\command{\link{RMgengneiting}} \tab generalized Gneiting model \cr
\command{\link{RMgneitingdiff}} \tab Gneiting model for tapering \cr
\command{\link{RMhyperbolic}} \tab generalized hyperbolic model \cr
\command{\link{RMlgd}} \tab Gneiting's local-global distinguisher\cr
\command{\link{RMlsfbm}} \tab locally stationary fractal Brownian motion \cr
\command{\link{RMpenta}} \tab penta model (see Chiles and Delfiner) \cr
\command{\link{RMpower}} \tab Golubov's model \cr
\command{\link{RMwave}} \tab cardinal sine \cr
}

\bold{Variogram models (stationary increments/intrinsically stationary)}

\tabular{ll}{
\command{\link{RMbcw}} \tab bridging model between
                       \command{\link{RMcauchy}} and \command{\link{RMgenfbm}} \cr
\command{\link{RMdewijsian}} \tab generalized version of the DeWijsian model \cr
\command{\link{RMgenfbm}} \tab generalized fractal Brownian motion \cr
\command{\link{RMflatpower}} \tab similar to fractal Brownian motion but
always smooth at the origin\cr
}

\bold{General composed models (operators)}

Here, composed models are given that can be of any kind (stationary/non-stationary), depending on the submodel.

\tabular{ll}{
% \command{\link{RMCauchy}} \tab Cauchy like transform -- TO BE PROGRAMMED (includes \code{ma1})) \cr
\command{\link{RMbernoulli}} \tab Correlation function of a binary field
based on a Gaussian field \cr
 \command{\link{RMexponential}} \tab exponential of a covariance model \cr
 \command{\link{RMintexp}} \tab integrated exponential of a covariance model (INCLUDES \code{ma2})\cr
 \command{\link{RMpower}} \tab powered variograms\cr
 \command{\link{RMqam}} \tab Porcu's quasi-arithmetic-mean model\cr
 \command{\link{RMS}} \tab details on the optional transformation
 arguments (\code{var}, \code{scale}, \code{Aniso}, \code{proj})
}

\bold{Stationary and isotropic composed models (operators)}

\tabular{ll}{
 \command{\link{RMcutoff}} \tab Gneiting's modification towards finite range\cr
 \command{\link{RMintrinsic}} \tab Stein's modification towards finite range\cr
 \command{\link{RMnatsc}} \tab practical range\cr
 \command{\link{RMstein}} \tab Stein's modification towards finite range\cr 
% \command{\link{RMtbm2}} \tab Turning bands operator in two (spatial)
% dimensions\cr % nicht an user exportiert
 \command{\link{RMtbm}}\tab Turning bands operator
}

\bold{Stationary space-time models}
\cr
See \link{RMmodelsSpaceTime}.

\bold{Non-stationary models}
\cr
See \link{RMmodelsNonstationary}.

\bold{Negative definite models that are not variograms}
\tabular{ll}{
\command{\link{RMsum}} \tab a non-stationary variogram model\cr
}


\bold{Models related to max-stable random fields (tail correlation
  functions)}
\cr
See \link{RMmodelsTailCorrelation}.


\bold{Other covariance models}
\tabular{ll}{
  \command{\link{RMcov}} \tab covariance structure given by a variogram\cr 
  \command{\link{RMfixcov}} \tab User defined covariance structure\cr 
 \command{\link{RMuser}} \tab User defined model \cr
}

\bold{Trend models}
\tabular{ll}{
  \code{\link[=RMS]{Aniso}} \tab for space transformation (not really
  trend, but similar)\cr
  \command{\link{RMcovariate}} \tab spatial covariates\cr
  \command{\link{RMprod}} \tab to model variability of the variance\cr
  \command{\link{RMpolynome}} \tab easy modelling of polynomial trends \cr
  \command{\link{RMtrend}} \tab for explicit trend modelling\cr
  \command{\link{R.models}} \tab for implicit trend modelling\cr
  \command{\link{R.c}} \tab for multivariate trend modelling \cr
}
 

\bold{Auxiliary models}\cr
 See \bold{\link{Auxiliary RMmodels}.}
}

\note{
  \itemize{
    \item
    Note that, instead of the named arguments, a single argument \code{k}
    can be passed. This is possible if all the arguments
    are scalar. Then \code{k} must have a length equal to the number of
    arguments.
    \item
    If an argument equals \code{NULL} the
    argument is not set (but must have a valid name).
    \item
    \code{Aniso} can be given also by \command{\link{RMangle}}
    or any other \command{\link{RMmodel}} instead of a matrix
    \item
    Note also that a completely different possibility exists to define a
    model, namely by a list. This format allows for easy flexible models
    and modifications (and some few more options, as well as some
    abbreviations to the model names, see \command{PrintModelList()}).
    Here, the argument \code{var}, \code{scale},
    \code{Aniso} and \code{proj} must be passed by the model
    \command{\link{RMS}}. 
    For instance,
    \itemize{
      \item
      \code{model <- RMexp(scale=2, var=5)}
      \cr
      is equivalent to
      \cr
      \code{model <- list("RMS", scale=2, var=5, list("RMexp"))}
      \cr
      The latter definition can be also obtained by
      
      \code{print(RMexp(scale=2, var=5))}
      
      \item
      \code{model <- RMnsst(phi=RMgauss(var=7), psi=RMfbm(alpha=1.5),
	scale=2, var=5)}
       \cr
       is equivalent to
        \cr
	\code{model <- list("RMS", scale=2, var=5,} \cr
	\code{list("RMnsst", phi=list("RMS", var=7, list("RMgauss")),} \cr
	\code{psi=list("RMfbm", alpha=1.5))
	)}.
    }
 
    All models have secondary names that stem from 
    \pkg{RandomFields} versions 2 and earlier and
    that can also be used as strings in the list notation.
    See \code{\link{RFgetModelNames}(internal=FALSE)} for
    the full list.     
 }
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
 \item
   Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and
    Strokorb, K. (2015) 
    Analysis, simulation and prediction of multivariate
    random fields with package \pkg{RandomFields}. \emph{
      Journal of Statistical Software}, \bold{63} (8), 1-25,
    url          = \sQuote{http://www.jstatsoft.org/v63/i08/}

    \sQuote{\href{../doc/multivariate_jss.pdf}{multivariate}}, the
    corresponding vignette.

 \item Yaglom, A.M. (1987) \emph{Correlation Theory of Stationary and
 Related Random Functions I, Basic Results.}
 New York: Springer. \item Wackernagel, H. (2003) \emph{Multivariate Geostatistics.} Berlin:
 Springer, 3nd edition.
 }
}

\seealso{\command{\link{RFformula}},
  \link{RM},
  \command{\link{RMmodels}},
  \command{\link{RMmodelsAuxiliary}}.

}

\author{Alexander Malinowski; \martin}
\keyword{spatial}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## a non-stationary field with a sharp boundary
## of the differentiabilities
x <- seq(-0.6, 0.6, len=50)
model <- RMwhittle(nu=0.8 + 1.5 * R.is(R.p(new="isotropic"), "<=", 0.5))
z <- RFsimulate(model=model, x, x, n=4)
plot(z)
\dontshow{FinalizeExample()}}
