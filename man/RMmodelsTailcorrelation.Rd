\name{Tail Correlation Funtioncs}
\alias{tcf}
\alias{tail correlation functions}
\alias{Tail correlation functions}
\alias{RMmodelsTailCorrelation}
\title{Covariance models valid for max-stable random fields}
\description{
  This page summarizes the models that can be use
  for tail correlation functions
}

\details{
  The following models are available
  
  \bold{Completely monotone function allowing for arbitray scale}
  \tabular{ll}{
    \command{\link{RMbcw}} \tab Model bridging stationary and
    intrinsically stationary processes for \code{alpha <= 1}
    and \code{beta < 0}\cr
    \command{\link{RMdagum}} \tab Dagum model with \eqn{\beta < \gamma}
    and \eqn{\gamma \le 1}\cr
    \command{\link{RMexp}} \tab exponential model \cr
    \command{\link{RMgencauchy}} \tab generalized Cauchy family with
    \eqn{\alpha \le 1} (and arbitrary \eqn{\beta> 0})\cr
    \command{\link{RMmatern}} \tab Whittle-Matern model with
    \eqn{\nu \le 1/2}\cr
    %multiquadric todo
    %sine power todo
    \command{\link{RMstable}} \tab symmetric stable family or powered
    exponential model with \eqn{\alpha \le 1}\cr
    \command{\link{RMwhittle}} \tab Whittle-Matern model, alternative
    parametrization with \eqn{\nu \le 1/2}\cr 
  }
  
  \bold{Other isotropic models with arbitray scale}
  \tabular{ll}{
    \command{\link{RMnugget}} \tab nugget effect model \cr
  }
  
  \bold{Compactly supported covariance functions}
    \tabular{ll}{
      \command{\link{RMaskey}} \tab Askey's model\cr
      \command{\link{RMcircular}} \tab circular model\cr
       \command{\link{RMconstant}}\tab identically constant \cr
      \command{\link{RMcubic}} \tab cubic model\cr
   \command{\link{RMgengneiting}} \tab Wendland-Gneiting model;
      differentiable models with compact support \cr
   \command{\link{RMgneiting}} \tab differentiable model with compact
    support \cr 
    \command{\link{RMspheric}} \tab spherical model \cr
  }
  
  \bold{Anisotropic models}
  \tabular{ll}{
    none up to now.
  }
  
  \bold{Basic Operators}
  \tabular{ll}{
    \command{\link{RMmult}}, \code{*} \tab product of covariance models \cr
    \command{\link{RMplus}}, \code{+} \tab sum of covariance models or variograms\cr
  }

  \bold{Operators related to process constructions}
    \tabular{ll}{
    \command{\link{RMbernoulli}} \tab correlation of binary fields\cr
    \command{\link{RMbrownresnick}}\tab tcf of a \link{Brown-Resnick} process\cr
    \command{\link{RMschlather}}\tab  tcf of an extremal Gaussian
    process / \link[=RMschlather]{Schlather} process \cr 
    \command{\link{RMm2r}}\tab M2 process with monotone shape function\cr
    \command{\link{RMm3b}}\tab M3 process with balls of random radius\cr
    \command{\link{RMmps}}\tab M3 process with hyperplane polygons\cr
%    \command{\link{}}\tab  \cr
    }

 \bold{See \link{RMmodels} for cartesian models.}
}

\references{
  \itemize{
    \item
    Strokorb, K.,  Ballani, F.,  and Schlather, M. (2015)
    Tail correlation functions of max-stable processes: Construction
    principles, recovery and diversity of some mixing max-stable processes
    with identical TCF. \emph{Extremes}, \bold{18}, 241-271
  }
}

\author{
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}

\seealso{
  \link{coordinate systems},
  \link{RM},
  \command{\link{RMmodels}},
  \command{\link{RMtrafo}}
}

\keyword{spatial}
\keyword{models}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
RFgetModelNames(type="tail")

## an example of a simple model
model <- RMexp(var=1.6, scale=0.5) + RMnugget(var=0) #exponential + nugget
plot(model)

\dontshow{FinalizeExample()}
}