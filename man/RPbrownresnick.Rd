\name{BrownResnick}
\alias{RPbrownresnick}
\alias{Brown-Resnick}
\alias{Brown-Resnick process}

\title{Brown-Resnick process}
\description{ 
  \command{RPbrownresnick} defines a Brown-Resnick process.
}

\usage{
  RPbrownresnick(phi, tcf, xi, mu, s)
}

\arguments{
  \item{phi}{specifies the covariance model or variogram, see
    \link{RMmodel} and \link{RMmodelsAdvanced}.
  }

  \item{tcf}{the extremal correlation function; either \code{phi} or
     \code{tcf} must be given.}

   \item{xi, mu, s}{the extreme value index, the location parameter and the
     scale parameter, respectively, of the generalized extreme value
     distribution. See Details.}
}

\details{
  The extreme value index \code{xi} is always a number, i.e. \eqn{\xi}
  is constant 
  in space. In contrast, \eqn{\mu} and \eqn{s} might be constant
  numerical value or given a \code{\link{RMmodel}}, in particular by a
  \code{\link{RMtrend}} model. The default values of \eqn{mu} and \eqn{s}
  are \eqn{1} and \eqn{z\xi}, respectively.

 The functions \code{RPbrorig}, \code{RPbrshifted} and \code{RPbrmixed}
  perform the simulation of a Brown-Resnick process, which is defined by
  \deqn{Z(x) = \max_{i=1}^\infty X_i \exp(W_i(x) - \gamma^2),
  }{Z(x) = max_{i=1, 2, ...} X_i * exp(W_i(x) - gamma^2),}
  where the \eqn{X_i} are the points of a Poisson point process on the
  positive real half-axis with intensity \eqn{x^{-2} dx}{1/x^2 dx},
  \eqn{W_i \sim W}{W_i ~ Y} are iid centered Gaussian processes with
  stationary increments and variogram \eqn{\gamma}{gamma} given by
  \code{model}.
   
 For simulation, internally, one of the methods
 \command{\link{RPbrorig}}, \command{\link{RPbrshifted}} and
 \command{\link{RPbrmixed}} is chosen automatically.
 }

 
\author{Marco Oesting, \email{oesting@math.uni-mannheim.de},
        Martin Schlather, \email{schlather@math.uni-mannheim.de}
        \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\references{
  \itemize{
    \item
    Brown, B.M. and Resnick, S.I. (1977).
    Extreme values of independent stochastic
    processes. \emph{J. Appl. Probab.} \bold{14}, 732-739.


    \item
    Buishand, T., de Haan , L. and Zhou, C. (2008).
    On spatial extremes: With application to
    a rainfall problem. \emph{Ann. Appl. Stat.} \bold{2}, 624-642.


    \item Kabluchko, Z., Schlather, M. and de Haan, L (2009)
    Stationary max-stable random fields associated to negative
    definite functions \emph{Ann. Probab.} \bold{37}, 2042-2065.
    
    \item  Oesting, M., Kabluchko, Z. and Schlather M. (2012)
    Simulation of {B}rown-{R}esnick Processes, \emph{Extremes},
    \bold{15}, 89-107.
  }
}

% TO DO: ueberall diese notes einfuegen
\note{Advanced options
   are \code{maxpoints} and \code{max_gauss}, see
   \command{\link{RFoptions}}.
   
  Further advanced options related to the
  simulation methods \command{\link{RPbrorig}},
  \command{\link{RPbrshifted}} and \command{\link{RPbrmixed}} can be
  found in the paragraph \sQuote{Specific method options for Brown-Resnick
    Fields} in \command{\link{RFoptions}}.}


\seealso{
 \command{\link{RPbrorig}},
 \command{\link{RPbrshifted}},
 \command{\link{RPbrmixed}},
 \command{\link{RMmodel}},
 \command{\link{RPgauss}},
 \command{\link{maxstable}},
 \command{\link{maxstableAdvanced}}
}

\keyword{spatial}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}
\dontshow{\dontrun{
model <- ~ RPbrownresnick(RMfbm(alpha=A), xi=0)
x <- seq(0, 10, 0.2)
z <- RFsimulate(model=model, x, x, n=4, A=0.9) # about 1 min on a fast machine
plot(z)

z <- RFsimulate(model, x=x, n=4, A=1.9)
plot(z)


## basic model in Buishand, de Haan, Zhou (2008)
model <- RMfbm(proj=1, alpha=1, var=0.5) + RMfbm(proj=2, alpha=1, var=0.5) 
x <- seq(0, 5, 0.05)
z <- RFsimulate(RPbrownresnick(model, xi=0), x, x, every=1000)
plot(z)


}}
## for some more sophisticated models see 'maxstamableAdvanced'

\dontshow{FinalizeExample()}
}