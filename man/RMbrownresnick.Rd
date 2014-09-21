\name{RMbrownresnick}
\alias{RMbrownresnick}
\alias{powered error function}
\alias{error function model}
\title{Tail correlation function of the Brown-Resnick process}
\description{
  \command{RMbrownresnick} defines
  the tail correlation function of the Brown-Resnick process.
  \deqn{C(h) = 2 - 2\Phi(\sqrt{\gamma(h)} / 2)}
  where \eqn{\phi} is the standard normal distribution function
  and \eqn{\gamma} is the \bold{semi-}variogram.
}
\usage{
RMbrownresnick(phi, var, scale, Aniso, proj);
}
\arguments{
 \item{phi}{variogram of class \code{\link[=RMmodel-class]{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\value{
 object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\note{
  \bold{In the paper Kabluchko et al (2009)
    the variogram instead of the semi-variogram is considered, so the
    formulae differ slightly.}

  \bold{In Version 3.0.33 a typo has been corrected.}

  Here, a definition is used that is consistend with the rest of the package.
}

\details{
 For a given \command{\link{RMmodel}} the function
 \code{\link{RMbrownresnick}(\link{RMmodel}())} 'returns' the tail correlation
 function of a Brown-Resnick process with variogram \command{\link{RMmodel}}.
}
\references{
  \itemize{
    \item Kabluchko, Z., Schlather, M. & de Haan, L (2009)
    Stationary max-stable random fields associated to negative definite
    functions 
    \emph{Ann. Probab.} \bold{37}, 2042-2065.

    \item Strokorb, K., Ballani, F., and  Schlather, M. (2014)
    Tail correlation functions of max-stable processes: Construction
    principles, recovery and diversity of some mixing max-stable
    processes with identical TCF. 
    \emph{Extremes}, \bold{} Submitted.    
  }
}
\seealso{
 \command{\link{RFsimulate}},
 \command{\link{RMm2r}}, \command{\link{RMm3b}}, \command{\link{RMmps}},
 \command{\link{RMmodel}}.
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}
}
\keyword{spatial}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

#plot covariance model of type RMbrownresnick
RMmodel <- RMfbm(alpha=1.5, scale=0.2)
plot(RMbrownresnick(RMmodel))

#simulate and plot corresponding Gaussian random field
x <- if (interactive()) seq(-5, 5, 0.05) else seq(-5, 5, 5)
z <- RFsimulate(RMbrownresnick(RMmodel), x=x, y=x)
plot(z)
\dontshow{FinalizeExample()}
}
