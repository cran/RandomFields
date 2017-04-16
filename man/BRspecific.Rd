\name{Brown-Resnick-Specific}
\alias{RPbrmixed}
\alias{RPbrorig}
\alias{RPbrshifted}
\alias{BRmethods}
\title{Simulation methods for Brown-Resnick processes}



\description{
 These models define the particular way to simulate Brown-Resnick
 processes
 }

\usage{
RPbrmixed(phi, tcf, xi, mu, s, meshsize, vertnumber, optim_mixed,
          optim_mixed_tol,lambda, areamat, variobound) 

RPbrorig(phi, tcf, xi, mu, s)

RPbrshifted(phi, tcf, xi, mu, s)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 specifies the covariance model to be simulated.}

 \item{tcf}{the extremal correlation function; either \code{phi} or
 \code{tcf} must be given.}

 \item{xi, mu, s}{the shape parameter, the location parameter and the
 scale parameter, respectively, of the generalized extreme value
 distribution. See Details.}

\item{lambda}{positive constant factor in the intensity of the Poisson
  point processused in the M3 representation, cf. Thm. 6 and Remark 7
  in Oesting et. al (2012); can be estimated by setting
  \code{optim_mixed} if unknown. Default value is 1.}

\item{areamat}{vector or matrix of values in \eqn{[0,1]} with odd
  length (odd number of rows and columns, respectively). Each value
  represents the portion of processes whose maximum is located at a
  specific location on a grid taken into account for the simulation
  of the shape function in the M3 representation. The center of
  \code{areamat} represents the value for the origin, the other entries
  belong to the corresponding locations on a 1D or 2D grid.
  \code{areamat} can be used for dimensions 1 and 2 only; can be
  optimized by setting \code{optim_mixed} if unknown. Default value is
  1.}
   
\item{meshsize, vertnumber, optim_mixed,
  optim_mixed_tol, variobound}{further arguments
  for simulation via the mixed moving maxima (M3) representation; see
  \code{\link{RFoptions}}}
}

\details{
  The argument \code{xi} is always a number, i.e. \eqn{\xi} is constant
  in space. In contrast, \eqn{\mu} and \eqn{s} might be constant
  numerical value or given a \code{\link{RMmodel}}, in particular by a
  \code{\link{RMtrend}} model. 

  The functions \code{RPbrorig}, \code{RPbrshifted} and \code{RPbrmixed}
  simulate a Brown-Resnick process, which is defined by
  \deqn{Z(x) = \max_{i=1}^\infty X_i \exp(W_i(x) - \gamma),
  }{Z(x) = max_{i=1, 2, ...} X_i * exp(W_i(x) - gamma),}
 where the \eqn{X_i} are the points of a Poisson point process on the
 positive real half-axis with intensity \eqn{x^{-2} dx}{1/x^2 dx},
 \eqn{W_i \sim W}{W_i ~ Y} are iid centered Gaussian processes with
 stationary increments and variogram \eqn{\gamma}{gamma} given by
 \code{model}. The functions correspond to the following ways of
 simulation:
 \describe{
  \item{\code{RPbrorig}}{simulation via using the original definition
  (method 0 in Oesting et al., 2012)}

  \item{\code{RPbrshifted}}{simulation using a random shift (similar to
  method 1 and 2)}
  
  \item{\code{RPbrmixed}}{simulation using M3 representation (method
  4)}
 }
}

\value{
 The functions return an object of class
 \code{\link[=RMmodel-class]{RMmodel}}
}

\references{
 \itemize{
   \item  Oesting, M., Kabluchko, Z. and Schlather M. (2012)
   Simulation of {B}rown-{R}esnick Processes, \emph{Extremes},
   \bold{15}, 89-107.
 }}

 \note{Advanced options for \code{RPbroriginal} and \code{RPbrshifted}
   are \code{maxpoints} and \code{max_gauss}, see \command{\link{RFoptions}}.}
       
  
\author{Marco Oesting, \email{oesting@math.uni-mannheim.de},
        Martin Schlather, \email{schlather@math.uni-mannheim.de}
        \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}

% TO DO
## currently does not work

\dontshow{\dontrun{
model <- RPbrshifted(RMfbm(alpha=1.5), xi=0)
x <- 0:10
z <- RFsimulate(model=model, x=x, y=x, n=4)
plot(z)
}}

\dontshow{FinalizeExample()}
}

\seealso{
 \command{\link{RPbrownresnick}},
 \command{\link{RMmodel}},
 \command{\link{RPgauss}},
 \command{\link{maxstable}},
 \command{\link{maxstableAdvanced}}
}

\keyword{methods}
