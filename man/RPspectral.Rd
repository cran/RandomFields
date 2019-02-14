\name{Spectral}
\alias{Spectral}
\alias{RPspectral}
\title{Spectral turning bands method}

\description{The spectral turning bands method is
 a simulation method for stationary
 Gaussian random fields (Mantoglou and Wilson, 1982). 
 It makes use of
 Bochners's theorem and the corresponding spectral measure
 \eqn{\Xi}{\Xi}
 for a given covariance function \eqn{C(h)}. For \eqn{x \in
 {\bf R}^d}{x in R^d},
 the field \deqn{Y(x)= \sqrt{2} cos(<V,x> + 2 \pi U)}{Y(x)=\sqrt{2}
 cos(<V,x> + 2 pi U)}
 with \eqn{V ~ \Xi } and \eqn{U ~ Ufo((0,1))} is a random field with
 covariance function \eqn{C(h)}.
 A scaled superposition of many independent realizations of \eqn{Y}{Y}
 gives a Gaussian field according to the central limit theorem. For details 
 see Lantuejoul (2002). The standard method
 allows for the simulation of 2-dimensional random
 fields defined on arbitrary points or arbitrary grids. 
 }

\usage{
RPspectral(phi, boxcox, sp_lines, sp_grid, prop_factor, sigma) 
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 specifies the covariance model to be simulated.}
% \item{loggauss}{see \command{\link{RPgauss}}.}
\item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  See \command{\link{RFboxcox}} for details.
 }
 \item{sp_lines}{
 Number of lines used (in total for all additive components of the
 covariance function).

 Default: \code{2500}.
 }
 \item{sp_grid}{Logical.
 The angle of the lines is random if
 \code{grid=FALSE}, 
 and \eqn{k\pi/}\code{sp_lines}
 for \eqn{k}{k} in \code{1:sp_lines},
 otherwise. This argument is only considered
 if the spectral measure, not the density is used.
 
 Default: \code{TRUE}.
 }
 \item{prop_factor}{ % to do: use RRrectangular
   positive real value.
   Sometimes, the spectral density must be sampled by MCMC.
   Let \eqn{p} be the average rejection rate. Then
   the chain is sampled every \eqn{n}th point where
   \eqn{n = |log(p)| *}\code{prop_factor}.

   Default: \code{50}.
 }
 \item{sigma}{real. Considered if the Metropolis
 algorithm is used. It gives the standard deviation of the
 multivariate normal distribution of the proposing
 distribution.
 If \code{sigma}
 is not positive then \code{RandomFields} tries to find a good
 choice for 
 \code{sigma} itself.
 
 Default: \code{0}.
 }
}
 
 

\value{
 \code{RPspectral} returns an object of class
 \code{\link[=RMmodel-class]{RMmodel}}. 

}

\references{
 \itemize{

 \item Lantuejoul, C. (2002)
 \emph{Geostatistical Simulation: Models and Algorithms.}
 Springer.
 \item Mantoglou, A. and J. L. Wilson (1982),
 \emph{The Turning Bands Method for simulation of random fields using
 line generation by a spectral method.}
 Water Resour. Res., 18(5), 1379-1394. 
 }}
\me

\seealso{\link{Gaussian},
 \link{RP},
 \command{\link{RPtbm}}.
}

\keyword{methods}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RPspectral(RMmatern(nu=1))
y <- x <- seq(0,10, len=400)
z <- RFsimulate(model, x, y, n=2)
plot(z)
\dontshow{FinalizeExample()}}
