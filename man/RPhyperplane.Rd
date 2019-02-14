\name{Hyperplane}
\alias{Hyperplane}
\alias{Hyperplanes}
\alias{RPhyperplane}
\title{Hyperplane method}
\description{

 
 The Hyperplane method is a simulation method for stationary, isotropic
 random fields with exponential covariance function. It is
 based on a tessellation of the space by
 hyperplanes. Each cell takes a spatially constant value
 of an i.i.d. random variable. The superposition of several
 such random fields yields approximatively a Gaussian random field.
 
}

\usage{
RPhyperplane(phi, boxcox, superpos, maxlines, mar_distr, mar_param ,additive)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 specifies the covariance function to be simulated;
 only exponential covariance functions and scale mixtures
 of it are allowed.
 }
% \item{loggauss}{see  \command{\link{RPgauss}}.}
\item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  See \command{\link{RFboxcox}} for details.
 }
 \item{superpos}{integer. number of superposed hyperplane tessellations.
 Default: \code{300}.}
 
 \item{maxlines}{integer.
 Maximum number of allowed lines.
 
 Default: \code{1000}.
 }
 
 \item{mar_distr}{integer.
 code for the marginal distribution used in the
 simulation:
 \describe{
	\item{\code{0}}{uniform distribution}
	\item{\code{1}}{Frechet distribution with form parameter
	 \code{mar_param}}
	\item{\code{2}}{Bernoulli distribution (Binomial with \eqn{n=1}) with
	 argument \code{mar_param}}
 }
 This argument should not be changed yet.
 
 Default: \code{0}.
 }
 \item{mar_param}{Argument used for the marginal
 distribution. The argument should not be changed yet.
 
 Default: \code{NA}.
 }

 \item{additive}{logical.
   If \code{TRUE} independent realizations are added, else
   the maximum is taken.

   Default: \code{TRUE}.
 }
 }


\value{
 \code{RPhyperplane} returns an object of class
 \code{\link[=RMmodel-class]{RMmodel}}. 

}

\references{
 \itemize{
 \item Lantuejoul, C. (2002)
 \emph{Geostatistical Simulation: Models and Algorithms.}
 Springer.
 }
}

\me

\seealso{ \link{Gaussian},
  \link{RP}.
}

\keyword{methods}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RPhyperplane(RMexp(s=2), superpos=1)
x <- seq(0, 3, 0.04)
z <- RFsimulate(x=x, x, model=model, n=1)
plot(z)

\dontshow{FinalizeExample()}}
