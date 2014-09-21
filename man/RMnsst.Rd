\name{RMnsst}
\alias{RMnsst}
\title{Non-Separable Space-Time model}
\description{
 \command{\link{RMnsst}} is a univariate stationary spaceisotropic covariance model
 whose corresponding covariance is given by
 \deqn{C(h,u)= (\psi(u)+1)^{-\delta/2} \phi(h /\sqrt(\psi(u) +1))}
}
\usage{
RMnsst(phi, psi, delta, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{is normal mixture \command{\link{RMmodel}},  cf.\cr
    \code{RFgetModelNames(monotone="normal mixture")}
  }
 \item{psi}{is a variogram \command{\link{RMmodel}}.}
 \item{delta}{a numerical value; must be greater than or equal to the spatial dimension of the field}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model is used for space-time modelling where the spatial
 component is isotropic.
}
\value{
 \command{\link{RMnsst}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Gneiting, T. (1997) Normal scale mixtures and dual probability
 densitites, \emph{J. Stat. Comput. Simul.} \bold{59}, 375-384.
 
 \item Gneiting, T. (2002) Nonseparable, stationary covariance
 functions for space-time data, \emph{JASA} \bold{97}, 590-600.
 
 \item Gneiting, T. and Schlather, M. (2001)
 Space-time covariance models.
 In El-Shaarawi, A.H. and Piegorsch, W.W.:
 \emph{The Encyclopedia of Environmetrics.} Chichester: Wiley.

% \item Zastavnyi, V. and Porcu, E. (2011)
% Caracterization theorems for the Gneiting class space-time
% covariances.
% \emph{Bernoulli}, \bold{??}.

 \item Schlather, M. (2010)
 On some covariance models based on normal scale mixtures.
 \emph{Bernoulli}, \bold{16}, 780-797.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMnsst(phi=RMgauss(), psi=RMfbm(alpha=1), delta=2)
x <- seq(0, 10, if (interactive()) 0.25 else 1) 
plot(model, dim=2)
plot(RFsimulate(model, x=x, y=x))
\dontshow{FinalizeExample()}
}
