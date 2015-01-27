\name{RMgennsst}
\alias{RMgennsst}
\title{Non-Separable Space-Time model}
\description{
 \command{\link{RMgennsst}} is a univariate stationary spaceisotropic covariance model
 whose corresponding covariance is given by
 \deqn{C(h,u)= (\psi(u)+1)^{-\delta/2} \phi(h /\sqrt(\psi(u) +1))}
}
\usage{
RMgennsst(phi, psi, var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{is normal mixture \command{\link{RMmodel}},  cf.\cr
    \code{RFgetModelNames(monotone="normal mixture")}
  }
 \item{psi}{is a variogram \command{\link{RMmodel}}.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model is used for space-time modelling where the spatial
 component is isotropic.
}
\value{
 \command{\link{RMgennsst}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{

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
 \command{\link{RMnsst}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}


\examples{
% to do
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

% RMgennsst( RMwhittle(nu=4, s=3), RMconstant(diag(3))


\dontshow{FinalizeExample()}
}
