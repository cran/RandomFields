\name{RMwhittlematern}
\alias{RMwhittle}
\alias{RMkbessel}
\alias{RMmatern}
\alias{Sobolev}
\title{Whittle-Matern Covariance Model}
\description{
 \command{\link{RMmatern}} is a stationary isotropic covariance model
 belonging to the Matern family. 
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0}
 between two points.

 The Whittle model is given by
 \deqn{C(r)=W_{\nu}(r)=2^{1- \nu}
 \Gamma(\nu)^{-1}r^{\nu}K_{\nu}(r)}{C(r)=W_{\nu}(r)=2^{1- \nu}
 \Gamma(\nu)^{-1}r^{\nu}K_{\nu}(r)} 
 where \eqn{\nu > 0}{\nu > 0} and \eqn{K_\nu}{K_\nu} is the modified
 Bessel function of second kind.

 

 The Matern model is given by
 \deqn{C(r) = \frac{2^{1-\nu}}{\Gamma(\nu)} (\sqrt{2\nu}r)^\nu
 K_\nu(\sqrt{2\nu}r)}{C(r) = 2^{1- \nu} \Gamma(\nu)^{-1} (\sqrt{2\nu}
 r)^\nu K_\nu(\sqrt{2\nu} r)} 

 }
\usage{
RMwhittle(nu, notinvnu, var, scale, Aniso, proj)

RMmatern(nu, notinvnu, var, scale, Aniso, proj)
}
\arguments{
 \item{nu}{a numerical value called \dQuote{smoothness parameter};
 should be greater than 0.} 
 \item{notinvnu}{logical, if not given the model is defined as above.
 (default). This argument should not be set
 by users. See the Notes.}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{

 \command{\link{RMwhittle}} and \command{\link{RMmatern}}
 are two alternative parametrizations of the same covariance function.
 The Matern model should be preferred as this model seperates the
 effects of scaling parameter and the shape parameter. 

 This is the normal scale mixture model of choice if the smoothness of a
random field is to be parametrized: the sample paths of a Gaussian
random field with this covariance structure are \eqn{m}{m} times
differentiable if and only if \eqn{\nu > m}{\nu > m} (see Gelfand et
al., 2010, p. 24).

Furthermore, the fractal dimension (see also \command{\link{RFfractaldim}})
\emph{D} of the Gaussian sample paths
is determined by \eqn{\nu}{\nu}: we have
\deqn{D = d + 1 - \nu, \nu \in (0,1)}{D = d + 1 - \nu, 0 < \nu < 1}
and \eqn{D = d}{D = d} for \eqn{\nu > 1}{\nu > 1} where \eqn{d}{d} is
the dimension of the random field (see Stein, 1999, p. 32).

If \eqn{\nu=0.5}{\nu=0.5} the Matern model equals \command{\link{RMexp}}.

For \eqn{\nu}{\nu} tending to \eqn{\infty}{\infty} a rescaled Gaussian
model \command{\link{RMgauss}} appears as limit of the Matern model. 

 For generalisations see section \sQuote{seealso}.
}

\note{
 The Matern model called by \eqn{C(r \sqrt(2))}{C(r \sqrt(2))} equals
 the Handcock-Wallis (1994) parametrisation. 

 The model allows further to be reparameterized by substituting
 \eqn{\nu}{\nu} for \eqn{\nu^{-1}}{\nu^{-1}} setting the argument
 \code{invnu=TRUE}. Note that the inversion of \eqn{\nu} does not really
 make sense for the Whittle model. Due to this fact, if the argument
 \code{invnu} is given, the Whittle model changes its definition and
 becomes identical to the Matern model.
}

\value{
 The function return an object of class \code{\link[=RMmodel-class]{RMmodel}}

}
\references{
 \itemize{
 \item Chiles,
 J.-P. and Delfiner, P. (1999)
 \emph{Geostatistics. Modeling Spatial Uncertainty.}
 New York: Wiley.
 
 \item Gelfand, A. E., Diggle, P., Fuentes, M. and Guttorp,
 P. (eds.) (2010) \emph{Handbook of Spatial Statistics.}
 Boca Raton: Chapman & Hall/CRL.
 
 \item Guttorp, P. and Gneiting, T. (2006) Studies in the
 history of probability and statistics. XLIX. On the Matern
 correlation family. \emph{Biometrika} \bold{93}, 989--995.
 
 \item Handcock, M. S. and Wallis, J. R. (1994) An approach to
 statistical spatio-temporal modeling of meteorological fields.
 \emph{JASA} \bold{89}, 368--378.

 \item Stein, M. L. (1999) \emph{Interpolation of Spatial Data --
   Some Theory for Kriging.} New York: Springer.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \itemize{
 \item \command{\link{RMexp}}, \command{\link{RMgauss}} for special
 cases of the model (for \eqn{\nu=0.5}{\nu=0.5} and
 \eqn{\nu=\infty}{\nu=\infty}, respectively)

 \item \command{\link{RMhyperbolic}} for a univariate
 generalization

 \item \command{\link{RMbiwm}} for a multivariate generalization
 
 \item \command{\link{RMnonstwm}}, \command{\link{RMstein}} for anisotropic (space-time) generalizations
 
% \item \command{\link{}} for 
 
 \item \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}} for general use.
 }

}


\keyword{spatial}
\keyword{models}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
x <- seq(0, 1, len=if (interactive()) 100 else 3)
model <- RMwhittle(nu=1, Aniso=matrix(nc=2, c(1.5, 3, -3, 4)))
plot(model, dim=2, xlim=c(-1,1))
z <- RFsimulate(model=model, x, x)
plot(z)
\dontshow{FinalizeExample()}
}
