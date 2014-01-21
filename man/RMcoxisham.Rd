\name{RMcoxisham}
\alias{RMcoxisham}
\title{Cox Isham Covariance Model}
\description{
 \command{\link{RMcoxisham}} is a stationary covariance model
 which depends on a univariate stationary isotropic covariance model
 \eqn{C_0}, which is a normal scale mixture.
 
 The corresponding covariance function only depends on the difference
 \eqn{(h,t) \in {\bf R}^{d+1}={\bf R}^d\times{\bf R}}{(h,t)} between two points in \eqn{d+1}-dimensional space and is given by
 \deqn{C(h,t)=|E + t^\beta D|^{-1/2} C_0([(h - t \mu)^T (E + t^\beta
 D)^{-1} (h - t \mu)]^{1/2})}
 Here \eqn{\mu \in {\bf R}^d}{\mu} is a vector in
 \eqn{d}-dimensional space;
 \eqn{E} is the \eqn{d \times d}{d x d}-identity matrix and \eqn{D} is
 a \eqn{d \times d}{d x d}-correlation matrix with \eqn{|D| > 0}.
 The parameter \eqn{\beta} is in \eqn{(0,2]}.
 Currently, the implementation is done only for \eqn{d=2}. 
}

\usage{
RMcoxisham(phi,mu,D,beta,var, scale, Aniso, proj)
}
\arguments{
 \item{phi}{a univariate stationary isotropic covariance model for random fields
 on \eqn{d}-dimensional space, which is moreover a normal scale
 mixture, that means an
 \command{\link{RMmodel}} whose \code{normalmix} equals
 TRUE and whose \code{maxdim} is at least 2.}
 \item{mu}{a vector in \eqn{d}-dimensional space}
 \item{D}{a \eqn{d \times d}{d x d}-correlation matrix with \eqn{|D| >
 0}}
 \item{beta}{a parameter in the interval \eqn{(0,2]}, default value is 2 }
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 This model stems from a rainfall model (cf. Cox, D.R., Isham,
 V.S. (1988)) and equals the following expectation
 \deqn{C(h,t)=\bold{E}_V C_0(h-Vt)}
 where the random wind speed vector \eqn{V} follows a \eqn{d}-variate
 normal distribution with expectation \eqn{mu} and covariance matrix \eqn{D/2}.
 (cf. See Schlather, M. (2010), Example 9).
}
\value{
 \command{\link{RMcoxisham}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Cox, D.R., Isham, V.S. (1988)
 A simple spatial-temporal model of rainfall.
 \emph{Proc. R. Soc. Lond. A}, \bold{415}, 317-328.

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
RFoptions(seed=0)
model <- RMcoxisham(RMgauss(), mu=1, D=1)
x <- seq(0, 10, if (interactive()) 0.3 else 1) 
plot(model, ylim=c(0, 1), dim=2)
plot(RFsimulate(model, x=x, y=x))
\dontshow{RFoptions(seed=NA)}
}