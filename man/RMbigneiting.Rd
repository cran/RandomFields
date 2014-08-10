\name{RMbigneiting}
\alias{RMbigneiting}
\alias{RMbiwendland}
\title{Gneiting-Wendland Covariance Models}
\description{
 \command{\link{RMbigneiting}} is a bivariate stationary isotropic covariance
 model family whose elements 
 are specified by seven parameters.

 Let \deqn{\delta_{ij} = \mu + \gamma_{ij} + 1.}
 Then, 
 \deqn{
 C_{n}(h) = c_{ij} (C_{n, \delta} (h / s_{ij}))_{i,j=1,2}
 }
 and \eqn{ C_{n, \delta} }
 is the generalised Gneiting model
 with parameters \eqn{n} and \eqn{\delta}, see
 \code{\link{RMgengneiting}}, i.e.,
 \deqn{C_{\kappa=0, \delta}(r) = (1-r)^\beta 1_{[0,1]}(r), \qquad \beta=\delta
 + 2\kappa + 1/2;}{ 
 C_{\kappa=0, \delta}(r) = (1 - r)^\beta 1_{[0,1]}(r), \beta=\delta + 2\kappa + 1/2;}
 \deqn{C_{\kappa=1, \delta}(r) = \left(1+\beta r \right)(1-r)^{\beta} 1_{[0,1]}(r),
 \qquad \beta = \delta + 2\kappa + 1/2;}{
 C_{\kappa=1, \delta}(r) = (1+ \beta r)(1-r)^\beta 1_{[0,1]}(r),
 \beta = \delta + 2\kappa + 1/2;} 
 \deqn{C_{\kappa=2, \delta}(r)=\left( 1 + \beta r + \frac{\beta^{2} -
 1}{3}r^{2} \right)(1-r)^{\beta} 1_{[0,1]}(r), \qquad
 \beta=\delta + 2\kappa + 1/2;}{
 C(_{\kappa=2, \delta}(r) = (1 + \beta r + (\beta^2-1) r^(2)/3)(1-r)^\beta
 1_{[0,1]}(r), \beta = \delta + 2\kappa + 1/2;} 
 \deqn{ C_{\kappa=3, \delta}(r)=\left( 1 + \beta r + \frac{(2\beta^{2}-3)}{5} r^{2}+
 \frac{(\beta^2 - 4)\beta}{15} r^{3} \right)(1-r)^\beta 1_{[0,1]}(r),
 \qquad \beta=\delta+2\kappa+1/2.}{
 C_{\kappa=3, \delta}(r) = (1 + \beta r + (2 \beta^2-3 )r^(2)/5+(\beta^2 - 4) \beta
 r^(3)/15)(1-r)^\beta 1_{[0,1]}(r), \beta=\delta + 2\kappa + 1/2.}
 }
\usage{
RMbigneiting(kappa, mu, s, sred12, gamma, cdiag, rhored, c, var, scale, Aniso, proj)
}
\arguments{
 \item{kappa}{argument that chooses between the four different covariance
 models and may take values \eqn{0,\ldots,3}{0,...,3}.
 The model is \eqn{k} times
 differentiable.}
 \item{mu}{\code{mu} has to be greater than or equal to
 \eqn{\frac{d}{2}}{d/2} where \eqn{d}{d} is the (arbitrary)
 dimension of the randomfield.}
 \item{s}{vector of two elements giving the scale of the models on the
 diagonal, i.e., the vector \eqn{(s_{11}, s_{22})}.
 }
 \item{sred12}{value in \eqn{[-1,1]}. The scale on the offdiagonals is
 given by \eqn{s_{12} = s_{21} =}
 \code{sred12 *}
 \eqn{\min\{s_{11},s_{22}\}}{min{s_{11}, s_{22}}}.
 }
 \item{gamma}{a vector of length 3 of numerical values; each entry is
 positive.
 The vector \code{gamma} equals
 \eqn{(\gamma_{11},\gamma_{21},\gamma_{22})}.
 Note that \eqn{\gamma_{12} =\gamma_{21}}.
 }
 \item{cdiag}{a vector of length 2 of numerical values; each entry
 positive; the vector \eqn{(c_{11},c_{22})}}
 \item{c}{a vector of length 3 of numerical values;
 the vector \eqn{(c_{11}, c_{21}, c_{22})}.
 Note that \eqn{c_{12}= c_{21}}.

 Either 
 \code{rhored} and \code{cdiag} or \code{c} must be given.
 }
 \item{rhored}{value in \eqn{[-1,1]}.
 See
 also the Details for the corresponding value of \eqn{c_{12}=c_{21}}.
 }
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 A sufficient condition for the
 constant \eqn{c_{ij}} is
 \deqn{c_{12} = \rho_{\rm red} \cdot m \cdot \left(c_{11} c_{22}
 \prod_{i,j=1,2}
 \left(\frac{\Gamma(\gamma_{ij} + \mu + 2\kappa + 5/2)}{b_{ij}^{\nu_{ij} +
 2\kappa + 1} \Gamma(1 + \gamma_{ij}) \Gamma(\mu + 2\kappa + 3/2)}
 \right)^{(-1)^{i+j}}
 \right)^{1/2}
 }{
 c_{ij} = \rho_r m (c_{11} c_{22})^{1/2}
 }
 where \eqn{\rho_{\rm red} \in [-1,1]}{\rho_r in [-1,1]}.

 The constant \eqn{m} in the formula above is obtained as follows:
 \deqn{m = \min\{1, m_{-1}, m_{+1}\}}{m = min\{1, m_{-1}, m_{+1}\}}
 Let
 \deqn{a = 2 \gamma_{12} - \gamma_{11} -\gamma_{22}}
 \deqn{b = -2 \gamma_{12} (s_{11} + s_{22}) + \gamma_{11} (s_{12} +
 s_{22}) + \gamma_{22} (s_{12} + s_{11})}
 \deqn{e = 2 \gamma_{12} s_{11}s_{22} - \gamma_{11}s_{12}s_{22} -
 \gamma_{22}s_{12}s_{11}}
 \deqn{d = b^2 - 4ae}
 \deqn{t_j =\frac{- b + j \sqrt d}{2 a} }{t_j =(-b + j \sqrt d) / (2 a) }
 If \eqn{d \ge0} and \eqn{t_j \not\in (0, s_{12})}{t_j in (0, s_{12})^c} then \eqn{m_j=\infty} else
 \deqn{
 m_j =
 \frac{(1 - t_j/s_{11})^{\gamma_{11}}(1 -
 t_j/s_{22})^{\gamma_{22}}}{(1 - t_j/s_{12})^{2 \gamma_{11}}
 }{
 m_j = (1 - t_j/s_{11})^{\gamma_{11}} (1 -
 t_j/s_{22})^{\gamma_{22}} / (1 - t_j/s_{12})^{2 \gamma_{11}}
 }
 }

 In the function \command{\link{RMbigneiting}}, either \code{c} is
 passed, then the above condition is checked, or \code{rhored} is passed
 then \eqn{c_{12}} is calculated by the above formula.
 
}
\value{
 \command{\link{RMgengneiting}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{
 \itemize{
 \item Bevilacqua, M., Daley, D.J., Porcu, E., Schlather, M. (2012)
 Classes of compactly supported correlation functions for multivariate
 random fields. Arxiv.

 
 \item Gneiting, T. (1999)
 Correlation functions for atmospherical data analysis.
 \emph{Q. J. Roy. Meteor. Soc} Part A \bold{125}, 2449-2464.

 \item Wendland, H. (2005) \emph{Scattered Data Approximation.}
 {Cambridge Monogr. Appl. Comput. Math.}
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMaskey}},
 \command{\link{RMbiwm}},
 \command{\link{RMgengneiting}},
 \command{\link{RMgneiting}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
%# gamma is mainly a scale effect
model <- RMbigneiting(kappa=2, mu=0.5, gamma=c(0, 3, 6), rhored=1)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}