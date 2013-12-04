\name{RMbiwm}
\alias{RMbiwm}
\title{Full Bivariate Whittle Matern Model}
\description{
 \command{\link{RMbiwm}} is a bivariate stationary isotropic covariance model
 whose corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between
 two points and is given for \eqn{i,j \in \{1,2\}}{i,j = 1,2} by
 \deqn{C_{ij}(r)=c_{ij} W_{\nu_{ij}}(r/s_{ij}).}{C_{ij}(r)=c_{ij} W_{\nu_{ij}}(r/s_{ij}).}
 Here \eqn{W_\nu} is the covariance of the
 \command{\link{RMwhittle}} model.
 For constraints on the constants see details.
}

\usage{
RMbiwm(nudiag, nured12, nu, s, cdiag, rhored, c, notinvnu, var,
 scale, Aniso, proj)
}
\arguments{
 \item{nudiag}{a vector of length 2 of numerical values; each entry
 positive; the vector \eqn{(\nu_{11},\nu_{22})}}
 \item{nured12}{a numerical value in the interval \eqn{[1,\infty)};
 \eqn{\nu_{21}} is calculated as \eqn{0.5 (\nu_{11} + \nu_{22})*\nu_{red}}.}
 \item{nu}{alternative to \code{nudiag} and \code{nured12}:
 a vector of length 3 of numerical values; each entry
 positive; the vector \eqn{(\nu_{11},\nu_{21},\nu_{22})}.
 Either 
 \code{nured} and \code{nudiag} or \code{nu} must be given.}
 \item{s}{a vector of length 2 of numerical values; each entry
 positive; the vector \eqn{(s_{11},s_{22})}}
 \item{cdiag}{a vector of length 2 of numerical values; each entry
 positive; the vector \eqn{(c_{11},c_{22})}}
 \item{rhored}{a numerical value; in the interval \eqn{[-1,1]}. See
 also the Details for the corresponding value of \eqn{c_{12}=c_{21}}.
 }
 \item{c}{a vector of
 length 3 of numerical values; 
 the vector \eqn{(c_{11},c_{21}, c_{22})}. Either 
 \code{rhored} and \code{cdiag} or \code{c} must be given.}
 \item{notinvnu}{logical or \code{NULL}.
 If not given (default) then the formula of the 
 (\command{\link{RMwhittle}}) model applies.
 If logical then the formula for the \command{\link{RMmatern}} model
 applies. See there for details.
 }
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{
 Constraints on the constants: 
 For the diagonal elements we have
 \deqn{\nu_{ii}, s_{ii}, c_{ii} > 0.}
 For the offdiagonal elements we have 
 \deqn{s_{12}=s_{21} > 0,}
 \deqn{\nu_{12} =\nu_{21} = 0.5 (\nu_{11} + \nu{22}) / \nu_{red}}
 for some constant \eqn{\nu_{red} \in (0,1]} and
 \deqn{c_{12} =c_{21} = \rho_{red} \sqrt{f m c_{11} c_{22}}}
 for some constant \eqn{\rho_{red}} in \eqn{[-1,1]}.
 
 The constants \eqn{f} and \eqn{m} in the last equation are given as follows:
 \deqn{f = (\Gamma(\nu_{11} + d/2) \Gamma(\nu_{22} + d/2)) / (\Gamma(\nu_{11}) \Gamma(\nu_{22})) * (\Gamma(\nu_{12}) / \Gamma(\nu_{12}+d/2))^2 * ( s_{12}^{2*\nu_{12}} / (s_{11}^{\nu_{11}} s_{22}^{\nu_{22}}) )^2}
 where \eqn{\Gamma} is the Gamma function and \eqn{d} is the dimension
 of the space.
 The constant \eqn{m} is 
 the infimum of the function \eqn{g} on \eqn{[0,\infty)} where
 \deqn{g(t) = (1/s_{12}^2 +t^2)^{2\nu_{12} + d} (1/s_{11}^2 + t^2)^{-\nu_{11}-d/2} (1/s_{22}^2 + t^2)^{-\nu_{22}-d/2}} 
 (cf. Gneiting, T., Kleiber, W., Schlather, M. (2010), Full 
 Bivariate Matern Model (Section 2.2))
 % For an alternative model see also \command{\link{RMbiwm}}. 
}
\value{
 \command{\link{RMbiwm}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Gneiting, T., Kleiber, W., Schlather, M. (2010)
 Matern covariance functions for multivariate random fields
 \emph{JASA}
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMparswm}},
 \command{\link{RMwhittle}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
set.seed(0)
x <- y <- seq(-10, 10, if (interactive()) 0.2 else 5)
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
                s=c(1, 1, 2))
simu <- RFsimulate(model, x, y, grid=TRUE)
plot(simu)
}