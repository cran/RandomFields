\name{RMflatpower}
\alias{RMflatpower}
\title{Variogram Model Similar to Fractal Brownian Motion}
\description{
 \command{\link{RMflatpower}} is an intrinsically stationary isotropic variogram
 model. The corresponding centered semi-variogram only depends on the
 distance \eqn{r \ge 0}{r \ge 0} between two points and is given by
 \deqn{\gamma(r) = r^2 / ( 1 + r^2)^\alpha}
 where \eqn{\alpha \in (0,1]}{0 < \alpha \le 1}.\cr
 
  For related models see \command{\link{RMgenfbm}}.
}
\usage{
RMflatpower(alpha, var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{numeric in \eqn{(0,1]}; refers to the fractal dimension of the
 process}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 variogram remains unmodified.}
}
\details{
  The model is always smooth at the origin. 
  
  The parameter \eqn{\alpha}{\alpha} only gives the tail behaviour and
  satisfies \eqn{\alpha \in (0,1]}{0
 < \alpha \le 1}.

The variogram is unbounded and belongs to a non-stationary process with
stationary increments. 



}
\value{
 \command{\link{RMflatpower}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
 \itemize{
 \item Oesting, M., Schlather, M., and Friederichs, P. (2014)
 Conditional Modelling of Extreme Wind Gusts
 by Bivariate {Brown-Resnick} Processes
 \emph{arxiv} \bold{1312.4584}.
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMgenfbm}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

model <- RMflatpower(alpha=0.5)
x <- seq(0, 10, 0.1)
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}