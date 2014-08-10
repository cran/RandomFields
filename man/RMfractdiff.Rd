\name{RMfractdiff}
\alias{RMfractdiff}
\title{Fractionally Differenced Process Model}
\description{
 \command{\link{RMfractdiff}} is a stationary isotropic covariance model.
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between two points and is given for integers
 \eqn{r \in {\bf N}}{r} by
 \deqn{C(r) = (-1)^r \frac{ \Gamma(1-a/2)^2 }{ \Gamma(1-a/2+r) \Gamma(1-a/2-r) } r \in {\bf N}}{C(r) = (-1)^r \Gamma(1-a/2)^2 / (\Gamma(1-a/2+r) \Gamma(1-a/2-r))}
 and otherwise linearly interpolated. Here \eqn{a \in [-1,1)}{-1 \le a < 1},
 \eqn{\Gamma}{\Gamma} denotes the gamma function.
 It can only be used for one-dimensional random fields.
}
\usage{
RMfractdiff(a, var, scale, Aniso, proj)
}
\arguments{
 \item{a}{ \eqn{-1 \le a < 1}}
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}
\details{The model is only valid for dimension \eqn{d = 1}{d = 1 }.
 It stems from time series modelling where the grid locations are
 multiples of the scale parameter.
}
\value{
 \command{\link{RMfractdiff}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}
%\references{
 % reference missing!
 %\itemize{
 % \item
 %}
%}

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
model <- RMfractdiff(0.5, scale=0.2)
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model)
plot(RFsimulate(model, x=x))
\dontshow{FinalizeExample()}
}