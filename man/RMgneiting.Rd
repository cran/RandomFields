\name{RMgneiting}
\alias{RMgneiting}
\title{Gneiting Covariance Model}
\description{
 \command{\link{RMgneiting}} is a stationary isotropic covariance model
 which is only valid up to dimension 3, or 5 (see the argument \code{orig}).
 The corresponding covariance function only depends on the distance
 \eqn{r \ge 0}{r \ge 0} between 
 two points and is given by
 \deqn{C(r) = (1 + 8 s r + 25 s^2 r^2 + 32 s^3 r^3)(1-s r)^8 }{C(r) = (1 + 8 s r + 25 s^2 r^2 + 32 s^3 r^3)(1-s r)^8 }
 if \eqn{0 \ge r \ge \frac{1}{s}}{0 <= r <= 1/s} and 
 \deqn{C(r)=0}{C(r)=0} otherwise. Here,
  \eqn{s=0.301187465825}{s=0.301187465825}.
 For a generalized model see also \command{\link{RMgengneiting}}.
}
\usage{
RMgneiting(orig, var, scale, Aniso, proj)
}
\arguments{
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
   \command{\link{RMmodel}}. If not passed, the above
   covariance function remains unmodified.}
 \item{orig}{logical.
   if \code{TRUE} the above model is used.
   Otherwise
   the \code{\link{RMgengneiting}} model C(s r) with \code{kappa=3} as above, but
   but with \code{mu = 2.683509} and \code{s=0.2745640815} is used.
   The latter has the advantage of being closer to the Gaussian model and
   it is valid up to dimension 5.

   Default: \code{TRUE}   
 }
}
\details{
 This isotropic covariance function is valid only for dimensions less than or equal to 3. 
It is 6 times differentiable and has compact support.

 This model is an alternative to \command{\link{RMgauss}} as its graph is hardly distinguishable from the 
graph of the Gaussian model, but possesses neither the mathematical nor the numerical disadvantages of the Gaussian model.

It is a special case of \command{\link{RMgengneiting}} for the choice \eqn{kappa=3, \mu=1.5}{kappa=3, mu=1.5}.

Note that, in the original work by Gneiting (1999), a numerical value slightly deviating from
the optimal one was used for \eqn{mu=1.5}: \eqn{s=\frac{10
    \sqrt(2)}{47}}{s=10 sqrt(2)/47}.



}
\value{
 \command{\link{RMgneiting}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.
}
\references{
  For the original version
 \itemize{
 \item Gneiting, T. (1999)
 Correlation functions for atmospherical data analysis.
 \emph{Q. J. Roy. Meteor. Soc} Part A \bold{125}, 2449-2464.
 }

 For the version (\code{orig=FALSE})
 \itemize{
   \item this package \pkg{RandomFields}
 }
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMbigneiting}},
 \command{\link{RMgengneiting}},
 \command{\link{RMgauss}},
 \command{\link{RMmodel}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}




\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

plot(RMgneiting(), model2=RMgneiting(orig=FALSE), model3=RMgauss(), 
     xlim=c(-3,3), maxchar=100)
plot(RMgneiting(), model2=RMgneiting(orig=FALSE), model3=RMgauss(), 
     xlim=c(1.5,2.5), maxchar=100)

model <- RMgneiting(orig=FALSE, scale=0.4)
x <- seq(0,10, if (interactive()) 0.1 else 5)
z <- RFsimulate(model, x=x, y=x, z=x, T=c(1,1,4), maxGB=3)
plot(z, MARGIN.slices=4, MARGIN.movie=3)

\dontshow{FinalizeExample()}
}
