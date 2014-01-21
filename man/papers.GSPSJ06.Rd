\name{GSPSJ06}
\alias{GSPSJ06}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{Fast and Exact Simulation of Large Gaussian Lattice Systems
 in R2}
\description{
 Here the code of the paper on \sQuote{Fast and Exact Simulation of Large Gaussian Lattice Systems
 in R2} is given.
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de};
}
\references{
 \itemize{ 
 \item
 Gneiting, T., Sevcikova, H., Percival, D.B., Schlather, M., Jiang,
 Y. (2006) Fast and Exact Simulation of Large Gaussian Lattice Systems
 in R2: Exploring the Limits. \emph{J. Comput. Graph. Stat.},
 \bold{15}, 483-501.

 }
 
 % Schlather, M. (2001) Simulation of stationary and isotropic random
% fields. \emph{R-News} \bold{1} (2), 18-20.


}

\examples{
RFoptions(seed=0)
\dontshow{if (.C("isAuthor", a=integer(1))$a) {

## Figure 1 

stabletest <- function(alpha, theta, size=if (interactive()) 512 else 20) {
  RFoptions(trials=1, tolIm = 1e-8, tolRe=0, force = FALSE,
            useprimes=TRUE, strategy=0, skipchecks=!FALSE,
            printlevel=2, storing=TRUE)
  model <- RMcutoff(diameter=theta, a=1, RMstable(alpha=alpha))
  RFcov(dist=0, model=model, dim=2)
  r <- RFgetModelInfo(RFcov, modelname="RMcutoff")$internalq[5] # theor R
  x <- c(0, r, r / (size - 1)) * theta
  err <- try(RFsimulate(x, x, grid=TRUE, model=RPcirculant(model), n=0))
  return(if (!is.numeric(err) || err != 0) NA else r)
}

alphas <- seq(1.52, 2.0, if (interactive()) 0.02 else 0.1) 
thetas <-
  if (interactive()) seq(0.05, 3.5, 0.05) else seq(0.1, 3.5, 0.2) 

m <- matrix(NA, nrow=length(thetas), ncol=length(alphas))
for (it in 1:length(thetas)) {
  theta <- thetas[it]
  for (ia in 1:length(alphas)) {
  alpha <- alphas[ia]
  cat("alpha=", alpha, "theta=", theta,"\n")
  print(unix.time(m[it, ia] <- stabletest(alpha=alpha, theta=theta)))
  if (is.na(m[it, ia])) break
  }
  if (any(is.finite(m))) image(thetas, alphas, m, col=rainbow(100))
}

}}
\dontshow{RFoptions(seed=NA)}
}

\keyword{spatial}

