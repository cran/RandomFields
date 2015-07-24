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
 \url{http://ms.math.uni-mannheim.de/de/publications/software};
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
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}

## Figure 1 (pretty time consuming)
stabletest <- function(alpha, theta, size=512) {
  RFoptions(trials=1, tolIm = 1e-8, tolRe=0, force = FALSE,
            useprimes=TRUE, strategy=0, skipchecks=!FALSE,
             storing=TRUE)
  model <- RMcutoff(diameter=theta, a=1, RMstable(alpha=alpha))
  RFcov(dist=0, model=model, dim=2, seed=0)
  r <- RFgetModelInfo(modelname="RMcutoff", level=3)$internalq[5] # theor R
  x <- seq(0, r, by= r / (size - 1)) * theta
  err <- try(RFsimulate(x, x, model=RPcirculant(model), n=0))
  return(if (class(err) == "try-error") NA else r)
}

alphas <- seq(1.52, 2.0, 0.02) 
thetas <- seq(0.05, 3.5, 0.05)
\dontshow{if (RFoptions()$internal$examples_reduced) {
   warning("reduced size of alphas and thetas")
   alphas <- seq(1.52, 2.0, 0.5) 
   thetas <- seq(0.1, 3.5, 2)
}
}
m <- matrix(NA, nrow=length(thetas), ncol=length(alphas))
for (it in 1:length(thetas)) {
  theta <- thetas[it]
  for (ia in 1:length(alphas)) {
  alpha <- alphas[ia]
  cat("alpha=", alpha, "theta=", theta,"\n")
  m[it, ia] <- stabletest(alpha=alpha, theta=theta)
  if (is.na(m[it, ia])) break
  }
  if (any(is.finite(m))) image(thetas, alphas, m, col=rainbow(100))
}

\dontshow{FinalizeExample()}
}

\keyword{spatial}

