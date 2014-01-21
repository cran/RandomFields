\name{RMaskey}
\alias{RMaskey}
\alias{RMtent}
\title{Askey model}
\description{
 Askey's model
 \deqn{C(x)= (1-x)^\alpha 1_{[0,1]}(x)}
 }
\usage{
RMaskey(alpha, var, scale, Aniso, proj)
RMtent(var, scale, Aniso, proj)
}
\arguments{
 \item{alpha}{a numerical value in the interval [0,1]}
 \item{var,scale,Aniso,proj}{optional parameters; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above covariance
 function remains unmodified.} 
}
\details{
 This covariance function is valid for dimension \eqn{d}{d} if
 \eqn{\alpha \ge (d+1)/2}. 
 For \eqn{\alpha=1} we get the well-known triangle (or tent)
 model, which is valid on the real line, only.
 }
\value{
 \command{\link{RMaskey}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}
}
\references{

 \itemize{
 \item Askey, R. (1973) \emph{Radial characteristic functions.}
 Technical report, Research Center, University of Wisconsin-Madison.

 \item
 Golubov, B. I. (1981) On Abel-Poisson type and Riesz means,
 \emph{Anal. Math.} 7, 161-184. 
 }


 Applications can be found in
 \itemize{
 \item Gneiting, T. (1999) Correlation functions for atmospheric data
 analysis. \emph{Quart. J. Roy. Meteor. Soc.}, 125:2449-2464.

 \item
 Gneiting, T. (2002) Compactly supported correlation functions.
 \emph{J. Multivar. Anal.}, 83:493-508.

 \item Strokorb, K., Ballani, F., and Schlather, M. (2014)
 In Preparation.
 
 \item
 Wendland, H. (1994) \emph{Ein Beitrag zur Interpolation mit radialen
 Basisfunktionen.} Diplomarbeit, Goettingen.

 \item
 Wendland, H.
 Piecewise polynomial, positive definite and compactly supported radial
 functions of minimal degree. Adv. Comput. Math., 4:389-396, 1995.

 }
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}},
 \command{\link{RMbigneiting}},
 \command{\link{RMgengneiting}},
 \command{\link{RMgneiting}},
 \command{\link{RFsimulate}},
 \command{\link{RFfit}}.
}


\keyword{spatial}
\keyword{models}
\examples{
RFoptions(seed=0)
model <- RMtent()
x <- seq(0, 10, if (interactive()) 0.02 else 1) 
plot(model, ylim=c(0,1))
plot(RFsimulate(model, x=x))
\dontshow{RFoptions(seed=NA)}
}