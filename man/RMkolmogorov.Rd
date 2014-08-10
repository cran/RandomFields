\name{RMkolmogorov}
\alias{RMkolmogorov}
\title{Identical Model}
\description{
 \command{RMkolmogorov} corresponds to a vector-valued random fields with
 covariance function
 \deqn{
 \gamma_{ij}(h) = \|h\|^{2/3}\left(\frac43 \delta_{ij} -
 \frac13\frac{h_ih_j}{\|h\|^2}\right)}{
 \gamma_{ij}(h) = |h|^{2/3}(4/3 * \delta_{ij} -
 1/3 * h_i h_j / \|h\|^2)
 } 
}

\usage{
 RMkolmogorov(var, scale, Aniso, proj)
}

\arguments{
 \item{var,scale,Aniso,proj}{optional arguments; same meaning for any
 \command{\link{RMmodel}}. If not passed, the above
 covariance function remains unmodified.}
}


\value{
 \command{\link{RMkolmogorov}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}
\references{
 The above formula is from eq. (6.32) of section 6.2 in
 
 Pope, S.B. (2011) \emph{Turbulent Flows.}
 \bold{Cambridge:} Cambridge University Press.
}


\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\seealso{
 \command{\link{RMmodel}}, \command{\link{RMcurlfree}},
 \command{\link{RMdivfree}},
 \command{\link{RMvector}},
}

\keyword{spatial}
\keyword{models}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
x <- y <- seq(-2, 2, len=if (interactive()) 20 else 2)
model <- RMkolmogorov()
plot(model, dim=3, MARGIN=1:2, fixed.MARGIN=1)

simu <- RFsimulate(model, x, y, z=0)
plot(simu, select.variables=list(c(1,2)), col=c("red"))
\dontshow{FinalizeExample()}
}




