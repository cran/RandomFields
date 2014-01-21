\name{Advanced Max-stable random fields}
\alias{maxstableAdvanced}

\title{Simulation examples of advanced Max-Stable Random Fields}
\description{
  Here, some more advanced examples are given
}

%\details{
%}
\references{
  Strokorb, K. (2013) Ph.D. thesis.
}
\seealso{
  \command{\link{RPmaxstable}}
}
\examples{
RFoptions(seed=0)

\dontshow{\dontrun{

auswertung <- function(simu, threshold=2) { 
 below <- simu <= threshold
 freq <- mean(below[1, ])
 Print(freq, exp(-1/threshold)) ## univariate kontrolle
 both <- t(below) & below[1, ]
 plot(x, 2-log(colMeans(both)) / log(freq))

 ## alle 3 Linien ergeben das Gleiche:
 lines(x, RFcov(model=RRbrownresnick(model), x), col="yellow")
 lines(x, RFcov(model=RMextremalgauss(RMbr2eg(model)), x), col="red", lty=2)
 lines(x, RFcov(model=RMbinarygauss(RMbr2bg(model)), x), col="blue", lty=3)
}


model <- RMgauss(var= 4)
x <- seq(0, 5, 0.1)
z <- seq(0, 10, 0.1)


binary.model <- RPbernoulli(model)
bild <- RFsimulate(RPschlather(binary.model), z, z, grid=TRUE)
plot(bild)


simu <- RFsimulate(RPschlather(model), x, grid=TRUE, n=5000, # 1000000, 
                    spC=FALSE)
auswertung(simu)



bild <- RFsimulate(RPschlather(RMbr2eg(model)), z, z, grid=TRUE)
plot(bild)
simu <- RFsimulate(RPschlather(RMbr2eg(model)), x, grid=TRUE, 
                   spC=FALSE, n=5000) # n= 1000000,
auswertung(simu)





binary.model <- RPbernoulli(RMgaussprocess(RMbr2bg(model)))
bild <- RFsimulate(RPschlather(binary.model), z, z, grid=TRUE)
plot(bild)
simu <- RFsimulate(RPschlather(RMbr2bg(model)), x, grid=TRUE,
                   n=5000, # 1000000, 
                   spC=FALSE, maxstable.maxGauss=5)
auswertung(simu)


model <- RMgauss()
BR <- RFsimulate(RPbrownresnick(model), x, x, grid=TRUE)

}}
\dontshow{RFoptions(seed=NA)}

}
\author{
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}
