\name{Max-stable random fields, advanced }
\alias{maxstableAdvanced}
\alias{RPmaxstableAdvanced}


\title{Simulation examples of advanced Max-Stable Random Fields}
\description{
  Here, an advanced example is given used to test whether the
  algorithms work correctly.
}

%\details{
%}
\references{
  Strokorb, K. (2013) Ph.D. thesis.
}
\seealso{
  \command{\link{RPmaxstable}}
}

\me
\keyword{spatial}


\examples{\dontshow{StartExample()}%\dontshow{StartExample(reduced=50)}
\dontshow{\dontrun{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

n <- 100

model <- RMexp(var=1.62 / 2) 
x <- seq(0, 5, 0.2)
y <- seq(0, 10, 0.2)


auswertung <- function(simu, model, threshold=2) {
  x <- as.vector(coordinates(simu))
  simu <- as.array(simu)
  below <- simu <= threshold
  freq <- rowMeans(below)
  meanfreq <- mean(freq)
  Print(freq, meanfreq, exp(-1/threshold)) ## univariate kontrolle
  both <- t(below) & below[1, ]
  ecf <-  2-log(colMeans(both)) / log(meanfreq)
  plot(x, ecf)

  ## alle 3 Linien ergeben das Gleiche:
  spC <- RFoptions()$general$spConform
  RFoptions(spConform = FALSE)
  lines(x, m1 <- RFcov(RMbrownresnick(model), x), col="yellow")
  lines(x, m2 <- RFcov(RMschlather(RMbr2eg(model)), x), col="red", lty=2) # OK
  m3 <- RFcov(RMbernoulli(RMbr2bg(model), centred=FALSE), x)
  lines(x, m3, col="blue", lty=3)
  RFoptions(spConform = spC)

  erfc <- function(x) 2 * pnorm(x, 0, sd=1/sqrt(2), lower=FALSE)
  lines(x, m4 <- erfc(0.45 * sqrt(1-exp(-x))), lty=4)
 
  ## theoretical curves correct?
  if (!all.equal(m1, m2) || !all.equal(m1, m3) || !all.equal(m1, m4))
    stop("calculation error")

  if ( (n <- ncol(simu)) >= 1000) {
    ## margins correct?
    mar.threshold <- 4 * 0.5 / sqrt(n)
    mmar.threshold <- 3 * 0.5 / sqrt(n)
    Print(abs(freq - exp(-1/threshold)), mar.threshold)
    if (abs(freq[sample(length(freq), 1)] - exp(-1/threshold)) > mar.threshold)
      stop("marginal distribution wrong? (single margin)")
    if (abs(meanfreq - exp(-1/threshold)) > mmar.threshold)
      stop("marginal distribution wrong? (mean margin)")
 
    ## extremal correlation function correct?
    meanthreshold <- 4 / sqrt(n)
    maxthreshold <- 2 * sqrt(nrow(simu)) / sqrt(n)
    Print(abs(ecf - m1), meanthreshold, maxthreshold)
    if (mean(abs(ecf - m2)) > meanthreshold)
      stop("ecf not correct? (mean deviation too large)")
    if (max(abs(ecf - m2)) > maxthreshold)
      stop("ecf not correct? (max deviation too large)")
  }
}


## BR currently not programmed
\dontshow{\dontrun{ ## to to
## Brown-Resnick
z <- RFsimulate(RPbrownresnick(model), y, y)
plot(z)
simu <- RFsimulate(RPbrownresnick(model), x, n=n,  max_gauss=5)
auswertung(simu, model)
}}

%\dontshow{\dontrun{ ## Windows hat Fehler. Unklar wieso
## Extremal Gaussian
z <- RFsimulate(RPschlather(RMbr2eg(model)), y, y)
plot(z)
simu <- RFsimulate(RPschlather(RMbr2eg(model)), x,  n=n)
auswertung(simu, model)

## Extremal Binary Gaussian
binary.model <- RPbernoulli(RMbr2bg(model))
z <- RFsimulate(RPschlather(binary.model), y, y)
plot(z)
simu <- RFsimulate(RPschlather(binary.model), x, n=n, max_gauss=5)
auswertung(simu, model)


\dontshow{\dontrun{ ## OK!
zaehler <- 0
repeat {
  Print(zaehler)
  zaehler <- zaehler + 1
  simu <- RFsimulate(RPschlather(RMbr2eg(model)), x, spConform=FALSE, n=n,
                    pch="")
  auswertung(simu, model)
}
}}

}}

\dontshow{FinalizeExample()}

}