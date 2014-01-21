\name{SBS14}
\alias{SBS14}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{Systematic co-occurrence of tail correlation functions
among max-stable processes
}
\description{
 Here the code of the paper on \sQuote{On some covariance models based on normal scale mixtures} is given.
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de};
}
\references{
 \itemize{ 
 \item
 Strokorb, K., Ballani, F. and Schlather, M. (2014)
 Systematic co-occurrence of tail correlation functions
 among max-stable processes. Work in progress.
 }

}

\examples{
RFoptions(seed=0)
\dontshow{if (!.C("isAuthor", a=integer(1))$a) {}}

%  library(RandomFields, lib="~/TMP")
RFoptions(every=1000, xi=0)
model <- RMfbm(alpha=1, var=1) ## not 8 as in the paper
x <- seq(0, 5, if (interactive()) 0.01 else 1)

# Fig. 2, Brown-Resnick
z <- RFsimulate(RPbrownresnick(model), x, x)
plot(z)
\dontshow{if (.C("isAuthor", a=integer(1))$a && interactive()) dev.copy2pdf(file="SBS_BRorig.pdf")}


# Fig. 2, Monotone functions
bm <- RMbrownresnick(model)
z <- RFsimulate(RPsmith(RMstrokorbMono(bm), xi=0), x, x, 0)
H <- c(0:9, LETTERS[1:6])
numbers <- rep(0:255, each=2)
col <- c(paste("#0000", H[1 + ((0:255) / 16)], H[1 + ((0:255) \%\% 16)], sep=""),
         paste("#", H[1 + (numbers / 16)], H[1 + (numbers \%\% 16)],
               H[1 + (numbers / 16)], H[1 + (numbers \%\% 16)], "FF", sep=""))
plot(z, col=col)
\dontshow{if (.C("isAuthor", a=integer(1))$a && interactive()) dev.copy2pdf(file="SBS_mono.pdf")}


# Fig. 2, balls
\dontshow{if (.C("isAuthor", a=integer(1))$a && interactive()) dev.copy2pdf(file="SBS_ball.pdf")}


# Fig. 2,
\dontshow{if (.C("isAuthor", a=integer(1))$a && interactive()) dev.copy2pdf(file="SBS_poly.pdf")}





\dontshow{\dontrun{

strokorb <- function(t) {
  ( 1 + 1 * t ) / (4 * (2 * pi)^1.5 * t^2.5) * exp(- 0.5 * t)
}


density <- function(x, u) {
  z <- apply(abs(x - u), 2, max)
#  Print(z, u, x)
  return(strokorb(min(z)))
}


RFoptions(xi=0)
model <- RMfbm(var=4, alpha=0.65) ## not 8 as in the paper
model <- RMfbm(var=4,alpha=1) ## not 8 as in the paper
model <- RMfbm(alpha=1) ## not 8 as in the paper
sm <- RMstrokorbMono(RMbrownresnick(model))

x <- seq(0, 10, 0.02)


z <- RFsimulate(RPsmith(RMstrokorbMono(RMbrownresnick(model)), xi=0),
                x, x, 0, n=1, print=10)
H <- c(0:9, LETTERS[1:6])
numbers <- rep(0:255, each=2)
col <- c(paste("#0000", H[1 + ((0:255) / 16)], H[1 + ((0:255) %% 16)], sep=""),
         paste("#", H[1 + (numbers / 16)], H[1 + (numbers %% 16)],
               H[1 + (numbers / 16)], H[1 + (numbers %% 16)], "FF", sep="")
         )
plot(z, type="p", pch=20, cex=0.4, col=col)

# dev.copy2pdf(file="strokorb.pdf")

getOption("device")()
y <- as.double(as.matrix(z@data))
Print(range(y))
hist(y, 200, freq=FALSE)
a <- seq(-4, 10, 0.1)
lines(a, exp(-a) * exp(-exp(-a)))

if (!interactive()) .C("sleepMicro", as.integer(3e6))

############# test

unix.time(
z <- RFsimulate(RPsmith(RMstrokorbMono(RMbrownresnick(model)), xi=0),
                seq(0, len=400, by=0.02), seq(0, len=400, by=0.02), 0, n=1,
                spC=FALSE, every=0)
)

###

### Hintergrundbild
z <- RFsimulate(RPsmith(RMstrokorbMono(RMbrownresnick(model)), xi=0),
                seq(0, 1920/100, 0.02), seq(0, 1280/100, 0.02), 0, n=1,
                spC=FALSE, every=0)
               
H <- c(0:9, LETTERS[1:6])
numbers <- rep(0:240, each=1)
col <- c(paste("#0000", H[1 + ((0:255) / 16)], H[1 + ((0:255) %% 16)], sep=""),
         paste("#", H[1 + (numbers / 16)], H[1 + (numbers %% 16)],
               H[1 + (numbers / 16)], H[1 + (numbers %% 16)], "FF", sep="")
         )
png(filename="~/Pictures/strokorb.png", width=1920, height=1280)#1930,1290
par(mar=rep(0,4))
image(z[,,1], col=col)
dev.off()

}}


\dontshow{RFoptions(seed=NA)}

}


\keyword{spatial}

