\name{SBS14}
\alias{SBS14}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{Systematic co-occurrence of tail correlation functions
among max-stable processes
}
\description{
 Here, the code of the paper on \sQuote{On some covariance models based on normal scale mixtures} is given.
}
 
\me

\references{
 \itemize{ 
 \item
 Strokorb, K., Ballani, F. and Schlather, M. (2014)
 Systematic co-occurrence of tail correlation functions
 among max-stable processes. Work in progress.
 }

}

\examples{\dontshow{StartExample()}
\dontshow{\dontrun{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

H <- c(0:9, LETTERS[1:6])
numbers <- rep(0:255, each=2)
col <- c(paste("#0000", H[1 + ((0:255) / 16)], H[1 + ((0:255) \%\% 16)],sep=""),
         paste("#", H[1 + (numbers / 16)], H[1 + (numbers \%\% 16)],
               H[1 + (numbers / 16)], H[1 + (numbers \%\% 16)], "FF", sep=""))

RFoptions(every=1000, xi=0)
x <- seq(0, 5, 0.01)
\dontshow{if(RFoptions()$internal$examples_red){warning("modified <x>"); x <- seq(0, 5, 1)}}
## Figure 2
model <- RMfbm(alpha=1, var=4) ## not 8 as in the paper, since semi-variogram
bm <- RMbrownresnick(model)


# Fig. 2, Brown-Resnick
z1 <- RFsimulate(RPbrownresnick(model), x, x)
plot(z1, col=col)
}}

\dontshow{\dontrun{
dev.copy2pdf(file="SBS_BRorig.pdf");  dev.print(device=jpeg,
file="SBS_BRorig.jpg", height=1000, width=1000, quality=100)
}}



\dontshow{\dontrun{
# Fig. 2, Monotone functions
z2 <- RFsimulate(RPsmith(RMm2r(bm), xi=0), x, x, 0)
plot(z2, col=col)
}}

\dontshow{\dontrun{
dev.copy2pdf(file="SBS_mono.pdf"); dev.print(device=jpeg,
file="SBS_mono.jpg", height=1000, width=1000, quality=100)
}}


\dontshow{\dontrun{
# Fig. 2, balls
z3 <- RFsimulate(RPsmith(RMm3b(bm)), x, x, 0, eps_zhou=0.1)
plot(z3, col=col)
}}
\dontshow{\dontrun{
  dev.copy2pdf(file="SBS_ball.pdf"); dev.print(device=jpeg,
file="SBS_ball.jpg", height=1000, width=1000, quality=100)
}}

\dontshow{\dontrun{
# Fig. 2,
z4 <- RFsimulate(RPsmith(RMmps(bm)), x, x, eps_zhou=0.1)
plot(z4, col=col)
}}
\dontshow{\dontrun{
dev.copy2pdf(file="SBS_poly.pdf"); dev.print(device=jpeg,
file="SBS_poly.jpg", height=1000, width=1000, quality=100)
}}


\dontshow{\dontrun{
## Figure 4
model <- RMexp(var=1.62) 
bm <- RMbrownresnick(model)
}}

\dontshow{\dontrun{
# Fig. 4, Brown-Resnick
z5 <- RFsimulate(RPbrownresnick(model), x, x)
plot(z5, col=col)
}}
\dontshow{\dontrun{
dev.copy2pdf(file="SBS_BRexp.pdf"); dev.print(device=jpeg,
file="SBS_BRexp.jpg", height=1000, width=1000, quality=100)
}}


\dontshow{\dontrun{
# Fig. 4, extremal Gaussian
z6 <- RFsimulate(RPschlather(bm), x, x)
plot(z6, col=col)
}}
\dontshow{\dontrun{
dev.copy2pdf(file="SBS_schlather.pdf");
dev.print(device=jpeg, file="SBS_schlather.jpg",
height=1000, width=1000, quality=100)
}}


\dontshow{\dontrun{
# Fig. 4, extremal binary 
z7 <- RFsimulate(RPschlather(RPbernoulli(bm)), x, x)
plot(z7, col=col)
}}
\dontshow{\dontrun{
dev.copy2pdf(file="SBS_binary.pdf");
dev.print(device=jpeg, file="SBS_binary.jpg",
height=1000, width=1000, quality=100)
}}



\dontshow{\dontrun{

# gamma_strokorb = 2 gamma_RandomFields
# -> ecf = erfc(sqrt(gamma_RandomFields / 4))

strokorb <- function(t) {
  ( 1 + 4 * t ) / ((2 * pi)^1.5 * (2 * t)^2.5) * exp(- 2 * t)
}


density <- function(x, u) {
  z <- apply(abs(x - u), 2, max)
#  Print(z, u, x)
  return(strokorb(min(z)))
}


RFoptions(xi=0)


model <- RMfbm(var=4, alpha=1) ## not 8 as in the paper as
## SEMI-variograms are used in RandomFields, not variograms
## as in several theoretical papers.

sm <- RMm2r(RMbrownresnick(model))

x <- seq(0, 10, 0.02)
 

z <- RFsimulate(RPsmith(RMm2r(RMbrownresnick(model)), xi=0),
                x, x, 0, n=1)
H <- c(0:9, LETTERS[1:6])
numbers <- rep(0:255, each=2)
col <- c(paste("#0000", H[1 + ((0:255) / 16)], H[1 + ((0:255) %% 16)], sep=""),
         paste("#", H[1 + (numbers / 16)], H[1 + (numbers %% 16)],
               H[1 + (numbers / 16)], H[1 + (numbers %% 16)], "FF", sep="")
         )
plot(z, type="p", pch=20, cex=0.4, col=col)

# dev.copy2pdf(file="strokorb.pdf") 

readline()
y <- as.double(as.matrix(z@data))
Print(range(y))
hist(y, 200, freq=FALSE)
a <- seq(-4, 10, 0.1)
lines(a, exp(-a) * exp(-exp(-a)))

############# test

system.time(
z <- RFsimulate(RPsmith(RMm2r(RMbrownresnick(model)), xi=0),
                seq(0, len=400, by=0.02), seq(0, len=400, by=0.02), 0, n=1,
                spConform=FALSE, every=0)
)

###

### Hintergrundbild
z <- RFsimulate(RPsmith(RMm2r(RMbrownresnick(model)), xi=0),
                seq(0, 1920/100, 0.02), seq(0, 1280/100, 0.02), 0, n=1,
                spConform=FALSE, every=0)
               
H <- c(0:9, LETTERS[1:6])
numbers <- rep(0:240, each=1)
col <- c(paste("#0000", H[1 + ((0:255) / 16)], H[1 + ((0:255) %% 16)], sep=""),
         paste("#", H[1 + (numbers / 16)], H[1 + (numbers %% 16)],
               H[1 + (numbers / 16)], H[1 + (numbers %% 16)], "FF", sep="")
         )
png(filename="~/Pictures/strokorb.png", width=1920, height=1280)#1930,1290
par(mar=rep(0,4))
image(z[,,1], col=col)


}}

\dontshow{FinalizeExample()}}


\keyword{spatial}

