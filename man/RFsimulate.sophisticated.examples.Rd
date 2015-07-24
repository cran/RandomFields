
\name{RFsimulate.sophisticated.examples}
\alias{RFsimulate.sophisticated.examples}
\title{Sophisticated Examples for the Simulation of Random Fields}

\description{
 This man page will give a collection of basic examples for the use of
 \code{\link{RFsimulate}}.

 For other kind of random fields (binary, max-stable, etc.) or
 more sophisticated approches see \link{RFsimulateAdvanced}.
}


\seealso{
 \command{\link{RFsimulate}},
 \command{\link{RFsimulateAdvanced}}
}


\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample()}

\dontshow{\dontrun{


#############################################################
##                                                         ##
## Example 1: Gaussian field approximated by Poisson fields##
##                                                         ##
#############################################################

for (mpp.intensity in c(1, 10, 100)){
  # mpp.intensity of 1 is much too small but illustrates
  # how the "Coins" method works

  z <- RFsimulate(x=x, model=RPcoins(RMspheric()),
  mpp.intensity = mpp.intensity)
  #getOption("device")()
  plot(z)
  readline("press return")
}

par(mfcol=c(2,1))
plot(z@data[,1:min(length(z@data), 1000)], type="l")
hist(z@data[,1])



#############################################################
##                                                         ##
## Example 2: A max-stable random field                    ##
##                                                         ##
#############################################################

### Smith's Gaussian extremal process
x <- GridTopology(0, .1, 500)
z <- RFsimulate(RPsmith(RMgauss()), x=x, n=10)
plot(z, nmax=3)

z <- as.vector(as.matrix(z@data))

par(mfcol=c(2,1))
plot(pmin(15, z[1:min(length(z), 1000)]), type="l")
hist(ylim=c(0,1), pmin(z,5), 200, freq=FALSE)
xx <- seq(0,4,len=1000)
lines(xx, exp(-1/xx) / xx^2)
 
 
## a more complicated mixed moving maximum process
model <- RPsmith(RMmppplus(RMgauss(), RMexp(), p=c(0.3, 0.7)))
z <- RFsimulate(model, x=x, n=10)
 
plot(z, nmax=1, ylim=c(0, 15))

z<-as.vector(as.matrix(z@data))
par(mfcol=c(2,1))
plot(pmin(15, z[1:min(length(z), 1000)]), type="l")
hist(ylim=c(0,1), pmin(z,5), 200, freq=FALSE)
xx <- seq(0,4,len=1000)
lines(xx, exp(-1/xx) / xx^2)

 
## there are different possibilities to define a max-stable process:
## * m[[1]] below is a detailed way of defining a model.
## * m[[2]] is the same as m[[1]] as only one component is given
## * m[[3]] uses the fact that the standard schlather model is based 
##          on a Gaussian random field. So, it suffices to pass the
##          covariance model

m <- list(RMmppplus(RPgauss(RMgauss())),
          RPgauss(RMgauss()),
          RMgauss())

 x <- GridTopology(0, .1, 500)

for (i in 1:3){ %# same seed always
  z <- RFsimulate(model=Schlather(m[[i]]),x=x, n=2, seed=0)

  plot(z, nmax=1, ylim=c(0, 15))
  z <- as.vector(as.matrix(z@data))

  par(mfcol=c(2,1))
  plot(pmin(10, z[1:min(length(z), 1000)]), type="l")
  hist(ylim=c(0,1), pmin(z,5), 200, freq=FALSE)
  xx <- seq(0,4,len=1000)
  lines(xx, exp(-1/xx) / xx^2)
  print(quantile(as.vector(z), probs=seq(0,1,0.05)))
}
 
 
## mixture of extremal Gaussian models:
x <- GridTopology(0, .03, 500)
model <- RMmppplus(RPgauss(RMgauss()), RPgauss(RMexp()),
                   p = c(0.7, 0.3))

z <- RFsimulate(model = Schlather(model), x=x, 
gauss.meth="sp", n=1)
plot(z)

z <- as.vector(as.matrix(z@data))
par(mfcol=c(2,1))
plot(pmin(1000, z[1:min(length(z), 1001)]), type="l")
hist(ylim=c(0,1), pmin(z, 5), 200, freq=FALSE)
xx <- seq(0,4,len=1000)
lines(xx, exp(-1/xx) / xx^2)
print(summary(z))



## non-separable space-time model applied for two space dimensions
## note that tbm method works in some special cases.
#% library(RandomFields, lib="~/TMP")
x <- y <- seq(0, 7, 0.05)
T <- c(1,32,1) * 10 ## note necessarily gridtriple definition
model <- RMnsst(aniso=diag(c(3, 3, 0.02)), delta=2,
                phi=RMgauss(), psi=RMgenfbm(alpha=1, delta=0.5))
z <- RFsimulate(x=x, y=y, T=T, model=model,
                method="ci", CE.strategy=1, CE.trials=4)
rl <- function() readline("Press return")
for (i in 1:dim(z)[3]) { image(z[,,i], zlim=range(z)); rl();}
for (i in 1:dim(z)[2]) { image(z[,i,], zlim=range(z)); rl();}
for (i in 1:dim(z)[1]) { image(z[i,,], zlim=range(z)); rl();}



#############################################################
## Example 3 shows the benefits from stored,               ##
## intermediate results: in case of the circulant          ##
## embedding method, the speed is doubled in the second    ##
## simulation.                                             ## 
#############################################################

RFoptions(storing=TRUE)
y <- x <- seq(0, 50, 0.1)
(p <- c(runif(3), runif(1)+1))
ut <- system.time(f <- RFsimulate(RPcirculant(RMexp())), x=x, y=y)
% method="circ", param=p))
plot(f) 
%hist(f)
%c( mean(as.vector(f)), var(as.vector(f)) )
cat("system time (first call)", format(ut,dig=3),"\n")

# second call with the same paramters can be much faster:
ut <- system.time(f <- RFsimulate()) 
plot(f) 

%hist(f)
%c( mean(as.vector(f)), var(as.vector(f)) )
cat("system time (second call)", format(ut,dig=3),"\n")

#############################################################
##                                                         ##
## Example 4: how the cutoff method can be set             ##
## explicitly using hypermodels                            ##
##                                                         ##
#############################################################

## NOTE: this feature is still in an experimental stage
## which has not been yet tested intensively;
## further: arguments and algorithms may change in
## future.

#% library(RandomFields, lib="~/TMP");source("~/R/RF/RandomFields/tests/source.R")
## simuation of the stable model using the cutoff method
x <- seq(0, 1, 1/24) # nicer pictures with 1/240
scale <- 1.0
model1 <- RPcutoff(RMstable(alpha=1, scale=scale))
rs <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
z1 <- RFsimulate(x, x, model=model1, n=1, storing=TRUE)
(size <- RFgetRegisterInfo(meth=c("cutoff", "circ"))$S$size)
cut.off.param <- RFgetRegisterInfo(meth=c("cutoff", "circ"),
                                   modelname="cutoff")$param
print(cut.off.param)
plot(z1)

## simulation of the same random field using the circulant
## embedding method and defining the hypermodel explicitely
model2 <- RMcutoff(scale = scale, diam=cut.off.param$diam, a=cut.off.param$a, 
                   RMstable(alpha=1.0))
		 
assign(".Random.seed", rs, envir=.GlobalEnv)
z2 <- RFsimulate(x, x, gridtriple=FALSE, model=model2,
                 meth="circulant", n=1, CE.mmin=size, Storing=TRUE)
image(x, x, z2)
Print(range(z1-z2)) ## essentially no difference between the fields!



#############################################################
## Example 5:                                              ##
## The cutoff method simulates on a torus and a (small)    ##
## rectangle is taken as the required simulation.          ##
##                                                         ##
## The following code shows a whole such torus.            ##
## The main part of the code sets local.dependent=TRUE and ##
## local.mmin to multiples of the basic rectangle lengths  ##
#############################################################

# definition of the realisation
RFoptions(circulant.useprimes=FALSE)
x <- seq(0, 2, len=4) # better 20
y <- seq(0, 1, len=5) # better 40
grid.size <- c(length(x), length(y))
model <- RMexp(var=1.1, aniso=matrix(nc=2, c(2, 1, 0.5, 1)))

# determination of the (minimal) size of the torus
InitRFsimulate(x, y, model=model, method="cutoff")
ce.info.size <- RFgetRegisterInfo(meth=c("cutoff", "circ"))$S$size
blocks <- ceiling(ce.info.size / grid.size / 4) *4
(size <- blocks * grid.size)

# simulation and plot of the torus 
z <- RFsimulate(x, y, model=model, method="cu",
                n=prod(blocks) * 2,
                local.dependent=TRUE, local.mmin=size)[,,c(TRUE, FALSE)]
height <- 8
ScreenDevice(height=height,
             width=height / blocks[2] / diff(range(y)) *
                   blocks[1] * diff(range(x))))

close.screen(all = TRUE)
sc <- matrix(nc=blocks[1], split.screen(rev(blocks)), byrow=TRUE)
sc <- as.vector(t(sc[nrow(sc):1, ]))

for (i in 1:prod(blocks)) {
  screen(sc[i])
  par(mar=rep(1, 4) * 0.0)
  image(z[,, i], zlim=c(-3, 3), axes=FALSE, col=rainbow(100)) 
}

% folgender Befehl muss unbedingt drin bleiben
close.screen(all = TRUE)

}}

\dontshow{FinalizeExample()}

}