#  source("technical.report.positive.def.fcts.R")
# cutoff, intrinsic, TBM2, spectral, add MPP

# q()

#library(SoPhy, lib="~/TMP")
#library(RandomFields, lib="~/TMP")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

fast <- !interactive()
fast <- TRUE

  p3 <- function(model, param, col=grey(100:0/100),
                 x=c(0, 10, if (fast) 1 else 10/127),
                yf = 1, zf = 1){
    par(mfcol=c(3,5), cex=0.5)

    if (!fast) Print(model)
    
    methods <- c("circ", "intrinsic", "TBM3", "direct", "coins")
    for (m in 1:length(methods)) {
      meth <- methods[m]
      Print(meth)
      xx <- x
      if (meth=="direct") xx[3] <- xx[3] * 12
                                        # if (meth=="intrinsic CE") break;
                                        #if (meth!="TBM3") next
      yy <- xx;   yy[2] = yy[2] * yf
      zz <- xx;   zz[2] = zz[2] * zf
      xa <- seq(xx[1], xx[2], xx[3])
      ya <- seq(yy[1], yy[2], yy[3])
      za <- seq(zz[1], zz[2], zz[3])
      z <- try(GaussRF(xx, yy, zz, grid=TRUE, gridtriple=TRUE, method=meth,
                       model=model, param=param, reg=m))
 #     str(GetRegisterInfo(m, TRUE), vec=15)
      if (is.numeric(z)) {
        par(mar=c(4, 4.5, 4.5, 0.5), cex=0.5)
        image(xa, ya, z[,,1],
              main=meth, col=col, zlim=c(-2.5, 2.5), xlim=range(xx),
              ylim=range(x)) # x not y !
        par(mar=c(0.5, 4.5, 4.5, 0.5))
        image(xa, za, z[,1,],
              main=meth, col=col, zlim=c(-2.5, 2.5), xlim=range(xx),
              ylim=range(x)) # x not y !
        par(mar=c(0.5, 4.5, 4.5, 0.5))
        image(ya, za, z[1,,],
              main=meth, col=col, zlim=c(-2.5, 2.5), xlim=range(xx),
              ylim=range(x)) # x not y !
     } else
      for (i in 1:3) plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), main=meth) 
    }
  }


model <- list(list(model="whittle", var=1, kappa=2,
                   aniso=matrix(nr=3, c(0.05, 0, 0, 0, 0.2, 0, 0, 0, 0.2) * 5 )))
if (FALSE)
p3(model=model)


## Rest muss noch umgeschrieben werden 


model <- list(list(model="whittle", var=1, kappa=2,
                   aniso=matrix(nr=3, c(0.2, 0, 0, 0, 0.2, 0, 0, 0, 0.05) * 5 )))
if (FALSE)
p3(model=model)


# system("sleep 3")


  p <- function(model, param, col=grey(100:0/100),
                x=c(0, 10, if (fast) 0.5 else 0.05),
                yf = 1){

    Print(model)
    
    par(mfcol=c(3,3))
    methods <- GetMethodNames()[1:9]
    mittel <- v <- numeric(length(methods))
    names(v) <- methods
    for (m in 1:length(methods)) {
      meth <- methods[m]
      print(meth)
      xx <- x
      if (meth=="direct matrix decomposition") xx[3] <- xx[3] * 8
                                        # if (meth=="intrinsic CE") break;
                                        #if (meth!="TBM3") next
      yy <- xx;
      yy[2] = yy[2] * yf
      xa <- seq(xx[1], xx[2], xx[3])
      ya <- seq(yy[1], yy[2], yy[3])
      if (!fast) Print(xx, yy, model)
      z <- try(GaussRF(xx, yy, grid=TRUE, gridtriple=TRUE, method=meth,
                       model=model, param=param, reg=m))
#      str(GetRegisterInfo(m, TRUE), vec=15)
      if (is.numeric(z)) {
        par(mar=c(0.5, 4.5, 4.5, 0.5))
        image(xa, ya, z, main=meth, col=col, zlim=c(-2.5, 2.5), xlim=range(xx),
              ylim=range(x)) # x not y !
        v[m] <- var(as.vector(z))
        mittel[m] <- mean(z)
      } else plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), main=meth) 
    }
    print(cbind(v, mittel))
    invisible(cbind(v, mittel))
  }

m <- function(model, param, col=grey(100:0/100), meth, ...){
  x <- c(0, 10, if (fast) 0.5 else 0.05)
  xx <- x
  if (meth=="direct matrix decomposition") xx[3] <- xx[3] * 8
  y <- seq(xx[1], xx[2], xx[3])
  z <- try(GaussRF(xx, xx, grid=TRUE, gridtriple=TRUE, method=meth,
                    model=model, param=param,...))
  if (is.numeric(z)) {
    image(y, y, z, main=meth, col=col, zlim=c(-2.5, 2.5))
  } else plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), main=meth)
  print(range(z))
}

if (interactive() || file.exists("../../../SOPHY/makefile"))
  do.call(getOption("device"), list(height=4.3, width=4.3))

RFparameters(CE.force=TRUE, TBMCE.force=TRUE, CE.trials=1, TBMCE.trials=1,
             Print=4 * 0, Storing=TRUE, spectral.lines=500)


gneitingdiff <-  function(p, op="*"){
   list(list(m="gneiting", v=p[2], s=p[6]*p[4]*10*sqrt(2)/47),
       op,
       list(m="whittle", k=p[5], v=1.0, s=p[4])
        )
}

if (FALSE)
  {
  x <- 0
  for (i in 1: 100) {
    print(i)
    x <- x + p(model="gauss", param=c(0, 1, 0, 1/0.875))
    print(c(i, as.double(x) / i), dig=2)
  }
# 1.78 und 1.45 als varianz fuertbm2 und 3
}

# str(RFparameters())

#model = list(list(model="stable", var=2, kappa=1.5, scale=1),
model = list(list(model="stable", var=2, kappa=1.5, aniso=c(1,0,0,1)))
#  "+"
#  list(model="stable", var=0, kappa=0.2, aniso=c(0,1,-1, 0))

##str(GetRegisterInfo(5), vec=20)
#system("sleep 3")

# DeleteAllRegisters()

## 3d bilder -- anisotropien richtig?


#system("sleep 3")
#q()



RFparameters(PracticalRange=11)
model <- list(list(model="whittle", var=1, kappa=2,
                   aniso=matrix(nr=2, c(0.05, 0, 0, 0.2))))
if (FALSE)
p(model=model, x=c(-10, 110, 1), yf=0.5)
RFparameters(PracticalRange=FALSE)


model <- list(list(model="whittle", var=0.25, kappa=2,
                   aniso=matrix(nr=2, c(0.05, 0, 0, 0.2))))
if (FALSE)
p(model=model, x=c(-10, 110, 1))


model = list(list(model="stable", var=2,kappa=2,aniso=matrix(nr=2,c(1,4,0,1.5))) ,
  "+",
  list(model="stable", var=1, kappa=1.5, aniso= matrix(nr=2, c(1,4,0,1.5) * 4)),
  "+",
  list(model="nugget", var=1, aniso=matrix(nr=2, c(1,1,1,1)))
  )
p(model=model)

model = list(list(model="stable", var=2, kappa=2, aniso=matrix(nr=2, c(1,4,0,1.5))) ,
  "+",
  list(model="stable", var=1, kappa=1.5, aniso=matrix(nr=2, c(1,4,0,1.5) * 4)),
  "+",
  list(model="nugget", var=10, aniso=matrix(nr=2, c(1,1,1,1)))
  )
p(model=model)

model = list(list(model="stable", var=2, kappa=2,
  aniso= matrix(nr=2, c(1,4,0,1.5))) ,
  "+",
  list(model="stable", var=1, kappa=1.5, aniso= matrix(nr=2, c(1,4,0,1.5) * 0.4))
  )
p(model=model)


model = list(list(model="stable", var=2, kappa=2, aniso=matrix(nr=2, c(1,4,0,1.5))) ,
  "+",
  list(model="stable", var=1, kappa=1.5, aniso=matrix(nr=2, c(1,4,0,1.5) * 4))
  )
p(model=model)

model = list(list(model="stable", var=2, kappa=1.5, aniso=matrix(nr=2, c(1,4,0,1.5))))
p(model=model)

model = list(list(model="stable", var=2, kappa=1.5, aniso=matrix(nr=2, c(0.5,0,0,4))))
p(model=model)

m(model="exp", param=c(0, 1, 0, 2/3), meth="hyper", hyper.superpos=1)
m(model="exp", param=c(0, 1, 0, 2/3), meth="hyper", hyper.superpos=2)
m(model="exp", param=c(0, 1, 0, 2/3), meth="hyper", hyper.superpos=10)
m(model="exp", param=c(0, 1, 0, 2/3), meth="hyper", hyper.superpos=100)


m(model="circ", param=c(0,1,0, 5/2), meth="coins", mpp.intens=1)
m(model="circ", param=c(0,1,0, 5/2), meth="coins", mpp.intens=10)
m(model="circ", param=c(0,1,0, 5/2), meth="coins", mpp.intens=100)
m(model="circ", param=c(0,1,0, 5/2), meth="coins", mpp.intens=1000)

p(model="nugget", param=c(0, 1, 0, 1/0.4))

p(model="whittle", param=c(0, 1, 0, 1/2.0, 1))
p(model="whittle", param=c(0, 1, 0, 1/2.75, 2))
p(model="whittle", param=c(0, 1, 0, 1/3.65, 4))

p(model="stable", param=c(0, 1, 0, 1/4.5, 0.5))
p(model="stable", param=c(0, 1, 0, 1/1.0, 1.5))
p(model="stable", param=c(0, 1, 0, 1/0.9, 1.9))

p(model="sph", param=c(0, 1, 0, 1/0.4))

p(model="gauss", param=c(0, 1, 0, 1/0.875))

p(model="exp", param=c(0, 1, 0, 2/3))

p(model="gneiting", param=c(0, 1, 0, 4.7 * 10*sqrt(2)/47))

p(model=gneitingdiff(c(0, 1, 0, 0.833, 4, 5)))
p(model=gneitingdiff(c(0, 1, 0, 0.889, 2, 5)))
p(model=gneitingdiff(c(0, 1, 0, 1, 1, 5)))
p(model=gneitingdiff(c(0, 1, 0, 1, 1, 5), "+"), col=rainbow(100))
p(model=gneitingdiff(c(0, 1, 0, 1.38, 2, 3)))
p(model=gneitingdiff(c(0, 1, 0, 2, 2, 2)))


p(model="cauchy", param=c(0, 1, 0, 1/10, 0.5))
p(model="cauchy", param=c(0, 1, 0, 1/1.25, 1.5))
p(model="cauchy", param=c(0, 1, 0, 1/0.575, 3.5))

p(model="bessel", param=c(0, 1, 0, 1/3.5, 4))
p(model="bessel", param=c(0, 1, 0, 1/3.5, 2))
p(model="bessel", param=c(0, 1, 0, 1/3.5, 1))
p(model="bessel", param=c(0, 1, 0, 1/3.5, 0.5))
p(model="bessel", param=c(0, 1, 0, 1/3.5, 0))


## anisotrope modelle
## additive modelle und multiplicative modelle

