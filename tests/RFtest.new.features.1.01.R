# source("RFtest.new.features.1.01.R")


if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

cat("\n landmark -2")


RFparameters(Print=6)
#RFparameters(Print=1)
RFparameters(Print=0)

#

# .C("niceFFTnumber",as.integer(18 * 365))  ## gives 6615

IL <- function(x,nu) besselI(x,nu)-.C("StruveL",as.double(x),as.double(nu),as.integer(0))[[1]]
ILE <- function(x,nu) besselI(x,nu)- exp(x) * .C("StruveL",as.double(x),as.double(nu),as.integer(1))[[1]]
H <- function(x,nu) .C("StruveH",as.double(x),as.double(nu))[[1]]
E <- function(x,nu) besselI(x,nu-0.5)-.C("StruveL",as.double(x),as.double(0.5-nu),as.integer(0))[[1]]
E(10,2.5)

IL(0.1,0)
ILE(0.1,0)
IL(5,0)
ILE(5,0)
IL(3.2,1)
ILE(3.2,1) 


#stop("")
# TIME
DeleteAllRegisters()

RFparameters(CE.trial=3)


cat("\n landmark 0\n")

m <- matrix(c(0.00134,0,0, 0,0.00134,0, 0,0,0.9347),ncol=3) 
m <- matrix(c(0.00134,0, 0,0.9347),ncol=2) 
k <- c(1,phi=1, 1.544,0.61,psi=1, 2/0.61)
model <- list(list(m="nsst",v=1,k=k,a=m))
#m <- matrix(c(1,0.3, 0.2,0.1),ncol=2) / ff; k <- NULL
#model <- list(list(m="gneiting",v=1,k=k,a=m))
RFparameters(CE.usep=TRUE)
x <- c(1,450,5)
T <- c(1, .C("niceFFTnumber",as.integer(18 * 365))[[1]], 1)
unix.time(z <- GaussRF(x=x, T=T,
             grid=TRUE, model=model,
             me="ci",gridtriple=TRUE))
image(x=seq(x[1],x[2],x[3]),y=seq(T[1],T[2],T[3]),z,zlim=c(-3,3))

#stop("")
cat("\n landmark 1\n")

y <-  x <-  1:256
ff <- 10
m <- matrix(c(1,0, 0,1),ncol=2) / ff
k <- c(1,phi=1,1,0.5,psi=1,dim=3)
model <- list(list(m="nsst",v=1,k=k,a=m))
#m <- matrix(c(1,0.3, 0.2,0.1),ncol=2) / ff; k <- NULL
#model <- list(list(m="gneiting",v=1,k=k,a=m))

z <- GaussRF(x=x, T=c(1,256,1),grid=TRUE, model=model,me="ci")
image(z,zlim=c(-3,3))

#stop("")
cat("\n landmark 2\n")

xx <- matrix(c(expand.grid(x,y),recursive=TRUE),ncol=2)
cd <- CovarianceFct(xx,model)
image(matrix(cd,nrow=length(x)),zlim=c(0,1),col=grey((30:0)/30))

cat("\n landmark 3\n")

y <- x <- seq(-10,10,0.1)
xx <- matrix(c(expand.grid(x,y),recursive=TRUE),ncol=2)
model <- list(list(m="nsst",v=1,k=k,a=m))
cc <- CovarianceFct(xx,model)
image(matrix(cc,nrow=length(x)),zlim=c(-0.5,1))
dev.set()

y <- x <- seq(-10,10,0.1)
xx <- matrix(c(expand.grid(x,y),recursive=TRUE),ncol=2)
model <- list(list(m="nsst",v=1,k=k,a=m))
cd <- CovarianceFct(xx,model)
image(matrix(cc,ncol=sqrt(length(cd))),zlim=c(-0.5,1))

image(matrix((abs(cd) -abs(cc))/cd, ncol=sqrt(length(cd))),zlim=c(-1,1),
      col=rainbow(100))

Time <- c(min(y), max(y), diff(y)[1])
z <- GaussRF(x=x,T=Time,grid=TRUE, model=model)
image(z,zlim=c(-3,3))


#stop("")
cat("\n landmark 4\n")


RFparameters(TBMCE.trials=3, TBMCE.maxmem=2e7); # TBMCE.trials=5
x <- (1:32)/2
y <- (1:32)/2
T <- c(1,32,1)/2
ff <- 1
m <- matrix(c(1,0,0, 0,1,0, 0,0,1),ncol=3) / ff
k <- c(1,phi=1,1,0.5,psi=1,dim=3)
model <- list(list(m="nsst",v=1,k=k,a=m))

z <- GaussRF(x=x,y=y,T=T,grid=TRUE, model=model, me="TBM3", Print=5)


print(z[,,1])
print(z[,,1]-z[,,2])
if (FALSE) {
  for (i in 1:dim(z)[3]) { image(z[,,i]); readline();}
  for (i in 1:dim(z)[1]) { image(z[i,,]); readline();}
  for (i in 1:dim(z)[2]) { image(z[,i,]); readline();}
}  

if (!file.exists("/home/schlather/TMP/RandomFields")) q()


cat("\n landmark 5\n")

#stop("")


x <- (1:32)/2
y <- (1:32)/2
T <- c(1,32,1)/16
ff <- 1
m <- matrix(c(1,0,0, 0,1,0, 0,0,1),ncol=3) / ff
k <- c(1,phi=1,1,0.5,psi=1,dim=2)
model <- list(list(m="nsst",v=1,k=k,a=m))


print(c1 <- CovarianceFct(cbind(x,0,0),model))
print(c2 <- CovarianceFct(x,"stable",param=c(0,1,0,ff,1)))
c1 - c2

print(c1 <- CovarianceFct(cbind(0,0,seq(T[1],T[2],T[3])),model))
print(c2 <-CovarianceFct(seq(T[1],T[2],T[3]),"gencauchy",param=c(0,1,0,ff,1,1)))
c1 - c2

z1 <- GaussRF(x=x,y=y,T=T,grid=TRUE,model=model, me="ci")
print(z[,,1])
print(z[,,1]-z[,,2])
if (FALSE) {
  for (i in 1:dim(z1)[3]) { image(z1[,,i]); readline();}
  for (i in 1:dim(z1)[1]) { image(z1[i,,]); readline();}
  for (i in 1:dim(z1)[2]) { image(z1[,i,]); readline();}
}
print(CovarianceFct(cbind(0:10,0,0),model))

                                        #stop("")
cat("\n landmark 6\n")


x <- (1:32)/2
y <- (1:32)/2
T <- c(1,32,1)/16
ff <- 1
m <- matrix(c(1,0,0, 0,1,0, 0,0,1),ncol=3) / ff
k <- c(1,phi=1,1,1,psi=1,dim=2)
k <- NULL
model <- list(list(m="gneiting",v=1,k=k,a=m))

z <- GaussRF(x=x,y=y,T=T,grid=TRUE,
             model=model, me="ci")
print(z[,,1])
print(z[,,1]-z[,,2])
if (FALSE)
  for (i in 1:dim(z)[3]) { image(z[,,i]); readline();}

print(CovarianceFct(cbind(0:10,0,0),model))

#stop("")

cat("\n landmark 7\n")

## extremes
x <- (0:100)/10
x <- (0:30)/10

n <- 100
mx <- my <- 1:n

ms <- MaxStableRF(mx, my, grid=TRUE,
                  model=list(list(m="gauss", v=1,
                    a=matrix(c(1,2,3,4),ncol=2)/40)),
                  maxstable="B")
image(mx,my,sqrt(ms))



cat("\n landmark 8\n")

ms <- MaxStableRF(mx, my, grid=TRUE,
                  model="gauss", param=c(0, v=1, 0, s=40),
                  maxstable="B")
image(mx,my,sqrt(ms))


ms <- MaxStableRF(mx, my, grid=TRUE,
                  model=list(list(m="gauss", v=1, s=40)),
                    maxstable="B")
image(mx,my,sqrt(ms))


x <- (0:100)/10
x <- (0:30)/10

n <- 100
mx <- my <- 1:n

ms <- MaxStableRF(mx, my, grid=TRUE, model="exponen",
                  param=c(0,1,0,40), maxstable="extr")
image(mx,my,sqrt(ms))


ms <- MaxStableRF(mx, my, grid=TRUE,
                 model=list(list(m="expo",v=1,a=matrix(c(1,0,0,0.2)/10,ncol=2)),
                    "+",
                    list(m="gauss", v=1,a=matrix(c(0.2,0,0,1),ncol=2)/15)
                    ),
                    me="ci",maxstable="extr")
image(mx,my,sqrt(ms))

cat("\n landmark 9\n")



#stop("");

x <- (0:100)/10
x <- (0:30)/10

n <- 100



## additive model

z <- GaussRF(x=x,y=x,grid=TRUE,model="ci",param=c(0,1,0,1));image(z,zlim=c(-3,3))
z <- GaussRF(x=x, y=x, grid=TRUE,model="ci",param=c(0,1,0,1),me="ad");
dev.set(); image(z,zlim=c(-3,3));dev.set();

##print(c(mean(as.double(z)),var(as.double(z))))
cat("\n landmark 10\n")

dev.set()
m <- matrix(c(1,2,3,4),ncol=2)
z1 <- GaussRF(x=x, y=x, grid=TRUE,
              model=list(
                list(m="ci",v=1,a=m),
                "+", list(m="gau", v=1, a=m)
                ),
              me="add", reg=0,n=1)
print(c(mean(as.double(z1)),var(as.double(z1))))
image(z1,zlim=c(-3,3))
cat("\n landmark 11")
#}

dev.set()
m <- matrix(c(1,2,3,4),ncol=2)
z2 <- GaussRF(x=x, y=x, grid=TRUE,
              model=list(
                list(m="ci",v=1,a=m),
                "+", list(m="gau", v=1, a=m)
                ),
              me="ci", reg=1,n=1)
image(z2,zlim=c(-3,3))
print(c(mean(as.double(z1)),var(as.double(z1))))
print(c(mean(as.double(z2)),var(as.double(z2))))
cat("\n landmark 12\n")


## TBM3

dev.set()
m <- matrix(c(1,2,3,4),ncol=2)/5
z1 <- GaussRF(x=x, y=x, grid=TRUE,
              model=list(
                list(m="power",v=1,k=2,a=m),
                "*", list(m="sph", v=1, a=m)
                ),
              me="TBM3", reg=0,n=1)
print(c(mean(as.double(z1)),var(as.double(z1))))
str(z1)
image(z1,zlim=c(-3,3))

dev.set()
z2 <- GaussRF(x=x, y=x, grid=TRUE,
              model=list(
                list(m="power",v=1,k=2,a=m),  "*",
                list(m="sph", v=1, a=m)
                ),
              me="ci", reg=1)
image(z2,zlim=c(-3,3))
print(c(mean(as.double(z1)),var(as.double(z1))))
print(c(mean(as.double(z2)),var(as.double(z2))))

cat("\n landmark 12\n")
    
#stop("")


dev.set()
z1 <- GaussRF(x=x, y=x, grid=TRUE,
              model=list(
                list(m="power",v=1,k=2,s=2),
                "*", list(m="sph", v=1, s=3)
                ),
              me="TBM3", reg=0,trend=0,n=1)
print(c(mean(as.double(z1)),var(as.double(z1))))
image(z1,zlim=c(-3,3))

dev.set()
z2 <- GaussRF(x=x, y=x, grid=TRUE,
              model=list(
                list(m="power",v=1,k=2,s=2),  "*",
                list(m="sph", v=1, s=3)
                ),
              me="ci", reg=1)
image(z2,zlim=c(-3,3))
print(c(mean(as.double(z1)),var(as.double(z1))))
print(c(mean(as.double(z2)),var(as.double(z2))))
cat("\n landmark 13\n")


dev.set()
z1 <- GaussRF(x=x, y=x, grid=TRUE, model="gauss", param=c(0,1,0,1),
              me="TBM3",reg=0)
image(z1,zlim=c(-3,3))

dev.set()
z1 <- GaussRF(x=x, y=x, grid=TRUE, model=list(list(m="gauss", v=1,s=2)),
              me="ci",reg=1)
image(z1,zlim=c(-3,3))

dev.set()
z1 <- GaussRF(x=x, y=x, grid=TRUE, model=list(list(m="gauss", v=1,s=1)),
              me="TBM3",reg=0)
image(z1,zlim=c(-3,3))
#stop("")

cat("\n landmark 14\n")

## TBM2
z1 <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="power",v=10,k=2,a=matrix(c(0.5,0,1,0.5),
                                                  ncol=2)),
               "+",
               list(m="sph", v=10, a = matrix(c(0.5,-1,0,0.5),ncol=2))
               ),
              me="TBM2", TBM2.num=TRUE)
print(c(mean(as.double(z1)),var(as.double(z1))))

image(z1,zlim=c(-3,3))

cat("\n landmark 15\n")

z1 <- GaussRF(x=x, y=x, grid=TRUE, model=list(list(m="power", v=1,k=2,s=10)),
              me="TBM2", TBM2.num=TRUE)
image(z1,zlim=c(-3,3))

## spectral TBM
dev.set()
z <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="gauss", v=1, a=matrix(c(0.5,0,1,0.5),ncol=2)),
               "+",
               list(m="wave", v=1, a = matrix(c(0.5,-1,0,0.5),ncol=2))
               ),me="sp")
image(z,zlim=c(-3,3)) 
dev.set()
cat("\n landmark 16\n")

z1 <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="gauss", v=1, a=matrix(c(0.5,0,1,0.5),ncol=2))),
              me="sp")
image(z1,zlim=c(-3,3))

z2 <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="wave", v=1, a = matrix(c(0.5,-1,0,0.5),ncol=2))
               ),
              me="sp")
image(z2,zlim=c(-3,3))
image(z1+z2,zlim=c(-3,3))


cat("\n landmark 17\n")
x <- (0:30)/10
n <- 100

## check if calling twice the same procedure (due to different anisotropies)
## works well
z <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="exponen", v=1, a=matrix(c(0.5,0,1,0.5),ncol=2)),
               "*",
               list(m="gneiting", v=1, a = matrix(c(0.5,-1,0,0.5),ncol=2)),
               "+",
               list(m="nugget", v=0.3, a=matrix(c(1,0,0,1),ncol=2))
               ,"+",
               list(m="nugget", v=1, a=matrix(c(0,0,1+runif(1),0),ncol=2))
               ))
image(z) #1

cat("\n landmark 18")

## direct
z <- GaussRF(x=x, y=x, grid=TRUE, model="exponen", param=c(0,1,0,4),
             method="direct")
image(z) 

z <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="exponen", v=1, a=matrix(c(0.5,0,1,0.5),ncol=2)),
               "*",
               list(m="gneiting", v=1, a = matrix(c(0.5,-1,0,0.5),ncol=2)),
               "+",
               ##list(m="nugget", v=0.3, a=matrix(c(1,0,0,1),ncol=2))
               list(m="nugget", v=1, a=matrix(c(0,0,1,0),ncol=2))
               ), meth="direct")
image(z) #2

cat("\n landmark 19\n")


z <- GaussRF(x=x, y=x, grid=TRUE, model="exponen", param=c(0,1,0,4))
image(z) #3

cat("\n landmark 20\n")

#######################################################################

cat("\n landmark 21-0 \n")

x <- (0:100)/10
x <- (0:30)/10
n <- 100

RFparameters(CE.trials=3, Print=2); DeleteAllRegisters()
z <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="exponen", v=1, a=matrix(c(0.5,0,1,0.5),ncol=2)))
               )
image(z) #4

cat("\n landmark 21\n")

## unspecified, general testing

RFparameters(Print=1)
z <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="exponen", v=1, a=matrix(c(0.5,0,1,0.5),ncol=2)),
               "*",
               list(m="gneiting", v=1, a = matrix(c(0.5,-1,0,0.5),ncol=2)),
               "+",
               ##list(m="nugget", v=0.3, a=matrix(c(1,0,0,1),ncol=2))
               list(m="nugget", v=1, a=matrix(c(0,0,1,0),ncol=2))
               ))
image(z) #5


cat("\n landmark 22\n")


z <- GaussRF(x=x, y=x, grid=TRUE,
             model=list(list(m="exponen", v=1, a=matrix(c(0.5,0,0,0.5),ncol=2))))
image(z) 

z <- GaussRF(x=x, y=x, grid=TRUE, model="gneiting", param=c(0,1,0,2))
image(z) #6

cat("\n landmark 23\n")

#stop("")


## nested model
x <- seq(0,5,0.1)
p1 <- runif(2)
p2 <- runif(2, max=10)
p3 <- c(runif(1),0)
n1 <- CovarianceFct(x,"expo",param=c(NA,p1[1],0,p1[2]))
n2 <- CovarianceFct(x,"expo",param=c(NA,p2[1],0,p2[2])) 
n3 <- CovarianceFct(x,"nugget",param=c(NA,p3[1],0,1))

n <- CovarianceFct(x,"expo",param=cbind(p1,p2,p3))
sum(abs(n - n1 - n2 - n3))

cat("\n landmark 24\n")

## test on anisotropic model and arbitrary combinations
dd <- CovarianceFct(x,
                    list(list(m="whittle", k=3, v=5, s=4),
                         "*",
                         list(m="expo", v=4, s=10),
                         "+",
                         list(m="nugget", v=8)
                         ),dim=2)
cat("\n landmark 24a\n")

cc <- CovarianceFct(cbind(x,0),
                    list(list(m="whittle", k=3, v=5, a=diag(2)/4),
                         "*",
                         list(m="expo", v=4, a=diag(2)/10),
                         "+",
                         list(m="nugget", v=8,  a=diag(2))
                         ))
cat("\n landmark 25\n")

w <- CovarianceFct(x, "whittle", param=c(NA,5,0,4,3))
e <- CovarianceFct(x,"exponen", param=c(NA,4,0,10))
n <- CovarianceFct(x,"nugget", param=c(NA,8,0,0))

print(cc)
print(dd) 
print(w * e + n)
print(w * e + n -cc)
print(w * e + n -dd)
print(w)
print(e)
#stop("")

gneitingdiff <-  function(p){
  list(list(m="gneiting", v=p[2], s= p[6] * p[4] * 10 * sqrt(2) / 47),
       "*",
       list(m="whittle", k=p[5], v=1.0, s=p[4]),
       "+",
       list(m="nugget", v=p[3]))
}

x <- seq(0,5,len=10)
param <- c(NA, 3, 0.1, 2, 5, 6)
if (exists("PrepareModel"))
  CovarianceFct(x,gneitingdiff(param)) else
CovarianceFct(x, "gneitingdiff", param=param)

cat("\n landmark 26")


str(PrepareModel(list(list(m="whittle", k=3, v=5, s=4),
                      "*",
                      list(m="expo", v=4, s=10),
                      "+",
                      list(m="nugget", v=8)
                      ), timespacedim=1))

str(PrepareModel(list(list(m="whittle", k=3, v=5, a=matrix(1:4,ncol=2)),
                      "*",
                      list(m="expo", v=4, a=matrix(11:14,ncol=2)),
                      "+",
                      list(m="nugget", v=8, a=c(0,0,0,1))
                      ), timespacedim=2))


str(PrepareModel( list(list(m="whittle", k=3, v=5, s=4),
                   "*",
                   list(m="expo", v=4, s=10),
                   "+",
                   list(m="nugget", v=8)
                   ), timespacedim=1))

str(PrepareModel("whittle", cbind(c(1,2,3), c(0.5,0,9), c(5,3,7)),
                 timespacedim=1))

str(PrepareModel("expon",c(NA,1,2,3), timespacedim=1))

str(PrepareModel("expon",as.matrix(c(1,3)), timespacedim=1))

str(PrepareModel("expon", cbind(c(1,2), c(0.5,0), c(5,3)), timespacedim=1))

        
RFparameters(Print=5, Storing=FALSE)
GaussRF(1:10, 1:10, model=list(list(model="exp", s=5, var=1),
                      "+", list(model="gauss", s=5, var=1)), grid=TRUE)
