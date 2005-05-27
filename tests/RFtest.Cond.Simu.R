# source("RFtest.Cond.Simu.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

## normal checking
#nn <- c(10,1); step <- 1; #graphics.off(); X11(); X11();

## light; for automatic checking
nn <- c(1,3); step <- 1

p <- 1:7            ## step, p, x must be chosen such that all the 
x <-  seq(0,7,step) ## point given by p are in x !!
dims <- c(2:1)

## for 3d:
## dims <- 3; p <- (0:3)*2+1; 

zlim <- c(-2.6,2.6)
colour <- rainbow(100)
param <- c(0,1,0,1)
model <- "exponential"
RFparameters(PracticalRange=FALSE)
total <- 0; ## total error


add.given <- list(NULL,rbind(c(1.1,3.1,3.7),c(2.7,2.9,2.5))) ## null: grid,
                                                     ## arbitrary otherwise
add <- list(NULL,rbind(c(2.1,2.1,5.3),c(3.7,4.9,4.1)))       ## null: grid
krige.meth <- c("O", "S")

debug <- FALSE #debug <- TRUE
if (debug) { ## not for the public
  add.given <- list(rbind(c(1.1,3.1,3.7),c(2.7,2.9,2.5)))#NULL)
  add <- list(rbind(c(2.1,2.1,5.3),c(3.7,4.9,4.1)))
  nn <- 10
  dims <- 1
  meth <- "O"
}

for (krige.method in krige.meth) {
  for (dimension in dims) {
    for (n in nn) {
      for (add.points in add) {
        for (add.given.points in add.given) {
          cat("given=",length(add.given.points),
              " add=",length(add.points),
              " n=",n,
              " dim=",dimension,
              " krige.method=",krige.method
              )
          nugget <- max(0,runif(1,-1,1))
          switch(dimension,
                 { # 1 dim
            points <- c(p,add.given.points[,1])
            xp <- c(x,add.points[,1])
            index <- apply(as.matrix(p), 1,
                           function(k) (1:length(x))[abs(k-x)<0.0001])
          },{ #2 d
            pp <- as.matrix(expand.grid(p,p))
            points <- rbind(pp, add.given.points[,c(1,2)])
            xx <- rbind(as.matrix(expand.grid(x,x)), add.points[,c(1,2)])
            if (is.null(add.points)) {
              xp <- cbind(x,x)
            } else {
              xp <- xx
            }  
            index <-
              apply(pp, 1,
                    function(k) (1:nrow(xx))[apply(abs(k-t(xx)),2,
                                                   sum)<0.0001])
          },{ # 3 d
            pp <- as.matrix(expand.grid(p,p,p))
             points <- rbind(pp, add.given.points)
            xx <- rbind(as.matrix(expand.grid(x,x,x)), add.points)
            if (is.null(add.points)) {
              xp <- cbind(x,x,x)
            } else {
              xp <- xx
            }  
            index <-
              apply(pp, 1,
                    function(k) (1:nrow(xx))[apply(abs(k-t(xx)),2,
                                                   sum)<0.0001])
          }
                 ) # switch

          ##   cat(" generate and visualise spatial data\n")
          data <- GaussRF(points, grid=FALSE, model=model, param=param)
          
          ## conditional simulation 
          cz <- CondSimu(krige.method=krige.method,
                         x=xp, grid=is.null(add.points),
                         model=model, param=c(0,1,0,1),
                         given=points, data=data,n=n,
                         err.model="nugget",
                         err.param=c(0,0,nugget,0))
         
          dev.set(2)
          indexsum <- 0;
          switch(dimension,
                 { # 1d
            plot(p, data[1:length(p)], xlim=range(x), ylim=zlim)
            dev.set(3)
            for (i in 1:n) {
              if (n==1) y <- cz else y <- cz[,i]
              indexsum <- indexsum + sum(abs(y[index]-data[1:length(p)]))
              plot(xp, y, xlim=range(x), ylim=zlim)
            }
          }, { #2d
            image(p, p, xlim=range(x), ylim=range(x),
                  matrix(data[1:(length(p)^2)],ncol=length(p)),
                  col=colour,zlim=zlim)
            dev.set(3)
            for (i in 1:n) {
              if (n==1) zz <- matrix(cz[1:(length(x)^2)],ncol=length(x)) else
              if (is.null(add.points)) {
                zz <- cz[,,i]
              } else {
                zz <- matrix((cz[,i])[1:(length(x)^2)],ncol=length(x))
              }
              indexsum <- indexsum + sum(abs(zz[index]-data[1:(length(p)^2)]))
              image(x,x,zz,col=colour,zlim=c(-2.6,2.6))
            }
          }, { #3d
             for (i in 1:n) {
              if (n==1) zz <- cz[1:(length(x)^3)] else
              if (is.null(add.points)) {
                zz <- cz[,,,i]
              } else {
                zz <- (cz[,i])[1:(length(x)^3)]
              }
              indexsum <- indexsum + sum(abs(zz[index]-data[1:(length(p)^3)]))
            }
          })
          cat(" li=",sum(is.finite(index))," error=",indexsum,"\n")
          if (nugget==0) total <- total + indexsum
        }
      }
    }
  }
}
cat("total error=",total,"\n")











