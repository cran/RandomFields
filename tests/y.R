# source("y.R")

if (EXTENDED.TESTING <- file.exists("source.R")) source("source.R")

if (FALSE) {
 model <-"gencauchy"
 param <- c(0, 1, 0, 1, 1, 2)
 estparam <- c(0, NA, 0, NA, NA, 2) ## NA means: "to be estimated"
 ## sequence in `estparam' is
 ## mean, variance, nugget, scale, (+ further model parameters)
 ## So, mean, variance, and scale will be estimated here.
 ## Nugget is fixed and equals zero.
 points <- 100
 x <- runif(points,0,3)
 y <- runif(points,0,3) ## 100 random points in square [0, 3]^2
 d <- GaussRF(x=x, y=y, grid=FALSE, model=model, param=param, n=10)



estparam <- c(0, NA, NA, NA, NA, NaN) 
parampositions(model=model, param=estparam)
f <- function(param) {param[3] <- param[2]; param}
z <- unlist(fitvario(x=cbind(x,y), data=d, model=model,
      param=estparam, transform=f,
      cross.me=NULL)$values)/unlist(fitvario(x=cbind(x,y), data=d, model=model,
      param=estparam, transform=f, cross.me=NULL,
      standard.style=TRUE)$values) -1
z [abs(z) < 1e-7] <- 0
z

# warum funktioniert 2te methode bei mle besser und bei lsq schlechter?!
 ###
 
}
