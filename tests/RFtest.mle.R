# source("RFtest.mle.R")

library(RandomFields)

# source("../R/MLE.R")
##load("random seeed")
#runif(1);save(file="random.seed",.Random.seed)
#load("random.seed")


## this an old function, not programmed very nicely
## I do not care -- just used for debugging reasons
XMLEtarget<-function(model,param,coord,data) {
  covnr <- .C("GetModelNr",as.character(model),nr=integer(1))$nr;
  storage.mode(covnr) <- "integer"
  if (covnr<0) stop("model not ok")
  ip<-.C("GetParameterIndices", MEAN= integer(1), VARIANCE= integer(1),
          NUGGET= integer(1), SCALE= integer(1), 
          KAPPA= integer(1), LASTKAPPA= integer(1),integer(1),
          SILL= integer(1), DUP = FALSE)

  lc <- nrow(coord)
  storage.mode(lc) <- "integer"
  dim <- ncol(coord); if (dim!=2) stop("Sorry. dim<>2 not programmed.")
  storage.mode(dim) <- "integer"
  x <- dist(coord)
  storage.mode(x) <- "double"
  storage.mode(param) <- "double"
  cov.matrix <- .C("CovarianceMatrix",x,lc,covnr,param,
     as.integer(length(param)),dim,cov.matrix=double(lc * lc),
                   DUP=FALSE)$cov.matrix
  cov.matrix <- .Fortran("chol", x = cov.matrix, lc, lc,
                         v = matrix(0, nr = lc, nc = lc),
                         info = integer(1), DUP = FALSE, PACKAGE = "base")
  if (cov.matrix$info) { 
      print("trying svd")
      cov.matrix <- svd(matrix(cov.matrix$v,nrow=lc))
      det <- prod(cov.matrix$d)
      if (any(cov.matrix$d<0)) {print("det<0!"); det<-0}
      cov.matrix <- t(cov.matrix$v %*% (t(cov.matrix$u) /cov.matrix$d))  
    }
    else {
      det <- prod(diag(cov.matrix$v))
      cov.matrix <- chol2inv(cov.matrix$v); 
    }
  logdet <- 2.0*log(det);  if (!is.finite(logdet)) logdet<-sum(abs(variab)+1);

  data <- as.matrix(data)  
  res <- logdet+(t(data-param[ip$MEAN+1]) %*% cov.matrix %*%
                 (data-param[ip$MEAN+1]))
  return(res)
}


model <- "expon"; standardparam <- c(0,1,0,1); factor <- 5; low <- up <- NULL
model<-"whittle"; standardparam<-c(0,1,0,1,1); factor<-5; low<-0.1; up<-20

use <- TRUE ##use <- FALSE
PrintLevel <- 0
points <- 100  
repet <- 1:20; ## 500 if being serious
same <-  2;    ## 20

####### for heavy testing, uncomment the following 3 lines
# repet <- 1:30;  model<-"whittle";param<-c(0,1,0,0.1,3);
# estparam <- c(NA,NA,0,NA,NA); low<-0.2;up<-10; points <- 300
# PrintLevel <- 5
# source("../R/MLE.R") #

param <- 0
est <- nest <- 1

for (i in repet) {
########## define model + simulation
  x <- runif(points,0,3)
  y <- runif(points,0,3)

  repeat {
    if (i %% same == 1) {
      cat(  " true  ",formatC(param,width=5),
          "\naverage",formatC(est / nest,width=5),"\n")
      nest <- est <- rep(0,length(standardparam))
      switch(1+sum(runif(1)>c(0.2,1)),
             param <- standardparam,
             param <- runif(length(standardparam)) * factor
             )
      cat(" ******* ",format(param)," ***********\n")
    }
    if ((!is.null(low) && (low>param[-1:-4])) ||
        (!is.null(up) && (up<param[-1:-4]))) next;

    f <- GaussRF(x=x,y=y,grid=FALSE,model=model,param=param)
    if (!is.null(f)) break;
  }

  estparam <- param
  estparam[runif(length(param))<0.5] <-  NA
                     ## NA means to be estimated,
                     ## sequ mean,variance,nugget,scale,further parameters
                     ## (as usual in RandomFields)

  # estparam <- param; estparam[4] <- NA;
  
  if (all(is.na(estparam[c(2,3)])) && (runif(1)<0.5)) {
    sill <- sum(param[c(2,3)])
  } else sill <- NA
   
  cat(formatC(i,width=3),
      "est:",formatC(round(estparam,dig=2),dig=2,width=4),
      " sill=",ifelse(is.na(sill),"NA ",formatC(sill,dig=2,width=3)))
  
  
  ut <- unix.time( result <- (mleRF(coord=cbind(x,y),data=f,model=model,
                                    param=estparam, sill=sill,
                                    lower.kappa=low,upper.kappa=up,
                                    PrintLevel=PrintLevel,use=use)))
  res <- formatC(round(result,dig=2),dig=2,width=5)
  res[!is.na(estparam)] <- "  .  "
  cat(" mle:",res,
      " [",format(XMLEtarget(model,result,cbind(x,y),f),dig=2),
      "] (time:", format(ut[c(1)]),")\n")
  est[is.na(estparam)] <- est[is.na(estparam)] + result[is.na(estparam)]
  nest[is.na(estparam)] <- nest[is.na(estparam)] + 1
  
}




