# source("RFtest.mle.R")

if (file.exists("source.R")) source("source.R")

##load("random seed")
#
runif(1);save(file="random.seed",.Random.seed)
#load("random.seed")

model <- "expon"; standardparam <- c(0,1,0,1); factor <- 5; low <- up <- NULL
#model<-"whittle"; standardparam<-c(0,1,0,1,1); factor<-5; low<-0.1; up<-20

use <- TRUE ##use <- FALSE
PrintLevel <- 1
points <- 200 
repet <- 1:20   ## 500 if being serious
same <-  10;    ## 20, 500

####### for heavy testing, uncomment the following 3 lines
# repet <- 1:30;  model<-"whittle";param<-c(0,1,0,0.1,3);
# estparam <- c(NA,NA,0,NA,NA); low<-0.2;up<-10; points <- 300
# PrintLevel <- 5


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
  
  #######estparam[c(3,4)] <- NA
  
  estparam[runif(length(param))<0.5] <-  NA
                     ## NA means to be estimated,
                     ## sequ mean,variance,nugget,scale,further parameters
                     ## (as usual in RandomFields)

  # estparam <- param;estparam[1] <- NA;
  
  if (all(is.na(estparam[c(2,3)])) && (runif(1)<0.5)) {
    sill <- sum(param[c(2,3)])
  } else sill <- NA
   
  cat(formatC(i,width=3),
      "est:",formatC(round(estparam, dig=2), dig=2, width=4),
      " sill=", ifelse(is.na(sill),"NA ", formatC(sill,dig=2,width=3)))
  
  ut <- unix.time( result <- (fitvario(x=cbind(x,y),data=f,model=model,
                                       param=estparam, sill=sill,
                                       lower=low, upper=up,
                                       PrintLevel=PrintLevel,use.nat=use,
                                       cross.m=NULL)))
  res <- formatC(round(result$vario$ml$param,dig=2),dig=2,width=5)
  res[!is.na(estparam)] <- "  . "
  cat("  mle:",res," [",format(result$value$ml,dig=2),
      "] (time:", format(ut[c(1)]),")\n")
  if (any(is.na( result$vario$ml$param))) sfasdf
  est[is.na(estparam)] <-
    est[is.na(estparam)] + result$vario$ml$param[is.na(estparam)]
  if (any(is.na( est))) sfasdffgcc
  nest[is.na(estparam)] <- nest[is.na(estparam)] + 1  
}
cat(  " true  ",formatC(param,width=5),
    "\naverage",formatC(est / nest,width=5),"\n")


