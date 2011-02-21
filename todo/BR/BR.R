
MaxStableRF.BR<- function (x, y = NULL, z = NULL, grid, model, param, maxstable,
                         method = NULL, n = 1, register = 0,
                         gridtriple = FALSE,
                         ...) {
  if (maxstable=="BrownResnick") {
      stopifnot(grid==TRUE)
      stopifnot(gridtriple==TRUE)
      step <- x[3]
      dim <- 4-is.null(y)-is.null(z)-is.null(T)
      eval(parse(text=paste("ab <- cbind(", paste(c("x","y","z","T")[1:dim],
                      collapse=","), ")",collapse="")))
      a <- ab[1,]
      b <- ab[2,]
      method <- as.integer(method)
      stopifnot(0 <= method && method <= 4)
      return(eval(parse(text=paste("meth",method,"(a=a,b=b,step=step,
                                   model=model,param=param,...)", sep=""))))
  } else {
    old.param <- RFparameters(no.readonly=TRUE)
    RFpar <- list(...)
    if (length(RFpar)>0) RFparameters(RFpar)
    if (delete <- n>1 && !RFparameters()$Storing) RFparameters(Storing=TRUE)
    on.exit({if (delete) DeleteRegister(register);
           RFparameters(old.param);
           })
    error <- InitMaxStableRF(x=x, y=y, z=z, grid=grid, model=model, param=param,
                           maxstable=maxstable, method=method,
                           register=register, gridtriple=gridtriple)
    if (error > 0) stop("InitMaxStable: error ", error, " occured")

    #print(delete); xxxx
  
    return(if (n<1) NULL else DoSimulateRF(n=n, reg=register, paired=FALSE))
  }
}


## * Struktur uebernehmen (z.B. nicht in MaxStableRF reinmogeln)
## * nach C runterschieben
## * vereinfachen, u.a do.call statt eval(...)
## * splitten in Init und Do
## * call-struktur identisch zu init und do
## * .C vermeiden, falls wrapper existiert, z.B. Variogram
## * method: ist die fuer das Gauss-Feld
##       * wahl ueber extra RFparameter
##       * automatic choice als default
## * .Rd schreiben
## emacs benutzen !
## make c : R CMD check RamdomFields
## local.dependent geschickt nutzen?

## BR methods

meth0 <- function(a, b, step, model, param, number, kmax, 
                  output=TRUE) {

   eucnorm <- function(x) sum(x^2)
   stopifnot(a<b)
   shift <- 0.5*(b+a)
   a <- a - shift
   b <- b - shift
   tdim <- length(a)
   param <- if (!missing(param)) param else NULL

   t <- t(expand.grid(seq(a[1],b[1],step),
         if(tdim > 1) seq(a[2],b[2],step) else 0,
         if(tdim > 2) seq(a[3],b[3],step) else 0,
         if(tdim > 3) seq(a[4],b[4],step) else 0))
   t <- t[1:tdim,,drop=FALSE]
   storage.mode(t) <- "double"
   lenges <- ncol(t)
   len <- apply(t, 1, function(x) length(unique(x)))
   .Call("CheckModelUser", PrepareModel(model, param)$model,
               as.integer(tdim), as.integer(tdim), FALSE, PACKAGE="RandomFields")
   trend <- double(lenges)
   .C("Variogram", t, NULL, as.integer(lenges), trend, PACKAGE="RandomFields",
               DUP = FALSE, NAOK=TRUE)
   dim(trend) <- len
   GaussRF(x=c(a[1],b[1],step), y=if(tdim > 1) c(a[2],b[2],step) else NULL,
           z=if(tdim > 2) c(a[3],b[3],step) else NULL,
           T=if(tdim > 3) c(a[4],b[4],step) else NULL,
           grid=TRUE, gridtriple=TRUE,
           model=model, param=param, n=0, Storing=TRUE)
   nullindex <- which(apply(t,2,eucnorm) == min(apply(t,2,eucnorm))) -1
   result <- matrix(-1e10, nrow=lenges, ncol=number)
   poisson <- rep(0,times=number)
   .C("meth0", as.double(poisson), as.integer(number),
          as.double(trend), as.integer(lenges), as.integer(kmax),
          as.integer(nullindex), result, as.integer(output),
          PACKAGE="RandomFields", DUP = FALSE)
   return(result)
}

meth1 <- function(a, b, step, model, param, h, number, kmax,
                  output=TRUE) {
   eucnorm <- function(x) sum(x^2)
   stopifnot(a<b)
   shift <- 0.5*(b+a)
   a <- a - shift
   b <- b - shift
   tdim <- length(a)
   h <- as.matrix(h-shift, nrow=tdim)
   n <- ncol(h)
   hmin <- apply(h,1,min)
   hmax <- apply(h,1,max)

   param <- if (!missing(param)) param else NULL
   t <- t(expand.grid(seq(min(a[1],hmin[1]),max(b[1],hmax[1]), step),
         if(tdim > 1) seq(min(a[2],hmin[2]),max(b[2],hmax[2]), step) else 0,
         if(tdim > 2) seq(min(a[3],hmin[3]),max(b[3],hmax[3]), step) else 0,
         if(tdim > 3) seq(min(a[4],hmin[4]),max(b[4],hmax[4]), step) else 0))
   t <- t[1:tdim,,drop=FALSE]
   proclen <- apply(t,1, function(x) length(unique(x)))
   trendt <- t(expand.grid(seq(min(t[1,1],a[1]-hmax[1]),
                               max(t[1,prod(proclen)],b[1]-hmin[1]), step),
                if(tdim>1) seq(min(t[2,1],a[2]-hmax[2]),
                               max(t[2,prod(proclen)],b[2]-hmin[2]), step) else 0,
                if(tdim>2) seq(min(t[3,1],a[3]-hmax[3]),
                               max(t[3,prod(proclen)],b[3]-hmin[3]), step) else 0,
                if(tdim>3) seq(min(t[4,1],a[4]-hmax[4]),
                               max(t[4,prod(proclen)],b[4]-hmin[4]), step) else 0))
   trendt <- trendt[1:tdim,,drop=FALSE]
   trendlen <- apply(trendt,1, function(x) length(unique(x)))
   storage.mode(trendt) <- "double"
   .Call("CheckModelUser", PrepareModel(model, param)$model,
               as.integer(tdim), as.integer(tdim), FALSE, PACKAGE="RandomFields")
   trend <- double(prod(trendlen))
   .C("Variogram", trendt, NULL, as.integer(prod(trendlen)), trend,
            PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)
   GaussRF(x=c(t[1,1], t[1,prod(proclen)], step),
           y=if(tdim >1) c(t[2,1], t[2,prod(proclen)], step) else NULL,
           z=if(tdim >2) c(t[3,1], t[3,prod(proclen)], step) else NULL,
           T=if(tdim >3) c(t[4,1], t[4,prod(proclen)], step) else NULL,
           grid=TRUE, gridtriple=TRUE, model=model, param=param, n=0)
   nullindex <- which(apply(t,2,eucnorm) == min(apply(t,2,eucnorm))) -1
   hindex <- apply(h, 2, function(h) which(apply(t-h,2,eucnorm) 
                                     == min(apply(t-h,2,eucnorm)))) - 1
   startindex <- which(apply(t-a,2,eucnorm) == min(apply(t-a,2,eucnorm))) - 1
   trendnullindex <- which(apply(trendt,2,eucnorm) == min(apply(trendt,2,eucnorm))) -1
   trendhindex <- apply(h, 2, function(h) which(apply(trendt-h,2,eucnorm) 
                                     == min(apply(trendt-h,2,eucnorm)))) - 1
   trendstartindex <- which(apply(trendt-a,2,eucnorm) == min(apply(trendt-a,2,eucnorm))) - 1
   reslen <- floor((b-a)/step) + 1
   result <- matrix(-1e10, nrow=prod(reslen), ncol=number)
   poisson <- rep(0,times=number*n)
   .C("meth1", as.double(poisson), as.integer(number),
        as.double(trend), as.integer(tdim), 
        as.integer(n), as.integer(kmax),
        as.integer(nullindex), as.integer(trendnullindex),
        as.integer(hindex), as.integer(trendhindex),
        as.integer(startindex), as.integer(trendstartindex), 
        as.integer(reslen), as.integer(proclen), as.integer(trendlen),
        result, as.integer(output), PACKAGE="RandomFields", DUP = FALSE)
   return(result)
}

meth2 <- function(a, b, step, model, param, lower, upper, number,
                  kmax, output=TRUE) {
   eucnorm <- function(x) sum(x^2)
   stopifnot(a<b)
   stopifnot(lower<upper)
   shift <- 0.5*(b+a)
   a <- a - shift
   b <- b - shift
   lower <- lower - shift
   upper <- upper - shift
   param <- if (!missing(param)) param else NULL
   tdim <- length(a)

   ## use: as.list(as.data.frame(cbind(a, b)))
   
   t <- t(expand.grid(seq(min(a[1],lower[1]),max(b[1],upper[1]), step),
         if(tdim > 1) seq(min(a[2],lower[2]),max(b[2],upper[2]), step) else 0,
         if(tdim > 2) seq(min(a[3],lower[3]),max(b[3],upper[3]), step) else 0,
         if(tdim > 3) seq(min(a[4],lower[4]),max(b[4],upper[4]), step) else 0))
   t <- t[1:tdim,,drop=FALSE]
   proclen <- apply(t,1,function(x) length(unique(x)))
   trendt <- t(expand.grid(seq(min(t[1,1], a[1]-upper[1]),
                               max(t[1,prod(proclen)], b[1]-lower[1])),
              if(tdim > 1) seq(min(t[2,1], a[2]-upper[2]),
                               max(t[2,prod(proclen)], b[2]-lower[2])) else 0,
              if(tdim > 2) seq(min(t[3,1], a[3]-upper[3]),
                               max(t[3,prod(proclen)], b[3]-lower[3])) else 0,
              if(tdim > 3) seq(min(t[4,1], a[4]-upper[4]),
                               max(t[4,prod(proclen)], b[4]-lower[4])) else 0))
  trendt <- t[1:tdim,,drop=FALSE]
  trendlen <- apply(trendt,1,function(x) length(unique(x)))
  .Call("CheckModelUser", PrepareModel(model, param)$model,
               as.integer(tdim), as.integer(tdim), FALSE, PACKAGE="RandomFields")
  trend <- double(prod(trendlen))
  .C("Variogram", trendt, NULL, as.integer(trendlen), trend,
            PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)
   GaussRF(x=c(t[1,1], t[1,prod(proclen)], step),
           y=if(tdim >1) c(t[2,1], t[2,prod(proclen)], step) else NULL,
           z=if(tdim >2) c(t[3,1], t[3,prod(proclen)], step) else NULL,
           T=if(tdim >3) c(t[4,1], t[4,prod(proclen)], step) else NULL,
           grid=TRUE, gridtriple=TRUE, model=model, param=param, n=0)
   nullindex <- which(apply(t,2,eucnorm) == min(apply(t,2,eucnorm))) -1
   startindex <- which(apply(t-a,2,eucnorm) == min(apply(t-a,2,eucnorm))) - 1
   trendnullindex <- which(apply(trendt,2,eucnorm) == min(apply(trendt,2,eucnorm))) -1
   trendstartindex <- which(apply(trendt-a,2,eucnorm) == min(apply(trendt-a,2,eucnorm))) - 1
   reslen <- floor((b-a)/step)+1
   result <- matrix(-1e10, nrow=prod(reslen), ncol=number)
   poisson <- rep(0,times=number)
   .C("meth2", as.double(poisson), as.double(t), as.integer(tdim), as.integer(number),
          as.double(trend), as.double(trendt), as.integer(proclen), as.integer(trendlen),
          as.integer(kmax), as.integer(trendnullindex),
          as.integer(startindex), as.integer(trendstartindex), as.double(lower),
          as.double(upper), as.integer(reslen), result, as.integer(output),
          PACKAGE="RandomFields", DUP = FALSE)
   return(result)
}

meth3 <- function(a, b, step, model, param, jmin, jmax, number, m=1,
                   kmax, output=TRUE) {

   eucnorm <- function(x) sum(x^2)
   stopifnot(a<=b && jmin <= jmax)
   if (a<jmin*step || jmax*step<b) 
                      warning("bounds may be too narrow")
   param <- if (!missing(param)) param else NULL
   if (model!="fractalB")
       warning("Method3 may not work with this model!")
   if (m%%2 == 0) m <- abs(m)+1
      #aus Symmetrigruenden ist nur m ungerade bei Diskretisierung
      # sinnvoll, denn mit p=step heisst - m*p/2 < argsup <= m*p/2
      # (kontinuierlich), dass |argsup| <= m%/%2 * p (diskret),
      # falls m ungerade
   shift <- 0.5*(b+a)
   a <- a - shift
   b <- b - shift
   jmin <- jmin - sapply(shift/step,round)
   jmax <- jmax - sapply(shift/step,round)
   tdim <- length(a)
   n <- jmax - jmin +1
   t <- t(expand.grid(jmin[1]:jmax[1],
          if (tdim>1) jmin[2]:jmax[2] else 0,
          if (tdim>2) jmin[3]:jmax[3] else 0,
          if (tdim>3) jmin[4]:jmax[4] else 0))*step
   t <- t[1:tdim,,drop=FALSE]
   trendt <- t(expand.grid(seq(jmin[1]-jmax[1],jmax[1]-jmin[1]),
                if(tdim>1) seq(jmin[2]-jmax[2],jmax[2]-jmin[2]) else 0,
                if(tdim>2) seq(jmin[3]-jmax[3],jmax[3]-jmin[3]) else 0,
                if(tdim>3) seq(jmin[4]-jmax[4],jmax[4]-jmin[4]) else 0))*step
   trendt <- trendt[1:tdim,,drop=FALSE]
   storage.mode(trendt) <- "double"
   .Call("CheckModelUser", PrepareModel(model, param)$model,
               as.integer(tdim), as.integer(tdim), FALSE, PACKAGE="RandomFields")
   trend <- double(prod(2*n-1))
   .C("Variogram", trendt, NULL, as.integer(prod(2*n-1)), trend,
               PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)
   GaussRF(x=c(t[1,1], t[1,prod(n)], step),
           y=if(tdim >1) c(t[2,1], t[2,prod(n)], step) else NULL,
           z=if(tdim >2) c(t[3,1], t[3,prod(n)], step) else NULL,
           T=if(tdim >3) c(t[4,1], t[4,prod(n)], step) else NULL,
           grid=TRUE, gridtriple=TRUE, model=model, param=param, n=0)
   trendnullindex <- which(apply(trendt,2,eucnorm) == min(apply(trendt,2,eucnorm))) -1
   startindex <- which(apply(t-a,2,eucnorm) == min(apply(t-a,2,eucnorm))) -1
   reslen <- floor((b-a)/step) +1
   result <- matrix(-1e10, nrow=prod(reslen), ncol=number)
   poisson <- rep(0,times=prod(n)*number)
   .C("meth3", as.double(poisson), as.integer(number),
               as.double(trend), as.integer(n), as.integer(tdim),
               as.integer(kmax), as.integer(trendnullindex), as.integer(startindex),
               as.integer(m), as.integer(reslen), result,
               as.integer(output), PACKAGE="RandomFields", DUP = FALSE)
   return(result)
}

meth4 <- function(a, b, step, model, param, lower, upper, number,
                  lambda, kmax, output=TRUE) {

   eucnorm <- function(x) sum(x^2)
   stopifnot(a<=b && lower <= upper)
   if (a<lower || upper<b) warning("bounds my be too narrow")
   shift <- (a+b)/2
   a <- a - shift
   b <- b - shift
   lower <- lower - shift
   upper <- upper - shift
   param <- if (!missing(param)) param else NULL
   if (model!="fractalB")
       warning("Method4 may not work with this model!")
   tdim <- length(a)
   t <- t(expand.grid(seq(a[1]+lower[1], b[1]+upper[1], step),
          if (tdim>1) seq(a[2]+lower[2], b[2]+upper[2], step) else 0,
          if (tdim>2) seq(a[3]+lower[3], b[3]+upper[3], step) else 0,
          if (tdim>3) seq(a[4]+lower[4], b[4]+upper[4], step) else 0))
   t <- t[1:tdim,,drop=FALSE]
   proclen <- apply(t,1,function(x) length(unique(x)))
   trendt <- t(expand.grid(seq(a[1]+lower[1]-upper[1], b[1]+upper[1]-lower[1], step),
              if (tdim>1) seq(a[2]+lower[2]-upper[2], b[2]+upper[2]-lower[2], step) else 0,
              if (tdim>2) seq(a[3]+lower[3]-upper[3], b[3]+upper[3]-lower[3], step) else 0,
              if (tdim>3) seq(a[4]+lower[4]-upper[4], b[4]+upper[4]-lower[4], step) else 0))
   trendt <- trendt[1:tdim,,drop=FALSE]
   trendlen <- apply(trendt,1,function(x) length(unique(x)))
   storage.mode(trendt) <- "double"
   .Call("CheckModelUser", PrepareModel(model, param)$model,
               as.integer(tdim), as.integer(tdim), FALSE, PACKAGE="RandomFields")
   trend <- double(prod(trendlen))
   .C("Variogram", trendt, NULL, as.integer(prod(trendlen)), trend,
             PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)
   GaussRF(x=c(t[1,1], t[1,prod(proclen)], step),
           y=if(tdim >1) c(t[2,1], t[2,prod(proclen)], step) else NULL,
           z=if(tdim >2) c(t[3,1], t[3,prod(proclen)], step) else NULL,
           T=if(tdim >3) c(t[4,1], t[4,prod(proclen)], step) else NULL,
           grid=TRUE, gridtriple=TRUE, model=model, param=param, n=0)
   nullindex <- apply(t,1, function(x) which(abs(unique(x))==min(abs(unique(x)))))
   lowerindex <- apply(t-lower,1, function(x) which(abs(unique(x))==min(abs(unique(x)))))
   startindex <- which(apply(t-a,2,eucnorm) == min(apply(t-a,2,eucnorm))) -1
   trendprocstartindex <- which(apply(trendt-a-lower,2,eucnorm)
                             == min(apply(trendt-a-lower,2,eucnorm))) -1
   intlen <- floor((upper-lower)/step) + 1
   reslen <- floor((b-a)/step) + 1
   result <- matrix(-1e10, nrow=prod(reslen), ncol=number)
   poisson <- rep(0,times=number)
   .C("meth4", as.double(lambda), as.double(poisson),
               as.integer(number), as.integer(tdim), as.double(trend),
               as.double(t), as.integer(proclen), as.integer(nullindex),
               as.integer(startindex), as.integer(trendprocstartindex),
               as.integer(trendlen), as.integer(lowerindex),
               as.double(lower), as.double(upper),
               as.integer(intlen), as.integer(reslen),
               as.integer(kmax), result, as.integer(output),
               PACKAGE="RandomFields", DUP = FALSE)
   return(result)
}

findlambda <- function(step, model, param, lower, upper, number,
                        kmax, output=TRUE) {

  param <- if (!missing(param)) param else NULL
  stopifnot(lower <0 && upper>0)

  result <- meth4(0, 0, step, model=model, param=param, lower=lower,
                  upper=upper, number=number, lambda=1, kmax=kmax,
                  output=output)
  lambda <- exp(0.557-mean(result))

}





