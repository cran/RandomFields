

##library(foreign)
##library(aux)
source("data/ismev.dat")

theta.hat <- function(Z, instances=length(Z)/4,
                      threshold=quantile(Z,p,na.rm=TRUE),
                      p=0.95, method="simple") {
  l.Z <- length(Z)
  theta <- rep(NA,instances)
  switch(pmatch(method, c("simple", "HT")),
         for (i in 1:instances) {
           X <- pmax(Z[-(1:i)], Z[1 : (l.Z - i)])
           K <- sum(X < threshold)
           S <- sum( 1 / X[X>=threshold] )
           theta[i] <- (length(X) - K) / (K/threshold + S)
         },
         {
           stop(" HT noch nicht vollstaendig programmiert")
           for (i in 1:instances) {
             X <- pmax(Z[-(1:i)], Z[1 : (l.Z - i)])
             theta[i] <-
             sum(X > threshold) / sum(1/pmax(threshold, X)) 
           }
         }
       )
  return(theta)
}


MMA <- function(n, kernel, normalise=TRUE) {
  if (normalise) kernel <- kernel / sum(kernel)
  l.k <- length(kernel)
  tot.l <- n + l.k - 1
  frechet <- -1/log(runif(tot.l))
  old.warn <- options()$warn
  options(warn=-1)
  m <- t(matrix(c(frechet,NA), nrow=tot.l, ncol=l.k))
  options(warn=old.warn)
  m.d <- apply(m * kernel, 2, max)              
  m.r <- apply(m * (-log(1-runif(tot.l * l.k))) *
               matrix(kernel,  nrow=l.k, ncol=tot.l), 2, max)
  return(list(uiv=frechet[n : 1], mm=m.d[tot.l : l.k],
              r=m.r[tot.l : l.k]))
}


MMA.model <- function(kernel, normalise=TRUE) {
  if (normalise) kernel <- kernel / sum(kernel)
  result <- rep(0,length(kernel))
  for (i in 1:length(kernel)) {
    result[i] <-
      sum(pmax(c(rep(NA, i-1), kernel), c(kernel, rep(NA, i-1)),na.rm=TRUE))
  }
  return(result)
}

ts.analyse <- function(l=length(data), kernel, instances, n=plotn, data, p,
                         plotn=0, lty=1, ylim=c(1,2)){
  text <- "Modell"
  col <- "blue"
  m.th <- s.th <- d.th <- NULL
  plot(0, 0, xlim=c(0,instances), ylim=ylim, xlab="t", ylab=expression(theta))
  if (n>0) {
    th <-  matrix(ncol=instances, nrow=n)
    text <- c(text, "Simu")
    col <- c(col, "grey")
    for (i in 1:n) {
      cat(i," ")
      ms <- MMA(l, kernel)$m
      th[i,] <- theta.hat(ms, instances, p=p)
      if (i<plotn) lines(1:instances, th[i,], col="grey", lty=lty)
    }
    s.th <- apply(th, 2, mean, na.rm=TRUE)
    text <- c(text, "Simu.mean")
    col <- c(col, "darkgreen")
    lines(1:instances, s.th, col="darkgreen", lty=lty)
    if (plotn<0) {
      lines(1:instances, apply(th, 2, max, na.rm=TRUE), col="brown", lty=lty)
      lines(1:instances, apply(th, 2, min, na.rm=TRUE), col="brown", lty=lty)
    } 
  }
  if (!(missing(data)||is.null(data))) {
    text <- c(text, "Daten")
    col <- c(col, "red")
    d.th <- theta.hat(data, instances, p=p)
    lines(1:instances, d.th, col="red",lty=lty)
  }
  m.th <- MMA.model(kernel)
  print(m.th)
  points(0:(length(kernel)-1), m.th, col="blue")
  legend(instances, 1, legend=text, xj=1, yj=0, lty=lty, col=col)
  return(list(m.th=m.th, s.th=s.th, d.th=d.th))
}

emp.mean.excess <-  function(data, max.p) {
  y <- seq(min(data), quantile(data,max.p), l=1000)
  n <- e <- rep(NA, length(y))
  for (i in 1:length(y)){
    over.threshold <- data>y[i]
    n[i] <- sum(over.threshold)
    e[i] <- sum(data[over.threshold]) /n[i]  - y[i]
  }
  plot(y,n/length(data))
  dev.set()
  plot(y,e)
  return(list(y=y,e=e))
}

to.frechet <- function(x) {
  l.x <- length(x)
  new <- -1/log( (order(order(x))-0.5) / l.x)
  stopifnot(all(is.finite(new)))
  sort.x <- sort(x)
  stopifnot(all(is.finite(sort.x)))  
  empirical <<- function(d) {
    result <- list()
    if (no.list <- !is.list(d)) d <- list(d)
    for (i in 1:length(d)) {
      u <- exp(-1/d[[i]]) * l.x + 0.5
      stopifnot(all(is.finite(u)))
      result[[i]] <- rep(0, length(d[[i]]))
      result[[i]][u>=l.x] <- sort.x[l.x]
      index <- (u>=1) & (u<l.x)
      ##print(index)
      ##print(u)
      u.int <- as.integer(u[index])
      result[[i]][index] <-
        sort.x[u.int] + (sort.x[u.int+1]-sort.x[u.int]) * (u[index]-u.int)
    }
    if (no.list) return(result[[1]]) else return(result)
  }
  return(new)
}
 
