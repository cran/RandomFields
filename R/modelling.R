
Kriging <- function(method, x, y=NULL, z=NULL, grid, gridtriple=FALSE,
                    model, param, given, data) {
  
  methodlist <- c("S","O")
  if (is.na(method.nr <- pmatch(method,methodlist)))
    stop(paste("kriging method not identifiable from the list.",
               "Possible values are",paste(methodlist,collapse=",")))
  ip<-.C("GetParameterIndices", MEAN= integer(1), VARIANCE= integer(1),
          NUGGET= integer(1), SCALE= integer(1), 
          KAPPA= integer(1), LASTKAPPA= integer(1),integer(1),
          SILL= integer(1), DUP = FALSE)
  ip <- lapply(ip,function(i){i+1})

  x <- CheckXYZ(x, y, z, grid, gridtriple)
  x <- cbind(x$x, x$y, x$z)
  y <- z <- NULL
  ncol.data <- ncol(data)
  is.matrix.data <- is.matrix(data) && (ncol.data>1)
  given <- as.matrix(given)
  nd <- nrow(given)
  covmatrix <- diag(nd) * (CovarianceFct(0,model,param) / 2)
  covmatrix[lower.tri(covmatrix)] <- CovarianceFct(dist(given),model,param)
  covmatrix <- covmatrix + t(covmatrix)

  tgiven <-  t(as.matrix(given))
  given <- NULL

  if (ncol(x)!=nrow(tgiven))
    stop("dimension of the kriged points and the given points do not match")
  if (grid) {
    dimension <- 1 + floor((x[2,]-x[1,])/x[3,])
    l <- length(dimension)
    if (l==1) x <- matrix(seq(x[1],x[2],x[3]),ncol=1)
    else {
     text <- paste("x <- expand.grid(",
                      paste("seq(x[1,",1:l,
                            "],x[2,",1:l,
                            "],x[3,",1:l,"])",collapse=","),
                      ")")
     eval(parse(text=text))
   }
  } 
   
  
  switch(method.nr,
         { ## simple kriging
           data <- as.matrix(solve(covmatrix,data-param[ip$MEAN]))
           res <- apply(x, 1,
                        function(z){CovarianceFct(sqrt(apply((tgiven-z)^2,2,
                                                             sum)),
                                                  model,param) %*% data
                                   }
                        ) + param[ip$MEAN]
         }, {
           ## ordinary kriging
           covmatrix <- rbind(cbind(covmatrix,1), c(rep(1,nd),0)) 
           data <- solve(covmatrix,rbind(as.matrix(data),0))
           res <-
             apply(x, 1,
                   function(z){c(CovarianceFct(sqrt(apply((tgiven-z)^2,2,sum)),
                                               model,param),1) %*% data
                             }
                   )
         }
         ) # switch
  x <- data <- NULL
  if (is.matrix.data) {
    if (grid) return(array(t(res),dim=c(dimension,ncol.data)))
    else return(t(res))
  } else {
    if (grid) return(array(res,dim=dimension))
    else return(res)
  }
}



CondSimu <- function(method, x, y=NULL, z=NULL, grid, gridtriple=FALSE,
                     n=1,
                     model, param, given, data, tol=1E-5) {
  x <- CheckXYZ(x, y, z, grid, gridtriple)
  y <- NULL   
  total <- x$total
  x <-  cbind(x$x, x$y, x$z)
  
  simu.grid <- grid
  given <- as.matrix(given)

  if (grid) {
    ind <- 1 + (t(given) - x[1,]) / x[3,]
    index <-  round(ind)
    dimension<- 1 + floor((x[2,]-x[1,])/x[3,])
    
    if (any(abs(ind-index)>tol) ||
        any(index<1) || any(index>dimension))
      {
      ## at least some data points are not on the grid
      ## simulate as if there is not grid
      simu.grid <- FALSE
      l <- length(dimension)
      if (l==1) xx <- matrix(seq(x[1],x[2],x[3]),nrow=1)
      else  eval(parse(text=paste("xx <-  t(as.matrix(expand.grid(",
                         paste("seq(x[1,",1:l,
                               "],x[2,",1:l,
                               "],x[3,",1:l,"])",collapse=","),
                         ")))")))
    } else {
      ## data points are all lying on the grid
      z <- GaussRF(x,model=model,param=param,grid=TRUE, n=n, gridtriple=TRUE)
      index <- t(index)
    }
  } else xx <- t(as.matrix(x)) ## not a grid
  
  if (!simu.grid) {
    one2ncol.xx <- 1:ncol(xx)
    tol <- tol * nrow(xx)
    index <- apply(as.matrix(given),1,function(z){
      i <- one2ncol.xx[apply(abs(xx - z),2,sum) < tol]
      if (length(i)==0) return(0)
      if (length(i)==1) return(i)
      return(NA)
    })
    if (any(is.na(index)))
      stop("identification of the given data points is not unique -- `tol' too large?")
    if (any(notfound <- (index==0))) {
      index[notfound] <- (ncol(xx)+1):(ncol(xx)+sum(notfound))
    }

    z <- GaussRF(rbind(t(xx), given[notfound,,drop=FALSE]),
                   model=model,param=param,grid=FALSE,n=n)
    xx <- NULL
  }
  if (is.null(z)) stop("random fields simulation failed")
  
  param[1] <- 0
  if (n==1) {
    zgiven <- z[index] ## z values at the `given' locations
    z <- z[1:total]
  } else {
    ## this is a bit more complicated since index is either a vector or
    ## a matrix of dimension dim(z)-1
    zgiven <- apply(z,length(dim(z)),function(m) m[index])
    z <- as.vector(apply(z,length(dim(z)),function(m) m[1:total]))
  }

  z + Kriging(method=method,
              x=x, grid=grid,
              model=model,param=param,
              given=given,data=as.vector(data)-zgiven,
              gridtriple=TRUE)
}  
  
