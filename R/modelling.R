
## to.do replace pure nugget effect model as error.model by arbitrarz model

Kriging <- function(krige.method, x, y=NULL, z=NULL, grid, gridtriple=FALSE,
                    model, param, given, data, pch=".") {
  krige.methodlist <- c("S","O")
  if (is.na(krige.method.nr <- pmatch(krige.method,krige.methodlist)))
    stop(paste("kriging method not identifiable from the list.",
               "Possible values are",paste(krige.methodlist,collapse=",")))
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
  if (ncol(x)!=nrow(tgiven))
    stop("dimension of the kriged points and the given points do not match")
  given <- NULL

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

  nd.step <-  ceiling(nrow(x) / 79)
  env <- environment()
  assign("enu",0,envir=env)
  
  covnr <- as.integer(.C("GetModelNr", as.character(model), nr = integer(1))$nr)
  pp <- param
  if (RFparameters()$PracticalRange) {
    pp[ip$SCALE] <- pp[ip$SCALE] / GetPracticalRange(model,pp[-1:(-ip$KAPPA+1)])
  }
  
  storage.mode(pp) <- "double"
  nn <- as.integer(ncol(tgiven))
  np <- as.integer(length(pp))

  switch(krige.method.nr,
         { ## simple kriging
           data <- as.matrix(solve(covmatrix,data-param[ip$MEAN]))
           res <- apply(x, 1,
                        function(z){
                          if ((enu %% nd.step) == 0) cat(pch);
                          assign("enu",enu+1,envir=env);
                          ## to.do replace the following by unchecked !
                          ## CovFct has been check already !
                          .C("UncheckedCovFct",
                             as.double(.C("onePdist",as.double(tgiven),
                                          as.integer(nrow(tgiven)),
                                          as.integer(ncol(tgiven)),
                                          as.double(z),
                                          dx=double(ncol(tgiven)),DUP=FALSE)$dx),
                             nn,covnr,pp,np,cov=double(nn),DUP=FALSE)$cov %*%
                               data}
                        ) + pp[ip$MEAN]
         }, {
           ## ordinary kriging
           covmatrix <- rbind(cbind(covmatrix,1), c(rep(1,nd),0)) 
           data <- solve(covmatrix,rbind(as.matrix(data),0))
           res <-
             apply(x, 1,
                   function(z){
                     if ((enu %% nd.step) == 0) cat(pch);
                     assign("enu",enu+1,envir=env);
                     c(.C("UncheckedCovFct",
                          as.double(.C("onePdist",as.double(tgiven),
                                       as.integer(nrow(tgiven)),
                                       as.integer(ncol(tgiven)),
                                       as.double(z),
                                       dx=double(ncol(tgiven)),DP=FALSE)$dx),
                          nn,covnr,pp,np,cov=double(nn),DUP=FALSE)$cov,
                       1) %*% data 
                   })
         }
         ) # switch

  if (pch!="") cat("\n")
  x <- data <- NULL
  if (is.matrix.data) {
    if (grid) return(array(t(res),dim=c(dimension,ncol.data)))
    else return(t(res))
  } else {
    if (grid) return(array(res,dim=dimension))
    else return(res)
  }
}



CondSimu <- function(krige.method, x, y=NULL, z=NULL, grid, model, param, 
                     method=NULL, n=1, register=0, gridtriple=FALSE,
                     err.model=NULL, err.param=NULL, err.method=NULL,
                     err.register=1,
                     given, data, tol=1E-5, pch=".") {
  if (is.character(method) && (is.na(pmatch(method, c("S","O")))))
    stop("Sorry. The parameters of the function `CondSimu' as been redefined. Use `krige.method' instead of `method'. See help(CondSimu) for details")

  ip<-.C("GetParameterIndices", MEAN=integer(1), VARIANCE=integer(1),
          NUGGET=integer(1), SCALE=integer(1), 
          KAPPA=integer(1), LASTKAPPA=integer(1),integer(1),
          SILL=integer(1), DUP=FALSE)
  ip <- lapply(ip,function(i){i+1})
  
  if (is.null(err.model)) {
    krige.model <- model
    krige.param <- param
  } else {
    if (.C("GetModelNr", as.character(err.model), nr=integer(1))$nr !=
        .C("GetModelNr", as.character("nugget"), nr=integer(1))$nr)
      stop("currently only pure nugget effect allowed as error.model")
    krige.model <- model  ## this has to be replaced by an additive model ->to.do
    krige.param <- param  ## dito
    krige.param[ip$NUGGET] <- krige.param[ip$NUGGET] +
      err.param[ip$VARIANCE] + err.param[ip$NUGGET]
  }
  
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
    outside.grid <- apply((abs(ind-index)>tol) | (index<1) |
                          (index>dimension),2,any)   
    if (any(outside.grid))
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
      eval(parse(text=paste("ll <- c(",
                   paste("length(seq(x[1,",1:l,"],x[2,",1:l,"],x[3,",1:l,"]))",
                         collapse=","),
                   ")")))
     # print(ll)
     # print(outside.grid)
     # print(index[,!outside.grid,drop=FALSE])
      new.index <- rep(0,ncol(index))
      if (!all(outside.grid)) {
        new.index[!outside.grid] <- 1 +
          apply((index[,!outside.grid,drop=FALSE]-1) *
                cumprod(c(1,ll[-length(ll)])), 2, sum)
      }
                                        # print(new.index)
      index <- new.index
      new.index <- NULL
    } else {
      ## data points are all lying on the grid
      z <- GaussRF(x, grid=TRUE, model=model, param=param,
                   method=method, n=n, register=register,
                   gridtriple=TRUE)
      index <- t(index)
    }
  } else {
    xx <- t(as.matrix(x)) ## not a grid
    ## the next step can be pretty time consuming!!!
    ## to.do: programme it in C
    one2ncol.xx <- 1:ncol(xx)
    index <- apply(as.matrix(given),1,function(z){
      i <- one2ncol.xx[apply(abs(xx - z),2,sum) < tol]
      if (length(i)==0) return(0)
      if (length(i)==1) return(i)
      return(NA)
    })
  }

  if (!simu.grid) {
    tol <- tol * nrow(xx)
    if (any(is.na(index)))
      stop("identification of the given data points is not unique -- `tol' too large?")
    if (any(notfound <- (index==0))) {
      index[notfound] <- (ncol(xx)+1):(ncol(xx)+sum(notfound))
    }

    xx <- rbind(t(xx), given[notfound,,drop=FALSE])
    z <- GaussRF(xx, grid=FALSE, model=model, param=param,
                 method=method, n=n, register=register)

    xx <- NULL
  }
  if (is.null(z)) stop("random fields simulation failed")

  if (n==1) {
    ## z values at the `given' locations
    zgiven <- z[index]
    z <- z[1:total]
  } else {
    ## this is a bit more complicated since index is either a vector or
    ## a matrix of dimension dim(z)-1
    zgiven <- apply(z,length(dim(z)),function(m) m[index])
    z <- as.vector(apply(z,length(dim(z)),function(m) m[1:total]))
  }

  if (!is.null(err.model)) {
     error <- GaussRF(given, grid=FALSE, model=err.model, param=err.param,
                     method=err.method, n=n, register=err.register)  
     if (is.null(error)) stop("error field simulation failed")
     zgiven <- zgiven + as.vector(error)
     error <- NULL
   }
 
  param[1] <- 0
  z + Kriging(krige.method=krige.method,
              x=x, grid=grid,
              model=krige.model,param=krige.param,
              given=given,data=as.vector(data)-zgiven,
              gridtriple=TRUE, pch=pch)
}  
  
