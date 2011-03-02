Kriging <- function(krige.method, x, y=NULL, z=NULL, T=NULL,
                    grid, gridtriple=FALSE,
                    model, param, given, data, trend=NULL,            
                    pch=".", return.variance=FALSE,
                    allowdistanceZero = FALSE, cholesky=FALSE) {
  ## currently only univariate Kriging possible
 
  ##  MINV <<- NULL
  SOLVE <- if (cholesky) {
    function(M, v) {
      ##    Print("cholesky!")
      if (missing(v)) chol2inv(chol(M)) else chol2inv(chol(M)) %*% v
    }
  } else solve
  
  silent = TRUE
  nuggetrange <- 1e-15 + 2 * RFparameters()$nugget.tol
  krige.methodlist <- c("S", "O", "U", "I")
  krige.method.nr <- pmatch(krige.method, krige.methodlist)
 
  given <- as.matrix(given)
  ngiven <- as.integer(nrow(given)) ## number of given points
  xdim <- as.integer(ncol(given))
  
  if (is.na(krige.method.nr))
    stop(paste("kriging method not identifiable from the list.",
               "Possible values are", paste(krige.methodlist,
                                            collapse = ",")))

  x <- CheckXT(x = x, y = y, z = z, T = T, grid = grid,
               gridtriple = gridtriple)
  y <- z <- NULL
  time.used <- (!is.null(x$T))
  if (x$spacedim + x$Time != xdim)
    stop("dimensions of the kriged points and the given points do not match")
  pm <- PrepareModel(model=model, param=param, trend=trend)
  vdim <- .Call("CheckModelIntern", pm$model, xdim, xdim,
                TRUE,  ## stationary !
                PACKAGE="RandomFields")

  data <- as.double(data)
  repet <- length(data) / (vdim *  ngiven)
  if (repet != as.integer(repet)) {
    Print(length(data), vdim, ngiven)
    stop("number of data not a multiple of the number of locations")
  }
  storage.mode(repet) <- "integer"
 
  dim(data) <- c(ngiven, vdim, repet)
  tgiven <- t(given)
  storage.mode(tgiven) <- "double"
 

  pos <- integer(ngiven)
  ## lexicographical ordering of vectors --- necessary to check
  ## whether any location is given twice, but with different value of
  ## the data
  .C("Ordering", tgiven, as.integer(ngiven), xdim,
     pos, PACKAGE = "RandomFields", DUP = FALSE)
  pos <- pos + 1

  
  ## are locations given twice with different data values?
  check <- tgiven[ , pos, drop = FALSE]
  if (any(dup <- c(FALSE, colSums(abs(check[, -1, drop = FALSE] -
                                      check[, -ngiven, drop = FALSE])) == 0))) {

  #  Print("dup")
    
    if (allowdistanceZero) {
      tgiven[1, dup] <- check[1, dup] + rnorm(n=sum(dup), 0, nuggetrange)
      data <- data[pos, , , drop = FALSE]
    } else {
      data <- data[pos, , , drop = FALSE] # pos is recysled if vdim > 1
      if (any(data[dup, , ] != data[c(dup[-1], FALSE), , ]))
        stop("duplicated conditioning points with different data values")
      tgiven <- tgiven[, !dup, drop = FALSE]
      data <- data[!dup, , drop = FALSE]
      ngiven <- as.integer(ncol(given)) ## may some point had been deleted
      stop("duplicated values not allowed")
    }
  }

  notna <- as.vector(apply(is.finite(data), 1:2, all))
#  Print(notna)
  data <- data[notna]
  dim(data) <- c(length(data) / vdim / repet, vdim, repet)

  
 # Print(length(data), dim(data))
  
  x$x <- as.matrix(x$x)
  if (grid) {
    zz <- as.matrix(cbind(x$x, x$T))
    xlist <- mapply(seq, zz[1,], zz[2,], zz[3,], SIMPLIFY=FALSE)
    dimension <- sapply(xlist, length)
    ## 'x' will be used in apply within kriging
    x <- as.matrix(do.call(expand.grid, xlist))
  }
  else if (time.used)
    x <- cbind(matrix(rep(t(x$x),times=length(seq(x$T[1],x$T[2],x$T[3]))),
                      ncol=ncol(x$x), byrow=TRUE),
               rep(seq(x$T[1],x$T[2],x$T[3]), each=nrow(x$x)))
  else x <- x$x

  nx <- nrow(x)
  storage.mode(x) = "double"
  dimgiven <- as.integer(dim(tgiven))
  ngvdim <- ngiven * vdim
  covmatrix <- double(ngvdim^2)

#  Print( ngiven , vdim, dimgiven, ngvdim)

#  Print(matrix(vd, nr=3))

  
  ## works also for variograms and returns -(gamma(x_i-x_j))_{ij}
  .C("CovMatrixIntern", tgiven, as.integer(FALSE), ngiven, covmatrix,
       PACKAGE="RandomFields", DUP = FALSE)
  dim(covmatrix) <- c(ngvdim, ngvdim)
  covmatrix <- covmatrix[notna, notna]


 #  Print(dim(covmatrix), covmatrix[1:10, 1:10], notna, ngvdim, ngiven, vdim)
 
 
  if(class(trend)=="formula")  {
    fo <- as.list(trend)[[length(as.list(trend))]]
    functions <- dummy <- functionargs <- NULL
    argnames <-c("x","y","z",paste("X",1:(xdim-time.used),sep=""),"T")
    searchargs <- lapply(argnames, function(x) grep(x, as.character(fo)))
    argsfound <- sapply(searchargs, function(x) length(x)>0)
    functionargs <- argnames[argsfound]
    functionargs <- paste(functionargs, collapse=",")
    while (as.list(fo)[[1]]=="+") {
      dummy <- as.list(fo)
      fo <- dummy[[2]]
      functions <- c(eval(parse(text=paste("2function(",functionargs,") ",
                                  as.character(dummy)[3],sep=""))), functions)
    }
    if (is.null(dummy))
      functions <-
        as.list(c(eval(parse(text=paste("function(", functionargs,") ",
                               as.character(as.list(trend))[3],sep=""))),
                  functions))
    else
      functions <-
        as.list(c(eval(parse(text=paste("function(",functionargs,") ",
                               as.character(dummy)[2],sep=""))), functions))
    trend <- functions
  }
  if(class(trend)=="list")  {
    nfct <- length(trend)
    functiontexts <- as.list(rep("", times = nfct))
    for (j in 1:nfct) {
      variablenames <- variables <- NULL
      if (any(!is.na(match(c("x","y","z"), names(formals(trend[[j]]))))))
        variablenames <-c(c("x","y","z")[1:(xdim - time.used)], if(time.used) "T")
      else variablenames <- c(paste("X",1:(xdim - time.used),sep=""),
                              if(time.used) "T")
      for(i in 1:ncol(x)) {
        variables <-
          cbind(variables, paste(variablenames[i],"=points[,",i,"]",sep=""))
      }
      argmatch <- match(variablenames, names(formals(trend[[j]])))
      args <- paste(variables[!is.na(argmatch)], collapse=",")
      functiontexts[[j]] <- paste("trend[[",j,"]](", args,")", sep="")
      ## functiontexts is a list with entries like
      ## "trend[[1]](x=points[,1],y=points[,2])"
    }
  }
  if(class(trend) == "numeric")  {
    if (as.integer(trend) != trend || trend < 0)
      stop("Trend should be a nonnegative integer.")
    Powers <- as.matrix(do.call("expand.grid", rep(list(0:trend), times=xdim)))
    Powers <- matrix(Powers[rowSums(Powers)<=trend, ], ncol = xdim)
    nfct <- nrow(Powers)
  }
  
  evaltrendfunction <- function(text, trend, functionnumber, points)  {
    points <- as.matrix(points)
    eval(parse(text=text[[functionnumber]]))
  }
  
  res <- double(nx * repet * vdim)
  switch(krige.method.nr, {
    ## simple kriging
    stopifnot(is.null(pm$trend))
    if (is.numeric(trend) && length(trend) == vdim) {
      mean <- trend
    } else {
      stopifnot(is.null(trend))
      mean <- rep(0, vdim)
    }
    mean <- as.double(rep(mean, each = ngiven))
    
    if (return.variance) {
      if (vdim > 1) stop("multivariate version not programmed yet.")
      if (!(is.numeric(try(invcov <- SOLVE(covmatrix), silent = silent)))) {
        stop("Covmatrix is singular")
      }
      
      sigma2 <- double(nx)
      .C("simpleKriging2", tgiven, x, as.double(data - mean), as.double(invcov),
         nx, ngiven, xdim, repet, mean, res, sigma2,
         PACKAGE = "RandomFields", DUP = FALSE)
    } else {
      dim(data) <- c(ngiven * vdim, repet)
      
      if (!(is.numeric(try(invcov <- SOLVE(covmatrix, data - mean),
                             silent = silent)))) {
          Print(covmatrix, pm ,mean, sum(data), dim(data), mean,
                dim(covmatrix), eigen(covmatrix)$values)
          stop("Covmatrix is singular.")
        
      } 
      
   #   Print(data, covmatrix, covmatrix, invcov, notna)
      
      if (!all(notna)) {        
        dummy <- invcov
        invcov <- numeric(ngiven)
        invcov[notna] <- dummy
        rm(dummy)
      }
      
##      MINV <<- invcov
      
      .C("simpleKriging", tgiven, x, as.double(invcov), as.integer(notna),
         nx, ngiven, xdim, repet, mean, res,
         PACKAGE = "RandomFields", DUP = FALSE)
    }
  }, {
    ## ordinary kriging
    stopifnot(is.null(pm$trend) || is.null(trend))
    covmatrix <- rbind(cbind(covmatrix, 1), c(rep(1, ngiven), 0))
    if (vdim > 1) stop("multivariate version not programmed yet.")
    dim(data) <- c(ngiven * vdim, repet)
    if (return.variance) {
      if (!(is.numeric(try(invcov <- SOLVE(covmatrix), silent = silent)))) {
                                        #       Print(covmatrix)
        stop("Covmatrix is singular!")
      }
      sigma2 <- double(nx)
      .C("ordinaryKriging2", tgiven, x, as.double(data), as.double(invcov),
         nx, ngiven, xdim, repet, res, sigma2,
         PACKAGE = "RandomFields", DUP = FALSE)
    }  else  {
#      Print(dim(covmatrix), dim(data), dim(as.matrix(data)),
#              dim( rbind(as.matrix(data), 0) ))
     if (!(is.numeric(try(invcov <- SOLVE(covmatrix, rbind(data, 0)),
                           silent = silent)))) {
        ##       Print(covmatrix)
         stop("Covmatrix is singular .") #
      }
      if (FALSE)
        if (RFparameters()$Print > 3) {
          Print(list("ordinaryKriging", tgiven,
                     x, as.double(invcov), nx,
                     ngiven, xdim, repet,
                     res, PACKAGE = "RandomFields", DUP = FALSE))
          Print("eigen(cov)", covmatrix, covmatrix, eigen(covmatrix),
                RFparameters()$Practical)
        }
      .C("ordinaryKriging", tgiven, x,  as.double(invcov), nx, ngiven, xdim,
         repet, res,
         PACKAGE = "RandomFields", DUP = FALSE)
    }
  }, {
    ## universal kriging
    if (vdim > 1) stop("multivariate version not programmed yet.")
    if(is.na(pmatch(class(trend), c("numeric", "list"))))
      stop(paste("Trend method not identifiable from the list.",
                 "Possible values are numeric, list, formula."))
    if (return.variance) {
      if (!(is.numeric(try(invcov <- SOLVE(covmatrix), silent = silent))))
        stop("covmatrix is singular!")
      sigma2 <- double(nx)
      if (class(trend) == "numeric")   {
        .C("universalKriging4", tgiven, x,
           as.double(data), as.double(invcov), nx, ngiven, xdim, repet, res,
           sigma2,  as.double(t(Powers)), as.integer(nfct),
           PACKAGE = "RandomFields", DUP = FALSE)
      } else {
        Fmatrix <- matrix(double(ngiven*nfct), nrow = ngiven)
        for (j in 1:nfct)
          Fmatrix[,j] <- evaltrendfunction(functiontexts, trend, j, given)
        f_expr <- quote(evaltrendfunction(functiontexts, trend, j,
                                          x[i,,drop=FALSE]))
        i <- integer(1) ##in- and output variable
        j <- integer(1) ##in- and output variable
        .C("universalKriging2", tgiven, x,
           as.double(data), as.double(invcov), nx, ngiven, xdim, repet, res,
           sigma2, as.double(t(Fmatrix)), as.integer(nfct), i, j, f_expr,
           new.env(), PACKAGE = "RandomFields", DUP = FALSE)
      }
    }
    else  {
     dim(data) <- c(ngiven * vdim, repet)
     data_null <- rbind(data, matrix(0.0, nrow=nfct, ncol=repet))
      Fmatrix <- matrix(ncol=nfct, nrow=ngiven)
      if (class(trend) == "numeric")   {
        for(j in 1:nfct)
          Fmatrix[,j] <- apply(as.matrix(tgiven^Powers[j,]), 2, prod)
        KMat <- rbind(cbind(covmatrix,Fmatrix), 
                      cbind(t(Fmatrix), matrix(0.0 ,nrow=nfct, ncol=nfct)))
        if (!(is.numeric(try(invcov <- SOLVE(KMat, data_null, silent=silent)))))
          stop("covmatrix is singular..")
        .C("universalKriging3", tgiven, x,
           as.double(invcov), nx, ngiven, xdim,
           repet, res, as.double(t(Powers)), 
           as.integer(nfct), PACKAGE = "RandomFields", DUP = FALSE)
      }
      else  {
        Fmatrix <- matrix(double(ngiven*nfct), nrow = ngiven)
        for (j in 1:nfct)
          Fmatrix[,j] <- evaltrendfunction(functiontexts, trend, j, given)
        KMat <- rbind(cbind(covmatrix,Fmatrix),
                      cbind(t(Fmatrix), matrix(0.0, nrow=nfct, ncol=nfct))) 
        if (!(is.numeric(try(invcov <- SOLVE(KMat, data_null, silent=silent)))))
          stop("covmatrix is singular...")
        f_expr <- quote(evaltrendfunction(functiontexts, trend, j,
                                          x[i,,drop=FALSE]))
        i <- integer(1)	##in- and output variable
        j <- integer(1) ##in- and output variable
        .C("universalKriging", tgiven, x, as.double(invcov), nx, ngiven, xdim,
           repet, res, as.double(t(Fmatrix)), as.integer(nfct), i, j, f_expr,
           new.env(), PACKAGE = "RandomFields", DUP = FALSE)
      }
    }
  }, {
    ## intrinsic kriging, notation see p. 169, Chiles
   
    if (vdim > 1) stop("multivariate version not programmed yet.")
    if (!(is.numeric(trend)))
      stop("trend should be a number/integer (order k of the intrinsic IRF(k) model)")
    if (return.variance)  {
      ## Inverting CovMatrix and initializing sigma 2
      if (!(is.matrix(try(invcov <- SOLVE(covmatrix), silent = silent))))
        stop("covmatrix is singular .")
      sigma2 <- double(nx)
      inv_expr <- quote(SOLVE(Q))
      ##now: Kriging
      .C("intrinsicKriging2", tgiven, x, as.double(data), as.double(invcov),
         nx, ngiven, xdim, repet, res, sigma2,
         as.double(t(Powers)), as.integer(nfct), inv_expr, new.env(),
         PACKAGE = "RandomFields", DUP = FALSE)
    }
    else  {
      ##Initializing Fmatrix
      Fmatrix <- matrix(double(ngiven*nfct), nrow=ngiven)
      for(j in 1:nfct)
        Fmatrix[ ,j] <- apply(tgiven^Powers[j,], 2, prod)
      KMat <- rbind(cbind(covmatrix,Fmatrix),
                    cbind(t(Fmatrix), matrix(0.0, ncol=nfct, nrow=nfct)))
      dim(data) <- c(ngiven * vdim, repet)
      data_null <- rbind(data, matrix(0.0, ncol=repet, nrow=nfct))
      if (!(is.numeric(try(invcov <- SOLVE(KMat, data_null, silent=silent)))))
        stop("covmatrix is singular !")
      .C("intrinsicKriging", tgiven, x, as.double(invcov), nx, ngiven, xdim,
         repet, res, as.double(t(Powers)), as.integer(nfct),
         PACKAGE = "RandomFields", DUP = FALSE)
    }
  })

  if (grid) {
    dim(res) <- c(dimension, if (vdim>1) vdim, if (repet>1) repet)
  } else {
    if (vdim > 1 || repet > 1)
      dim(res) <- c(nx, if (vdim>1) vdim, if (repet>1) repet)   
  }
  
  if (pch != "") cat("\n")

  return(if (return.variance)
         list(estim = res,
              var = if (grid) array(sigma2, dim = dimension) else sigma2
              )
         else res)
}


CondSimu <- function(krige.method, x, y=NULL, z=NULL, T=NULL,
                     grid, gridtriple=FALSE,
                     model, param, method=NULL,
                     given, data, trend=NULL,
                     n=1, register=0, 
                     err.model=NULL, err.param=NULL, err.method=NULL,
                     err.register=1, 
                     tol=1E-5, pch=".",
                     paired=FALSE,
                     na.rm=FALSE
                     ) {
  op.list <- c("+","*")  
  stopifnot(is.character(krige.method))
  if (is.character(method) && (!is.na(pmatch(method, c("S","O")))))
    stop("Sorry. The parameters of the function `CondSimu' as been redefined. Use `krige.method' instead of `method'. See help(CondSimu) for details")

  x  <- CheckXT(x, y, z, T, grid, gridtriple)
  total <- x$total
  y <- z <- NULL   
  pm <- PrepareModel(model, param, trend, method=method)
  givendim <- ncol(as.matrix(given))
 # Print(givendim, given, x$spacedim, x$Time) 
  if (givendim != x$spacedim + x$Time)
      stop("dimension of coordinates does not match dimension of given data")
 
  ##  krige.mean <- 0 ## pm$mean + err$mean ??
  krige <- pm$model
  krige.trend <- pm$trend
  if (!is.null(krige.trend) && length(krige.trend)!=1 &&!is.numeric(krige.trend)){
    stop("not programmed yet")
  }

  vdim <- .Call("CheckModelIntern", krige, givendim, givendim, TRUE,
                PACKAGE="RandomFields")
  if (vdim > 1) stop("multivariate version not programmed yet")
  

  if (!is.null(err.model)) {
    err <- PrepareModel(err.model, err.param, method=err.method)$model

    if (xor(is.null(pm$method), is.null(err$method)))
      stop("the simulation method must be defined for the data model and the error model or for none of them")

    if (!is.null(err$trend)) stop("trend for the error term not allowed yet")

    krige <-
      if (err[[1]]=="+") {
        if (krige[[1]]=="+") c(krige, err[-1]) else c(err, list(krige))
      } else {
        if (krige[[1]]=="+") c(krige, list(err)) else list("+", err, krige)
      }

 #   Print("CheckModelIntern", krige, givendim, givendim, TRUE,
#                  PACKAGE="RandomFields", err, err.param)
    
    ## just check whether model common model is ok
    vdim <- .Call("CheckModelIntern", krige, givendim, givendim, TRUE,
                  PACKAGE="RandomFields")
    
  }

 # Print(krige)

  simu.grid <- grid
  given <- as.matrix(given)
  if (nrow(given)!=length(data)) {
    cat("dimension of 'given' does not match 'data'")
    return(NA)
  }
  if (na.rm && any(data.idx <- is.na(data))) {
    data <- data[data.idx]
    given <- given[data.idx, , drop=FALSE]
  }
  
  txt <- "kriging in space time dimensions>3 where not all the point ly on a grid is not possible yet"
  ## if 4 dimensional then the last coordinates should ly on a grid

  ## now check whether and if so, which of the given points belong to the
  ## points where conditional simulation takes place
  if (grid) {
    zz <-  cbind(x$x, x$T)
    ind <- 1 + (t(given) - zz[1, ]) / zz[3, ]
    index <-  round(ind)
    endpoints <- 1 + floor((zz[2, ] - zz[1, ]) / zz[3, ])
    outside.grid <- apply((abs(ind-index)>tol) | (index<1) |
                          (index>endpoints), 2, any)   
    if (any(outside.grid)) {
      ## at least some data points are not on the grid
      ## simulate as if there is no grid
      simu.grid <- FALSE
      ll <- NULL ##  otherwise check will give a warning
      l <- ncol(zz)
      if (l>3) stop(txt)
      if (l==1) xx <- matrix(seq(x$x[1], x$x[2], x$x[3]), nrow=1)
      else  eval(parse(text=paste("xx <-  t(as.matrix(expand.grid(",
                         paste("seq(zz[1,", 1:l,
                               "],zz[2,", 1:l,
                               "],zz[3,", 1:l, "])", collapse=","),
                         ")))")))
      eval(parse(text=paste("ll <- c(",
                   paste("length(seq(zz[1,", 1:l, "],zz[2,", 1:l, "],zz[3,", 1:l,
                         "]))",
                         collapse=","),
                   ")")))

      new.index <- rep(0,ncol(index))
      ## data points that are on the grid, must be registered,
      ## so that they can be used as conditioning points of the grid
      if (!all(outside.grid)) {
        new.index[!outside.grid] <- 1 +
          colSums((index[, !outside.grid, drop=FALSE]-1) *
                  cumprod(c(1, ll[-length(ll)])))
      }
      index <- new.index
      new.index <- NULL
    } else {
      ## data points are all lying on the grid     
      z <- GaussRF(x=x$x, T=x$T,  grid=TRUE, model=model, param=param,
                   method=method, n=n, register=register,
                   gridtriple=TRUE, paired=paired)
      ## for all the other cases of simulation see, below
      index <- t(index)
    }
  } else {
    xx <- t(as.matrix(x$x)) ## not a grid
    
    if (!is.null(x$T)) {
      if (ncol(xx) > 2) stop(txt)
      T <- seq(x$T[1], x$T[2], x$T[3])
      ## multiple the xx structure by length of T,
      ## then add first component of T to first part ... last component of T
      ## to last part
      xx <- rbind(matrix(xx, nrow=nrow(xx), ncol=ncol(xx) * length(T)),
                  as.vector(t(matrix(T, nrow=length(T),ncol(xx)))))
    }
    
    ## the next step can be pretty time consuming!!!
    ## to.do: programme it in C
    ##
    ## identification of the points that are given twice, as points to
    ## be simulated and as data points (up to a tolerance distance !)
    ## this is important in case of nugget effect, since otherwise
    ## the user will be surprised not to get the value of the data at
    ## that point
    one2ncol.xx <- 1:ncol(xx)
    index <- apply(as.matrix(given), 1, function(z){
      i <- one2ncol.xx[colSums(abs(xx - z)) < tol]
      if (length(i)==0) return(0)
      if (length(i)==1) return(i)
      return(NA)
    })
  }  
  
  if (!simu.grid) {
    ## otherwise the simulation has already been performed (see above)
    tol <- tol * nrow(xx)
    if (any(is.na(index)))
      stop("identification of the given data points is not unique - `tol' too large?")
    if (any(notfound <- (index==0))) {
      index[notfound] <- (ncol(xx) + 1) : (ncol(xx) + sum(notfound))
    }
    xx <- rbind(t(xx), given[notfound, , drop=FALSE])

    ## hier fehler
    z <- GaussRF(x=xx, grid=FALSE, model=model, param=param,
                 method=method, n=n, register=register, paired=paired)
    xx <- NULL
  }
  if (is.null(z)) stop("random field simulation failed")
  
  if (n==1) {
    ## z values at the `given' locations
    zgiven <- z[index]
    z <- as.vector(z[1:total]) # as.vector is necessary !! Otherwise
    ##                          is not recognized as a vector
    
  } else {
    ## this is a bit more complicated since index is either a vector or
    ## a matrix of dimension dim(z)-1
    zgiven <- matrix(apply(z, length(dim(z)), function(m) m[index]), ncol=n)
    z <- as.vector(apply(z, length(dim(z)), function(m) m[1:total]))
  }

  if (!is.null(err.model)) {
     error <- GaussRF(given, grid=FALSE, model=err.model, param=err.param,
                      method=err.method, n=n, register=err.register,
                      paired=paired)
     if (is.null(error)) stop("error field simulation failed")
     zgiven <- zgiven + as.vector(error)
     error <- NULL
   }

  if (FALSE)  Print(krige.method=krige.method,
              x=x$x, grid=grid,
              model=krige,
              given=given,data=as.vector(data)-zgiven,
              gridtriple=TRUE, pch=pch)

  # zgiven is matrix
  if (FALSE)
    Print(z, krige.method=krige.method,
          x=x$x, grid=grid,
          model=krige,
          given=given,data=as.vector(data)-zgiven,
          gridtriple=TRUE, pch=pch,
          z, is.vector(z), length(z),                    
          Kriging(krige.method=krige.method,
              x=x$x, grid=grid,
              model=krige,
              given=given,data=as.vector(data)-zgiven,
              gridtriple=TRUE, pch=pch))


 # Print(krige, x$x, grid, given, as.vector(data)-zgiven);
  
  z +  Kriging(krige.method=krige.method,
              x=x$x, grid=grid,
              model=krige,
              given=given,data=as.vector(data)-zgiven,
              gridtriple=TRUE, pch=pch)
}
  
