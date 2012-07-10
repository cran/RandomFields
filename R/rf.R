# source("rf.R")
# see getNset.R for the variable .methods

# readline <- function(...) return("") 

.onAttach <- function (lib, pkg) {
##  cat("As all other package writers of R, I am happy if you cite RandomFields\n",
##      "whenever you use it. To get the correct reference please type in\n",
##      "citation(\"RandomFields\")    or   \n",
##      "help(RandomFields)   and see the Reference section for a complete list\n",
##      "Many thanks,  Martin Schlather\n")
##  cat("\n\n\n\n\n ACHTUNG!! wird bei DO richtig erkannt, wann res[..]=0 gesetzt werden muss??\n\n")
    packageStartupMessage("##  Note that several updates of RandomFields are expected during 2011. ##\n##  Please see help(\"changings\") for important changes.                 ##\n")
 
  
}

.onUnload <- function(lib, pkg){
   }
#Implementierung von Cox & Isham's non-separable model

#individuelle Praeferenzliste:
#  1) erste praeferenz
#  2) if # pkte < 10 / 50 -> Gauss??
#  4) if nugget/variance > 0.01 -> CE (and not spectral, e.g.)
#  3) if C(maxh)/C(0) > 0.05  -> TBM else CE


.onUnload <- function(lib, pkg){
  DeleteAllRegisters()
}


Covariance <- function(x, y=NULL, model, param=NULL,
                       dim = ifelse(is.matrix(x), ncol(x), 1),
                       Distances,
                       fctcall=c("Cov", "Variogram", "CovMatrix")) {
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic Covariance expects distances as values of x
  ##
  ## here, in contrast to Covariance, nonstatCovMatrix needs only x

  stopifnot(xor(missing(x), missing(Distances)))
  if (are.dist <- missing(x)) x <- Distances

  if (ismatrix.x <- is.matrix(x)) stopifnot(dim == ncol(x))

  xy <- list(x=x, y=y)
  for (i in 1:2) {
    if (!is.null(xy[[i]]) && !is.numeric(xy[[i]])) {
      if (ismatrix.x) stop("values of x are not numeric")
      if (RFparameters()$Print>0) 
        warning(paste(names(xy)[i],
                      "not being numeric it is converted to numeric"))
      assign(names(xy)[i], as.numeric(xy[[i]]))
    }
  }
  remove(xy)


  x <- as.matrix(x)## C code expects the coordinates columwise
  fctcall <- match.arg(fctcall)
  
  stopifnot(length(x)!=0, 
            all(is.finite(x)),
            length(dim) == 1,
            is.finite(dim))
  xdim <- as.integer(ncol(x)) ## may differ from dim if dim>1 and x gives distance,
  ##                 instead of distance vectors
  
  if (!is.null(y)) {
    stopifnot(length(x) == length(y), all(dim == ncol(y)))
    y <- t(y)  ## C code expects the coordinates columwise
    storage.mode(y) <- "double"
    stopifnot(fctcall != "CovMatrix")
 }
  
  storage.mode(dim) <- "integer"
  vdim <- .Call("CheckModelUser", PrepareModel(model, param)$model,
                as.integer(dim), as.integer(xdim), is.null(y),
                PACKAGE="RandomFields")

  len <- nrow(x)
  x <- t(x)
  storage.mode(x) <- "double"
  dim <- as.integer(dim)
 
  
  if (fctcall == "CovMatrix") {    
    if (are.dist) {
      stopifnot(is.null(y))
      len.old <- len
      len <- (1 + sqrt(1 + 8 * len)) / 2 # since the x's give
      ##                       the upper triangle of a quadratic matrix
      stopifnot(len * (len-1) / 2 == len.old)# only stationary case currently
    }
    size <- len * vdim
    result <- double(size^2)     
    .C(fctcall, x, as.integer(are.dist), as.integer(len), result, 
       PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)

    ##  Print(result, size, len, vdim, x, as.integer(are.dist)); xxxx
    
    dim(result) <- rep(size, 2)
  } else {
    if (are.dist)
      stop("Distances are only allowed when calculating covariance matrices")    
    result <- double(len * vdim * vdim)
    .C(fctcall, x, if (is.null(y)) double(0) else y,
       as.integer(is.null(y)),
       as.integer(len), result, 
     PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)
    if (vdim > 1) dim(result) <- c(len * vdim , vdim) ###
  }
      
  return(result)
}


CovarianceFct <- function(...) Covariance(...)

CovMatrix <- function(...) {
  Covariance(..., fctcall="CovMatrix")
}

DeleteAllRegisters <- function () {
    .C("DeleteAllKeys", PACKAGE="RandomFields")
    invisible()
}

DeleteRegister <- function (nr = 0) {
   stopifnot( length(nr) == 1, is.finite(nr) ) 
   .C("DeleteKey", as.integer(nr), PACKAGE="RandomFields")
   invisible()
}


DoSimulateRF <- function (n = 1, register = 0, paired=FALSE, trend=NULL) {
  stopifnot(length(n) == 1,
            n>0, is.finite(n),
            length(register) == 1,
            is.finite(register)
            )
  storage.mode(register) <- "integer"
  DoNameList <- c("DoSimulateRF", "DoSimulateRF", "DoMaxStableRF")  
  assign(".p", GetrfParameters(TRUE))
  keyinfo <- .C("GetKeyInfo", register, 
                total = integer(1), len = integer(.p$maxdim["simu"]),
                spatialdim = integer(1), timespacedim = integer(1),
                grid = integer(1), distr = integer(1),
                maxdim = as.integer(.p$maxdim["simu"]), vdim=integer(1),
                PACKAGE="RandomFields", DUP=FALSE)
  if (paired && (n %% 2 != 0)) stop("if paired, then n must be an even number")
  if (keyinfo$total <= 0)
    stop("register ", register, " does not look initialized")
  keyinfo$len <- keyinfo$len[1:keyinfo$timespacedim]
  error <- integer(1)

  if (!is.null(trend) && (keyinfo$distr<=1)) {
      x <- GetRegisterInfo(register)$loc
  }
  result <- double(keyinfo$total * n * keyinfo$vdim)
  .C(DoNameList[1+keyinfo$distr], register, as.integer(n),
     as.integer(paired), result, error, PACKAGE="RandomFields", DUP=FALSE)
  if (error) stop("error ", error); ## Hallo Marco,
  ## Du hattest hier auskommentiert, weiss nicht wieso
  dimres <- c(keyinfo$vdim,
              if (keyinfo$grid) keyinfo$len else
              if (keyinfo$spatialdim == keyinfo$timespacedim) keyinfo$total
              else keyinfo$len[c(1,length(keyinfo$len))], #
              if (n>1) n)
  dim(result) <- dimres
  if (!is.null(trend) && (keyinfo$distr<=1)) {
   if (keyinfo$vdim == 1) {
     trend<-list(trend)
   } else {
     if ((!is.list(trend) || (length(trend) != keyinfo$vdim)))
       stop(paste("Due to the covariance model a ", keyinfo$vdim, "-dimensional random
        field is created. Then trend should be a list of trends for each component."))
   }
 
   time.used <- keyinfo$timespacedim - keyinfo$spatialdim
   spacedim <- keyinfo$spatialdim
   if (x$grid) {
      zz <- as.matrix(cbind(matrix(x$xgr, nrow=3), x$T))
      xlist <- mapply(seq, zz[1,], zz[2,], zz[3,], SIMPLIFY=FALSE)
      dimension <- sapply(xlist, length)
      ## 'x' will be used in apply within kriging
      x <- as.matrix(do.call(expand.grid, xlist))
   } else if (time.used) {
     x$x <- matrix(x$x, nrow=spacedim)
     x <- as.matrix(cbind(
          matrix(rep(x$x,times=length(seq(x$T[1],x$T[2],x$T[3]))),
	  ncol=nrow(x$x), byrow=TRUE),
          rep(seq(x$T[1],x$T[2],x$T[3]), each=ncol(x$x))))
   } else x <- t(matrix(x$x, nrow=spacedim))

   len <- prod(dimres[-1])
   for (k in 1:keyinfo$vdim) {
     idx <- ((k-1) * len + 1) : (k * len)
    
     if (class(trend[[k]])=="numeric") {
       if (length(trend[[k]])==1) result[idx] <- result[idx] + trend[[k]]
       else if (length(trend[[k]]) != ncol(x) + 1)
          stop("if trend is a vector, its length should be d+1 with d being the
                dimension of the points to be simulated at (space + time)")
        else   {
	   dim(trend[[k]]) <- c(length(trend[[k]]), 1)
           result[idx] <- result[idx] + (cbind(1,x) %*% trend[[k]])[,1]
        }
     } else if (class(trend[[k]])=="function")  {
       variablenames <- variables <- NULL
       if (any(!is.na(match(c("x","y","z"), names(formals(trend[[k]]))))))
            variablenames <- c(c("x","y","z")[1:(ncol(x)-time.used)], if(time.used) "T")
       else variablenames <- c(paste("X",1:(ncol(x)-time.used),sep=""), if(time.used) "T")
       for(i in 1:ncol(x))
          variables <- cbind(variables, paste(variablenames[i],"=x[,",i,"]",sep=""))
       argmatch <- match(variablenames, names(formals(trend[[k]])))
       args <- paste(variables[!is.na(argmatch)], collapse=",")
       result[idx] <- result[idx] + eval(parse(text=paste("trend[[k]](",args,")", sep="")))
     } else if(!is.null(trend[[k]])) stop("trend should be an integer, vector or function")
   }
  }

  if (dimres[1] == 1) dim(result) <- dimres[-1]
  if (length(dim(result))==1) result <- as.vector(result)
  
  return(result)
}


InitSimulateRF <- function (x, y = NULL, z = NULL, T=NULL,
                            grid=!missing(gridtriple), model, param,
                            trend, method = NULL,
                            register = 0, gridtriple, distribution=NA) {

  InitNameList <- c("InitSimulateRF","InitSimulateRF","InitMaxStableRF")
  if (is.na(distribution)) {
    stop("This function is an internal function.\nUse `GaussRF', `InitGaussRF', `MaxStableRF', etc., instead.\n")    
  }
  
  distrNr <- .C("GetDistrNr", as.character(distribution), as.integer(1),
                nr=integer(1), PACKAGE="RandomFields")$nr
  if (distrNr<0)
    stop("Unknown distribution -- do not use this internal function directly")
  else {
    if ((distrNr==.C("GetDistrNr", as.character("Poisson"), as.integer(1),
           nr=integer(1), PACKAGE="RandomFields")$nr) && !exists("PoissonRF")) 
      stop("Sorry. Not programmed yet.")  ##
    InitName <- InitNameList[distrNr + 1]
  }

 
  neu <- CheckXT(x, y, z, T, grid, gridtriple)
  p <- PrepareModel(model, param, trend=trend, method=method)

# Print(p)
#  xxxxx
  
  dim <- as.integer(ncol(neu$x) + neu$Time)
  stopifnot(length(dim) > 0)

  
  .Call("CheckModelSimu", p$model, dim, dim, FALSE, PACKAGE="RandomFields")  

  ## "StoreTrend must be called before .C(InitName, ...)" !!
  ## see extrems.cc for example
 
  error <- .C("StoreTrend", as.integer(register), as.integer(0),
              as.character(""), as.integer(nchar("")),
              as.double(0), as.integer(1),
              err=integer(1), PACKAGE="RandomFields")$err
  if (error) stop("trend not correct -- error nr", error)

  
  .C(InitName, as.double(neu$x), as.double(neu$T),
              as.integer(neu$spacedim),
              as.integer(neu$l), 
              as.integer(neu$grid),
              as.integer(neu$Time),
              as.integer(distrNr),
              as.integer(register),
              as.integer(1),
              error=integer(1),
              PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)

   return(error)
}


InitGaussRF <- function(x, y = NULL, z = NULL, T=NULL, grid=!missing(gridtriple),
                        model, param,
                        trend=NULL, method = NULL,
                        register = 0, gridtriple) {
  return(InitSimulateRF(x=x, y=y, z=z,  T=T, 
                        grid=grid, model=model, param=param, trend=trend,
                        method=method, register=register,
                        gridtriple=gridtriple, distribution="Gauss"))
}


GaussRF <- function (x, y = NULL, z = NULL, T=NULL,
          grid=!missing(gridtriple), model, param, trend=NULL, method = NULL, 
          n = 1, register = 0, gridtriple,
          paired=FALSE, ...) {
  RFpar <- list(...)
  old.param <- NULL
  delete <- n>1 && !(old.param <- RFparameters(no.readonly=TRUE))$Storing
  modifiedParam <- length(RFpar)>0 || delete
  if (modifiedParam) {
    if (is.null(old.param))
      old.param <- RFparameters(no.readonly=TRUE)
    if (length(RFpar)>0) RFparameters(RFpar)
    if (delete) RFparameters(Storing=TRUE)
  }
  on.exit(if (modifiedParam) {
    RFparameters(old.param);
    if (delete) DeleteRegister(register);
  })

  p <- PrepareModel(model, param, trend=trend, method=method)

  error <- InitSimulateRF(x=x, y=y, z=z, T=T, grid=grid, model=model,
                          param=param,
                          trend=trend, method=method, register=register,
                          gridtriple=gridtriple, distribution="Gauss")

  if (error > 0)
    stop("Simulation could not be initiated.",
         if (RFparameters()$PrintLevel < 3) "\nRerun with higher value of RFparameters()$PrintLevel for more information. (Or put debug=TRUE if you are using Showmodels.)\n\n")
  res <- (if (n<1) NULL
           else DoSimulateRF(n=n, register=register, paired=paired,
                             trend=p$trend))
  return(res)
}


Variogram <- function (x, model, param=NULL, dim=ifelse(is.matrix(x),ncol(x),1))
  CovarianceFct(x=x, model=model, param=param, dim=dim, fctcall="Variogram")

