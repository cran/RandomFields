## source("rf.R")

"CovarianceFct" <-
  function (x, model, param, dim = ifelse(is.matrix(x), ncol(x), 1),
            fctcall="Covariance") 
{
  ## note: * if x is a matrix, the points are expected to be given row-wise
  ##       * if not anisotropic CovarianceFct expects distances as values of x
  ##
  if (ismatrix.x <- is.matrix(x)) stopifnot(dim == ncol(x))
  x <- as.matrix(x)
  stopifnot(length(x)!=0,
            all(is.finite(x)),
            length(dim) == 1,
            is.finite(dim),
            !is.null(pmatch(fctcall, c("Covariance", "Variogram",
                                       "CovarianceMatrix")))
            )
  xdim <- ncol(x) ## may differ from dim if dim>1 and x gives distance,
  ##                 instead of distance vectors
  p <- PrepareModel(model, param, xdim)
  if (!ismatrix.x && p$anisotropy)
    stop("x must be a matrix for anisotropic models")
  
  storage.mode(x) <- "double"
  x <- t(x)  ## C code extect the coordinates columwise
  if (fctcall %in% c("CovarianceMatrix")) {
    len <- (1 + sqrt(1 + 8 * ncol(x))) / 2 # since the x's give
    ##                       the upper triangle of a quadratic matrix     
    reslen <- len^2
  } else {
    reslen <- len <- ncol(x)
  }
  
  result <- .C(fctcall, as.double(x),
               as.integer(len),
               as.integer(p$covnr),
               as.double(p$param), 
               as.integer(length(p$param)), ## control
               as.integer(dim), as.integer(xdim),
               as.integer(length(p$covnr)),
               as.integer(p$anisotropy),
               as.integer(p$op),            
               res=double(reslen),
               PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)$res
  if (fctcall %in% c("CovarianceMatrix")) return(matrix(result, ncol=len))
  else return(result)
}


"DeleteAllRegisters" <-
function () 
{
    .C("DeleteAllKeys", PACKAGE="RandomFields")
    invisible()
}


"DeleteRegister" <-
function (nr = 0) 
{
   stopifnot( length(nr) == 1, is.finite(nr) ) 
   .C("DeleteKey", as.integer(nr), PACKAGE="RandomFields")
   invisible()
}


"DoSimulateRF" <-
function (n = 1, register = 0) 
{
    
  stopifnot(length(n) == 1,
            n>0, is.finite(n),
            length(register) == 1,
            is.finite(register)
            )
  DoNameList <- c("DoSimulateRF", "DoSimulateRF", "DoMaxStableRF")
  MAXDIM <- integer(1)
  .C("GetrfParameters", integer(1), integer(1), integer(1), integer(1),
     integer(1), integer(1), MAXDIM, integer(1),
     PACKAGE="RandomFields", DUP=FALSE)
  dimensions <- .C("GetKeyInfo", as.integer(register), 
                   total = integer(1), len = integer(MAXDIM),
                   spatialdim = integer(1),
                   timespacedim = integer(1),
                   grid = integer(1), distr = integer(1),
                   MAXDIM = as.integer(MAXDIM),
                   PACKAGE="RandomFields", DUP=FALSE)
  if (dimensions$total <= 0)
    stop(paste("register", register, "does not look initialized"))
  dimensions$len <- dimensions$len[1:dimensions$timespacedim]
  storage.mode(register) <- "integer"
  error <- integer(1)
  
  ## trend not used yet!
  modus <- integer(1)
  lt <- integer(1)
  ll <- integer(1)
  lx <- integer(1)
  .C("GetTrendLengths", register, modus, lt, ll, lx,
     PACKAGE="RandomFields", DUP=FALSE)
  if(modus<0) stop(paste("Error in trend storage", -modus))
  ## modus==0 : no trend, only mean
  if (modus>0) {
    stop("not programmed yet!!")
    trend <- .C("GetTrend", register, trend=paste(rep(" ",lt+1),collapse=""),
                lambda=double(ll), error=integer(1), PACKAGE="RandomFields")
    if (trend$error) return(NULL)
    lambda <- trend$lambda
    trend <- eval(parse(text=trend$trend))
    x <- double(lx)
    .C("GetCoordinates", register, x, error, PACKAGE="RandomFields", DUP=FALSE)
    if (error) return(NULL)
  }
  
  if (n == 1) {
    result <- double(dimensions$total)

    .C(DoNameList[1+dimensions$distr], register,
       result, error, PACKAGE="RandomFields", DUP=FALSE)
    if (error) return(NULL)
    else {
      if (dimensions$grid) 
        result <- array(result, dim = dimensions$len)
      else if (dimensions$spatialdim<dimensions$timespacedim) {
        result <- array(result, dim=dimensions$len[c(1,length(dimensions$len))])
      }
    }
  } else { ## n>1
    param <- RFparameters()[c("Storing", "pch")]
    if (!param$Storing) RFparameters(Storing=TRUE)
    if (dimensions$grid) {
      result <- array(NA, dim = c(dimensions$len, n))
      s <- paste(rep(",", dimensions$timespacedim), collapse="")
    }
    else if (dimensions$spatialdim==dimensions$timespacedim) {
      result <- matrix(nrow = dimensions$total, ncol = n)
      s <- ","
    } else {
      result <- array(NA, dim=c(dimensions$len[c(1,length(dimensions$len))], n))
      s <- ",,"
    }
    s <- paste("result[", s, "i] <- res")
    res <- double(dimensions$total)
     for (i in 1:n) {
      cat(param$pch)
      .C(DoNameList[1+dimensions$distr], register,
         res, error, PACKAGE="RandomFields", DUP = FALSE, NAOK = TRUE)
      if (error) break   
      eval(parse(text = s))
    }
    if (!param$Storing) {
      RFparameters(Storing=FALSE)
      DeleteRegister(register)
    }
    if (param$pch!="") cat("\n")
    if (error) return(NULL)
  }
  ## to do: add trend !!
  return(result)
}


".First.lib" <- function (lib, pkg) 
{
  if (RNGkind()[1]!="Mersenne-Twister") {
    cat("Random number generator changed to 'Mersenne-Twister'\n")
    RNGkind(kind = "Mersenne-Twister")
  }
  if (file.exists("/home/schlather/bef/x")) {
    ## to do list -- since my surname is rare, the message should 
    ## appear only on computers I have a login
    cat("To-Do List\n==========\n")
    print("include winddata")
    print("mleRF.Rd/wind.Rd/GaussRF.Rd examples fertig machen")
    
    print("")
    print("check documentation and readability of programs")
    print("docu of TBM2.linesimustep not understandable")
    print("entartete Felder")
    print("fractGauss -- was ist los?")
    print("2d/3d fractal testen!")
    print("clean up tbm.cc: trennung circ embed und direct")
    cat("\n\n")
    
    print("MLE: naturalscaling in anisotropic case")
    print("critical odd unused, see RFcircembed.cc; see also addodd in RFgetNset.cc")
    print("implement trend; implement REML")
    print("MPP.cc: anisotropies, time")
    print("MLE:  message innerhalb der Liste, dass irgendwas nicht in Ordnung falls eine Grenze angenommen wurde (error nr und error message, error message auch anzeigen, falls printlevel hoch genug)")
    cat("individuelle Praeferenzliste:\n",
        "   1) erste praeferenz\n",
        "   2) if # pkte < 10 / 50 -> Gauss??\n",
        "   4) if nugget/variance > 0.01 -> CE (and not spectral, e.g.)\n",
        "   3) if C(maxh)/C(0) > 0.05  -> TBM else CE\n")
    cat("spectral and other space sensitive procedure work curently only for",
        "zonal isotropy of arbitrary dim (if final dim==2) and arbitrary",
        "transformation in 2dim, but not for arbitrary transformation with",
        "final dim == 2\n")
    cat("register inhalt auf platte speichen\n")
    cat("Aufpassen, dass kein Chaos produziert wird, wenn InitPoissonRF und",
        "dann DoGaussRF fabriziert wird. dies wird abgefangen auf R-Ebene; ",
        "aber zu ueberdenken waere, ob man ein 'Mischen' nicht auch zulaesst\n")
    cat("init poisson : check mean=variance; in rf.R: either one is NaN, or equal, set both equal\n")
    cat("interface, such that user can add its own covariance function, written in R\n\n")

  }
  library.dynam("RandomFields", pkg, lib)
}


"InitGaussRF" <-
  function(x, y = NULL, z = NULL, T=NULL, grid, model, param,
           trend, method = NULL,
           register = 0, gridtriple = FALSE)
  {
  return(InitSimulateRF(x=x, y=y, z=z,  T=T, 
                        grid=grid, model=model, param=param, trend=trend,
                        method=method, register=register,
                        gridtriple=gridtriple, distribution="Gauss"))
}


"InitSimulateRF" <-
function (x, y = NULL, z = NULL, T=NULL, grid, model, param,
          trend, method = NULL,
          register = 0, gridtriple = FALSE, distribution=NA)  
{
  modus <-  0
  if (!missing(trend) && (!is.numeric(trend) || (length(trend)!=1))) {
    stop("sorry; not programmed yet.")
    if (is.list(trend)) {
      lambda <- trend$lambda
      stopifnot(is.numeric(lambda))
      trend <- trend$trend
      modus <- 3
    } else if (is.numeric(trend)) {
      lambda <- trend
      trend <-  ""
      modus <- 1
    } else {
      lambda <- 0
      trend <- as.character(substitute(trend))
      trend <- trend[length(trend)]
      modus <- 2
    }
    error <-  .C("StoreTrend",
                 as.integer(register),
                 as.integer(modus),
                 ## 0 : mean
                 ## 1 : linear trend
                 ## 2 : function
                 ## 3 : function with parameters
                 as.character(trend), as.integer(nchar(trend)),
                 as.double(lambda), as.integer(length(lambda)),
                 err=integer(1), PACKAGE="RandomFields")
    if (error) stop(paste("trend not correct -- error nr", error))
  } else {
    ## trend is indeed the mean
    ## dummy storage; indeed only a flag is set that StoreTrend had
    ## been called
    .C("StoreTrend",as.integer(register), as.integer(0), as.character(""),
       as.integer(0), double(0), as.integer(0), err=integer(1),
       PACKAGE="RandomFields")
  }

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
  new <- CheckXT(x, y, z, T, grid, gridtriple)
  p <- PrepareModel(model, param, new$spacedim+new$Time, trend=trend,
                    method=method)
  if (any(is.na(p$param))) stop("some model parameters are NA") 
  
  error <- max(0,.C(InitName, as.double(new$x), as.double(new$T),
                  as.integer(new$spacedim),
                  as.integer(new$l), 
                  as.integer(new$grid),
                  as.integer(new$Time),
                  as.integer(p$covnr),
                  as.double(p$param),
                  as.integer(length(p$param)),
                  as.double(p$mean),
                  as.integer(length(p$covnr)),
                  as.integer(p$anisotropy),
                  as.integer(p$op),
                  as.integer(p$method),
                  as.integer(distrNr),
                  as.integer(register),
                  error=integer(1),
                  PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)$error)
   return(error)
}



"GaussRF" <-
function (x, y = NULL, z = NULL, T=NULL,
          grid, model, param, trend, method = NULL, 
          n = 1, register = 0, gridtriple = FALSE) 
{
  if (n>1 && !(storing <- RFparameters()$Storing)) RFparameters(Storing=TRUE)
  error <- InitSimulateRF(x=x, y=y, z=z, T=T, grid=grid, model=model,
                          param=param,
                          trend=trend, method=method, register=register,
                          gridtriple=gridtriple, distribution="Gauss") > 0
  if (n>1 && !storing) RFparameters(Storing=FALSE)
  if (error) return(NULL)
  return(DoSimulateRF(n=n, reg=register))
}


"Variogram" <-
function (x, model, param, dim=ifelse(is.matrix(x),ncol(x),1))
  CovarianceFct(x, model, param, dim, fctcall="Variogram")


pokeTBM <- function(Out, In) {
  .C("pokeTBM", as.integer(Out), as.integer(In), err=integer(1),
     PACKAGE="RandomFields")$err
}

