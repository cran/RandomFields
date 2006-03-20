## source("rf.R")
# see getNset.R for the variable .methods

# readline <- function(...) return("") 


.onLoad <- function (lib, pkg) {
  if (file.exists("/home/schlather/bef/x")) {
    ## to do list -- since my surname is rare, the message should 
    ## appear only on computers I have a login
    cat("To-Do List\n==========\n")

    print("RFdirect: zusaetzlicher Parameter, der regelt ob SVD Zerlegung ueberprueft ob wirklich OK")
    print("RFparameters: set aufrufe zu einem + interne environment darstellung")
    print("Showmodels: die mdoelle die nicht funktioniern in schwaecherer Farbe")
    print("nsst/nsst2 als hypermodel formulieren")
    print("C-1-Berechnung bei MLE durch FFT ersetzen; mit Erweiterung von C auf circulaere Matrix")
    print("cutoff/intrinsic fuer Produkte / zonale Anisotropie bei Produkten")
    print("Bsp in GaussRF und RFparameters muessen von Hand durchgeschaut werden,
ob sie den richtigen Effekt zeigen!")
    print("environment statt der einzelnen Variablen; vereinfachung von RFparameters und der internen initGaussRD und DoSimulate");
    print("hyper: was ist die truetimespacedim???")
    print("hyper: 'Ma'-Modelle; TBM3/2;");
    print("cross: verbesserungen (geswchwindigkeit)")
    print("GENERAL_PRECISION: einbinden + ueber .Mschine$precision definieren (mal Faktor 50")
    print("chlo2inv in MLE: ")
    print("Empirical Variogram: allow for NAs")
    print("include winddata")
    print("fitvario.Rd/wind.Rd/GaussRF.Rd examples fertig machen")
    
    print("check documentation and readability of programs")
    print("MLE: naturalscaling in anisotropic case")
    print("critical odd unused, see RFcircembed.cc; see also addodd in RFgetNset.cc")
    print("implement trend")
    print("MPP.cc: anisotropies, time")
    cat("init poisson : check mean=variance; in rf.R: either one is NaN, or equal, set both equal\n")
    cat("interface, such that user can add its own covariance function, written in R\n\n")
  }
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
  x <- t(x)  ## C code expects the coordinates columwise
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


DoSimulateRF <- function (n = 1, register = 0, paired=FALSE) {
  stopifnot(length(n) == 1,
            n>0, is.finite(n),
            length(register) == 1,
            is.finite(register)
            )
  storage.mode(register) <- "integer"
  DoNameList <- c("DoSimulateRF", "DoSimulateRF", "DoMaxStableRF")
  assign(".p",
         .C("GetrfParameters", covmaxchar=integer(1), methodmaxchar=integer(1),
            distrmaxchar=integer(1),
            covnr=integer(1), methodnr=integer(1), distrnr=integer(1),
            maxdim=integer(1), maxmodels=integer(1),
            PACKAGE="RandomFields"))    
  keyinfo <- .C("GetKeyInfo", register, 
                total = integer(1), len = integer(.p$maxdim),
                spatialdim = integer(1), timespacedim = integer(1),
                grid = integer(1), distr = integer(1),
                maxdim = as.integer(.p$maxdim),
                PACKAGE="RandomFields", DUP=FALSE)
  if (paired && (n %% 2 != 0)) stop("if paired, then n must be an even number")
  if (keyinfo$total <= 0)
    stop(paste("register", register, "does not look initialized"))
  keyinfo$len <- keyinfo$len[1:keyinfo$timespacedim]
  error <- integer(1)
    
  result <- double(keyinfo$total * n)
  .C(DoNameList[1+keyinfo$distr], register, as.integer(n),
     as.integer(paired), result, error, PACKAGE="RandomFields", DUP=FALSE)
  
  if (error) stop(paste("error", error));
  if (n==1) {
    if (keyinfo$grid) {
      if (length(keyinfo$len)>1) dim(result) <- keyinfo$len
    } else if (keyinfo$spatialdim<keyinfo$timespacedim)
      dim(result) = keyinfo$len[c(1,length(keyinfo$len))]
  } else {
    dim(result) <-
      c(if (keyinfo$grid) keyinfo$len else
        if (keyinfo$spatialdim==keyinfo$timespacedim) keyinfo$total else
        keyinfo$len[c(1,length(keyinfo$len))]
        , n)
  }
  
  return(result)
}



"InitSimulateRF" <-
function (x, y = NULL, z = NULL, T=NULL, grid, model, param,
          trend, method = NULL,
          register = 0, gridtriple = FALSE, distribution=NA)  
{
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
  
  ## modus:
  ## 0 : mean
  ## 1 : linear trend
  ## 2 : function
  ## 3 : function with parameters
  if (is.null(p$trend)) {
    lambda <- p$mean
    trend <- ""    
    modus <- 0
  } else {
    stop("sorry; not programmed yet.")
    if (is.list(p$trend)) {
      stopifnot(is.numeric(lambda <- p$trend$lambda))
      trend <- p$trend$trend
      modus <- 3
    } else
    if (is.numeric(trend)) {
      lambda <- p$trend
      trend <- ""    
      modus <- 1
    } else { # function
      lambda <- 0
      trend <- as.character(substitute(trend))
      trend <- trend[length(trend)]
      modus <- 2
    }
  }
  
  ## "StoreTrend must be called before .C(InitName, ...)" !!
  ## see extrems.cc for example
  error <- .C("StoreTrend", as.integer(register), as.integer(modus),
              as.character(trend), as.integer(nchar(trend)),
              as.double(lambda), as.integer(length(lambda)),
              err=integer(1), PACKAGE="RandomFields")$err
  if (error) stop(paste("trend not correct -- error nr", error))

  error <- .C(InitName, as.double(new$x), as.double(new$T),
                  as.integer(new$spacedim),
                  as.integer(new$l), 
                  as.integer(new$grid),
                  as.integer(new$Time),
                  as.integer(p$covnr),
                  as.double(p$param),
                  as.integer(length(p$param)),
                  as.integer(length(p$covnr)),
                  as.integer(p$anisotropy),
                  as.integer(p$op),
                  as.integer(p$method),
                  as.integer(distrNr),
                  as.integer(register),
                  error=integer(1),
                  PACKAGE="RandomFields", DUP = FALSE, NAOK=TRUE)$error
  return(error)
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


"GaussRF" <-
function (x, y = NULL, z = NULL, T=NULL,
          grid, model, param, trend, method = NULL, 
          n = 1, register = 0, gridtriple = FALSE,
          paired=FALSE, ...) 
{
  old.param <- RFparameters(no.readonly=TRUE)
  RFpar <- list(...)
  if (length(RFpar)>0) RFparameters(RFpar)
  if (delete <- n>1 && !RFparameters()$Storing) RFparameters(Storing=TRUE)
  on.exit({if (delete) DeleteRegister(register);
           RFparameters(old.param);
           })

  error <- InitSimulateRF(x=x, y=y, z=z, T=T, grid=grid, model=model,
                          param=param,
                          trend=trend, method=method, register=register,
                          gridtriple=gridtriple, distribution="Gauss")
#str(GetRegisterInfo(register))
  
  if (error > 0)
    stop(paste("Simulation could not be initiated.",
               if (RFparameters()$PrintLevel >= 2) "\nRerun with higher value of RFparameters()$PrintLevel for more information. (Or put debug=TRUE if you are using Showmodels.)\n\n"))
  return(if (n<1) NULL else DoSimulateRF(n=n, reg=register, paired=paired))
}


"Variogram" <-
function (x, model, param, dim=ifelse(is.matrix(x),ncol(x),1))
  CovarianceFct(x, model, param, dim, fctcall="Variogram")

