"CovarianceFct" <-
function (x, model, param, dim = 1) 
{
  if (!is.finite(param[1])) {
        param[1] <- 0
    }
  stopifnot(length(x)!=0,
            all(is.finite(x)),
            is.character(model),
            all(is.finite(param)),
            length(dim) == 1,
            is.finite(dim)
            )
  
  covnr <- as.integer(.C("GetModelNr", as.character(model), nr = integer(1))$nr)
    if (covnr < 0) {
        .C("PrintModelList")
        stop("given model cannot be (uniquely) identified from the above list")
    }
    kappas <- as.integer(.C("GetNrParameters", covnr, k = integer(1),DUP=FALSE)$k)
    if (length(param) < 4 + kappas) 
        stop("not enough parameters for the covariance function")
    storage.mode(x) <- "double"
    storage.mode(param) <- "double"
    storage.mode(dim) <- "integer"
    return(.C("Covariance", x, as.integer(length(x)), covnr, param, 
              as.integer(length(param)), dim, res=double(length(x)), DUP = FALSE)$res
           )
}
"DeleteAllRegisters" <-
function () 
{
    .C("DeleteAllKeys")
    return(NULL)
}
"DeleteRegister" <-
function (nr = 0) 
{
   stopifnot( length(nr) == 1, is.finite(nr) ) 
    .C("DeleteKey", as.integer(nr))
    return(NULL)
}
"DoSimulateRF" <-
function (n = 1, register = 0) 
{
  stopifnot(length(n) == 1,
            is.finite(n),
            length(register) == 1,
            is.finite(register)
            )
    DoNameList <- c("DoSimulateRF","DoSimulateRF","DoMaxStableRF")    
    MAXDIM <- 3
    dimensions <- .C("GetKeyInfo", as.integer(register), 
        total = integer(1), len = integer(MAXDIM), dim = integer(1), 
        grid = integer(1), distr = integer(1), DUP=FALSE)
    param <- RFparameters()
    if (dimensions$total <= 0) {
        stop(paste("register", register, "does not look initialized"))
    }
    if (n == 1) {
      x <- .C(DoNameList[1+dimensions$distr], as.integer(register),
                 result=double(dimensions$total), error=integer(1),DUP=FALSE)
      if (x$error) 
        return(NULL)
      else if (dimensions$grid) 
        return(array(x$result, dim = dimensions$len[1:dimensions$dim])) 
      else return(x$result)
    }
    else {
        if (!param$Storing) 
            stop("Use option `RFparameters(Storing=TRUE)' if n>1")

        if (dimensions$grid) {
            result <- array(NA, dim = c(dimensions$len[1:dimensions$dim], 
                n))
            s <- paste(rep(",",dimensions$dim),collapse="")
          }
        else {
          result <- matrix(nrow = dimensions$total, ncol = n)
          s <- ","
        }
        s <- paste("result[", s, "i] <- dummy$res")
        pch <- param$pch
        for (i in 1:n) {
           cat(pch)
           dummy <- .C(DoNameList[1+dimensions$distr], as.integer(register),
                       res=double(dimensions$total),
                       error=integer(1), DUP = FALSE, NAOK = TRUE)
            if (dummy$error) 
                break
            eval(parse(text = s))
        }
        if (pch!="") cat("\n")
        if (dummy$error) 
            return(NULL)
        else return(result)
    }
}
"EmpiricalVariogram" <-
function (x, y = NULL, z = NULL, data, grid, bin, gridtriple = FALSE) 
{
    stopifnot(all(is.finite(data)),
            length(bin)>=2,
            all(is.finite(bin)),
              bin[1] <= 0
            )
  
    new <- CheckXYZ(x, y, z, grid, gridtriple)
    
    repet <- as.integer(length(data)/new$total)
    if (length(data) != new$total * repet) 
        stop("number of data does not match coordinates")
    
    centers <- pmax(0,(bin[-1] + bin[-length(bin)])/2)
    n.bins <- length(bin) - 1 
    emp.vario <- double(n.bins)
    n.bin <- integer(n.bins)
    
    .C("empiricalvariogram",
       as.double(new$x), as.double(new$y), as.double(new$z),
       as.integer(new$dim), as.integer(new$l), 
       as.double(data), as.integer(repet), as.integer(grid), 
       as.double(bin), as.integer(n.bins), as.integer(0), 
       emp.vario, n.bin, DUP = FALSE)
    emp.vario[is.na(emp.vario) & (centers==0)] <- 0
    return(list(centers=centers, emp.vario=emp.vario, n.bin=n.bin))
}
".First.lib" <- function (lib, pkg) 
{
    library.dynam("RandomFields", pkg, lib)
    library(mva)
}
"InitGaussRF" <-
  function(x, y = NULL, z = NULL, grid, model, param, method = NULL, 
    register = 0, gridtriple = FALSE) 
{
  return(InitSimulateRF(x=x, y=y, z=z, grid=grid, model=model,
                        param=param, method=method, register=register,
                        gridtriple=gridtriple, distribution="Gauss"))
}

CheckXYZ <- function(x, y, z, grid, gridtriple){
  stopifnot(length(x) != 0,
            is.logical(grid),
            is.logical(gridtriple)
            )
  if (is.data.frame(x)) {
    if (ncol(x)==1) x <- as.vector(x) else x <- as.matrix(x)
  }
  stopifnot(all(is.finite(x))) 
  if (is.matrix(x)) {
        if (!is.null(y) || !is.null(z)) {
            stop("If x is a matrix, then y and z may not be given")
        }
        dim <- ncol(x)
        if (dim > 3) {
            dim <- 3
            warning("only the first three columns are considered as coordinates")
        }
        if (dim > 1) {
            y <- x[, 2]
            if (dim > 2) 
                z <- x[, 3]
        }
        x <- x[, 1]
        l <- length(x)
    }
    else {
        l <- length(x)
        storage.mode(x) <- "double"
        if (is.null(y)) {
            dim <- 1
            if (!is.null(z)) {
                stop("y is not given, but z")
            }
        }
        else {
            if ((l != length(y)) && (!grid || gridtriple)) {
                stop("x and y differ in length")
            }
            stopifnot(all(is.finite(y)))
            storage.mode(y) <- "double"
            if (is.null(z)) {
                dim <- 2
            }
            else {
                dim <- 3
                if (l != length(z) && (!grid || gridtriple)) {
                  stop("x and z differ in length")
                }
                stopifnot(all(is.finite(z))) 
                storage.mode(z) <- "double"
            }
        }
    }
    if (l == 1) 
        stop("Use grid=FALSE if only a single point is simulated")
    else if (grid) {
        if (gridtriple) {
            if (l != 3) {
                stop("In case of simulating a grid with option gridtriple, exactly 3 numbers are needed for each direction")
            }
         lx <- length(seq(x[1],x[2],x[3]))
         ly <- lz <- 1
         x[2] <- x[1] + (lx - 0.999) * x[3]
         if (dim > 1) {
             ly <- length(seq(y[1],y[2],y[3]))
             y[2] <- y[1] + (ly - 0.999) * y[3]
             if (dim > 2) {
               lz <- length(seq(z[1],z[2],z[3]))
               z[2] <- z[1] + (lz - 0.999) * z[3]
              }
           }
         total <- lx * ly * lz
       if (total==0) stop("incorrect grid specification (one or no points)")
        
        }
        else {
            eqdist <- function(x) {
                step <- diff(x)                
                if (sum(abs(step - step[1]))/step[1] > 1e-10) {
                  stop("Grid must have equal distances in each direction.")
                }
                return(c(x[1], x[length(x)]+0.001*step[1], step[1]))
            }
            total <- length(x)
            x <- eqdist(x)            
            if (dim > 1) {
                total <- total * length(y)
                y <- eqdist(y)
                if (dim > 2) {
                  total <- total * length(z)
                  z <- eqdist(z)
                }
            }
            l <- 3
        }
    }
    else {
        total <- length(x)
    }
  return(list(x=x, y=y, z=z, total=total, l=l, dim=dim))
}

"InitSimulateRF" <-
function (x, y = NULL, z = NULL, grid, model, param, method = NULL, 
          register = 0, gridtriple = FALSE,distribution=NA)  
{
   distributionList <- c("Gauss","Poisson","MaxStable")
   InitNameList <- c("InitSimulateRF","InitSimulateRF","InitMaxStableRF")
   if (is.na(distribution)) {
     stop("This function is an internal function.\nUse `GaussRF', `InitGaussRF', `MaxStableRF', etc., instead.\n")    
   }
   distrNr <- pmatch(distribution,distributionList)   
   if (is.na(distrNr)) stop("Unknown distribution -- do not use this internal function directly")
   else {
     if ((distrNr==2) && !exists("PoissonRF")) stop("Sorry. Not programmed yet.")
     InitName <- InitNameList[distrNr]
     distrNr <- distrNr - 1;
   }
   stopifnot(is.character(model)) 
    covnr <- as.integer(.C("GetModelNr", as.character(model), 
        nr = integer(1))$nr)
    if (covnr < 0) {
        .C("PrintModelList")
        stop("given model cannot be (uniquely) identified from the above list")
    }
    if (is.null(method)) {
        method <- -1
    }
    else {
        if (!is.character(method)) {
            stop("method must be NULL or a string")
        }
        method <- as.integer(.C("GetMethodNr", method, nr = integer(1))$nr)
        if (method < 0) {
          .C("PrintMethods")
          stop("given method cannot be (uniquely) identified from the above list")
        }
    }
    kappas <- as.integer(.C("GetNrParameters", covnr, k = integer(1))$k)
    if (length(param) < 4 + kappas) 
        stop("not enough parameters for the covariance function")
    storage.mode(param) <- "double"

    new <- CheckXYZ(x, y, z, grid, gridtriple)
   
    return(
           max(0,.C(InitName, as.double(new$x), as.double(new$y),
                    as.double(new$z), as.integer(new$dim), as.integer(new$l), 
              as.integer(grid), covnr, param, as.integer(length(param)), 
              as.integer(method), as.integer(distrNr),
                    as.integer(register), error=integer(1),
              DUP = FALSE)$error)
    )
}
"PrintModelList" <-
function () 
{
    .C("PrintModelList")
    return(NULL)
}
"PrintMethodList" <-
function () 
{
    .C("PrintMethods")
    return(NULL)
}
"RFparameters" <-function (...) {
  RFparameters.default(...)
}
"RFparameters.default" <-
function (Storing = storing, PrintLevel = printlevel,
          PracticalRange = practicalrange, 
    CE.force = ce.force, CE.mmin = ce.mmin, CE.tolRe = ce.tolRe, 
    CE.tolIm = ce.tolIm, CE.trials = ce.trials,
    direct.checkprecision = directcheckprecision, 
    direct.maxvariables = directmaxvariables,
    direct.method = directmethod, 
    direct.requiredprecision = directrequiredprecision,
    spectral.lines = spectrallines, 
    spectral.grid = spectralgrid, TBMCE.force = tbmceforce, TBMCE.mmin = tbmcemmin, 
    TBMCE.tolRe = tbmcetolre, TBMCE.tolIm = tbmcetolim, TBMCE.trials = tbmcetrials, 
    TBM2.lines = tbm2lines, TBM2.linesimufactor = tbm2linesimufactor, 
    TBM2.linesimustep = tbm2linesimustep, TBM3D2.lines = tbm3D2lines, 
    TBM3D2.linesimufactor = tbm3D2linesimufactor, TBM3D2.linesimustep = tbm3D2linesimustep, 
    TBM3D3.lines = tbm3D3lines, TBM3D3.linesimufactor = tbm3D3linesimufactor, 
    TBM3D3.linesimustep = tbm3D3linesimustep,
    MPP.approxzero=mppapproxzero, add.MPP.realisations=addmpprealisations,
    MPP.radius=mppradius,
    maxstable.maxGauss=maxstablemaxGauss,
    pch = pchx
 ) {
    x <- .C("SetParam", as.integer(1), storing = integer(1), 
            printlevel = integer(1), practicalrange = integer(1),
            pch="  ")
    storing <- x$storing
    printlevel <- x$printlevel
    practicalrange <- as.logical(x$practicalrange)
    if (!is.finite(Storing + PrintLevel + PracticalRange)) 
        stop("some parameters are not finite")
    pchx <- x$pch
    if (!is.character(pch)) stop("pch is not a character")
   
    ## do not allow integer values for users!
    PracticalRange <- as.logical(PracticalRange)
    .C("SetParam", as.integer(0), as.integer(Storing), as.integer(PrintLevel), 
        as.integer(PracticalRange), pch)
    
    x <- .C("SetParamCircEmbed", as.integer(1), force = integer(1), 
        tolRe = double(1), tolIm = double(1), trials = integer(1), 
        mmin = integer(1))
    ce.force <- x$force
    ce.tolRe <- x$tolRe
    ce.tolIm <- x$tolIm
    ce.trials <- x$trials
    ce.mmin <- x$mmin
    if (!is.finite(CE.force + CE.mmin + CE.tolRe + CE.tolIm + 
        CE.trials)) 
        stop("some parameters are not finite")
    .C("SetParamCircEmbed", as.integer(0), as.integer(CE.force), 
        as.double(CE.tolRe), as.double(CE.tolIm), as.integer(CE.trials), 
        as.integer(CE.mmin))
    x <- .C("SetParamTBMCE", as.integer(1), tbmceforce = integer(1), 
        tbmcetolre = double(1), tbmcetolim = double(1), tbmcetrials = integer(1), 
        tbmcemmin = integer(1))
    tbmceforce <- x$tbmceforce
    tbmcetolre <- x$tbmcetolre
    tbmcetolim <- x$tbmcetolim
    tbmcetrials <- x$tbmcetrials
    tbmcemmin <- x$tbmcemmin
    if (!is.finite(TBMCE.force + TBMCE.mmin + TBMCE.tolRe + TBMCE.tolIm + 
        TBMCE.trials)) 
        stop("some parameters are not finite")
    .C("SetParamTBMCE", as.integer(0), as.integer(TBMCE.force), 
        as.double(TBMCE.tolRe), as.double(TBMCE.tolIm), as.integer(TBMCE.trials), 
        as.integer(TBMCE.mmin))
    x <- .C("SetParamTBM2", as.integer(1), tbm2lines = integer(1), 
        tbm2linesimufactor = double(1), tbm2linesimustep = double(1))
    tbm2lines <- x$tbm2lines
    tbm2linesimufactor <- x$tbm2linesimufactor
    tbm2linesimustep <- x$tbm2linesimustep
    if (!is.finite(TBM2.lines + TBM2.linesimufactor + TBM2.linesimustep)) 
        stop("some parameters are not finite")
    .C("SetParamTBM2", as.integer(0), as.integer(TBM2.lines), 
        as.double(TBM2.linesimufactor), as.double(TBM2.linesimustep))
    x <- .C("SetParamTBM3D2", as.integer(1), tbm3D2lines = integer(1), 
        tbm3D2linesimufactor = double(1), tbm3D2linesimustep = double(1))
    tbm3D2lines <- x$tbm3D2lines
    tbm3D2linesimufactor <- x$tbm3D2linesimufactor
    tbm3D2linesimustep <- x$tbm3D2linesimustep
    if (!is.finite(TBM3D2.lines + TBM3D2.linesimufactor + TBM3D2.linesimustep)) 
        stop("some parameters are not finite")
    .C("SetParamTBM3D2", as.integer(0), as.integer(TBM3D2.lines), 
        as.double(TBM3D2.linesimufactor), as.double(TBM3D2.linesimustep))
    x <- .C("SetParamTBM3D3", as.integer(1), tbm3D3lines = integer(1), 
        tbm3D3linesimufactor = double(1), tbm3D3linesimustep = double(1))
    tbm3D3lines <- x$tbm3D3lines
    tbm3D3linesimufactor <- x$tbm3D3linesimufactor
    tbm3D3linesimustep <- x$tbm3D3linesimustep
    if (!is.finite(TBM3D3.lines + TBM3D3.linesimufactor + TBM3D3.linesimustep)) 
        stop("some parameters are not finite")
    .C("SetParamTBM3D3", as.integer(0), as.integer(TBM3D3.lines), 
        as.double(TBM3D3.linesimufactor), as.double(TBM3D3.linesimustep))
    x <- .C("SetParamSpectral", as.integer(1), spectrallines = integer(1), 
        spectralgrid = integer(1))
    spectrallines <- x$spectrallines
    spectralgrid <- x$spectralgrid
    if (!is.finite(spectral.lines + spectral.grid)) 
        stop("some parameters are not finite")
    .C("SetParamSpectral", as.integer(0), as.integer(spectral.lines), 
        as.integer(spectral.grid))
    x <- .C("SetParamDirectGauss", as.integer(1), directmethod = integer(1), 
        directcheckprecision = integer(1), directrequiredprecision = double(1), 
        directmaxvariables = integer(1))
    directmethod <- x$directmethod
    directcheckprecision <- x$directcheckprecision
    directrequiredprecision <- x$directrequiredprecision
    directmaxvariables <- x$directmaxvariables
    if (!is.finite(direct.checkprecision + direct.maxvariables + 
        direct.method + direct.requiredprecision)) 
        stop("some parameters are not finite")
    .C("SetParamDirectGauss", as.integer(0), as.integer(direct.method), 
        as.integer(direct.checkprecision), as.double(direct.requiredprecision), 
        as.integer(direct.maxvariables))

    x <- .C("SetMPP", as.integer(1),
            mppapproxzero=double(1), addmpprealisations=double(1),
            mppradius=double(1))
    mppapproxzero <- x$mppapproxzero;
    addmpprealisations <- x$addmpprealisations;
    mppradius <- x$mppradius;
    if (!is.finite(MPP.approxzero+add.MPP.realisations+MPP.radius)) 
        stop("some parameters are not finite")
    .C("SetMPP", as.integer(0),
       as.double(MPP.approxzero),as.double(add.MPP.realisations),
       as.double(MPP.radius))
  
    
    x <- .C("SetExtremes", as.integer(1),
            maxstablemaxGauss=double(1))
    maxstablemaxGauss <- x$maxstablemaxGauss;
    if (!is.finite(maxstable.maxGauss)) 
        stop("some parameters are not finite")
    .C("SetExtremes", as.integer(0),
       as.double(maxstable.maxGauss))
     
    if (length(as.list(match.call()))>1) return(NULL)
    else return(list(Storing = as.logical(Storing), PrintLevel = PrintLevel, 
        PracticalRange = as.logical(PracticalRange),
        CE.force = as.logical(CE.force), CE.mmin = CE.mmin,
        CE.tolRe = CE.tolRe, CE.tolIm = CE.tolIm, CE.trials = CE.trials, 
        direct.checkprecision = as.logical(direct.checkprecision), 
        direct.maxvariables = direct.maxvariables,
        direct.method = direct.method,
        direct.requiredprecision = direct.requiredprecision, 
         spectral.lines = spectral.lines, 
        spectral.grid = as.logical(spectral.grid),
        TBMCE.force = as.logical(TBMCE.force), TBMCE.mmin = TBMCE.mmin, 
        TBMCE.tolRe = TBMCE.tolRe, TBMCE.tolIm = TBMCE.tolIm, 
        TBMCE.trials = TBMCE.trials, 
        TBM2.lines = TBM2.lines, TBM2.linesimufactor = TBM2.linesimufactor, 
        TBM2.linesimustep = TBM2.linesimustep, TBM3D2.lines = TBM3D2.lines, 
        TBM3D2.linesimufactor = TBM3D2.linesimufactor, TBM3D2.linesimustep = TBM3D2.linesimustep, 
        TBM3D3.lines = TBM3D3.lines, TBM3D3.linesimufactor = TBM3D3.linesimufactor, 
        TBM3D3.linesimustep = TBM3D3.linesimustep,
        MPP.approxzero=MPP.approxzero, add.MPP.realisations=add.MPP.realisations,
        MPP.radius=MPP.radius,
        maxstable.maxGauss=maxstable.maxGauss,
        pch=pch
        ))
}
"GaussRF" <-
function (x, y = NULL, z = NULL, grid, model, param, method = NULL, 
    n = 1, register = 0, gridtriple = FALSE) 
{
    if (InitSimulateRF(x=x, y=y, z=z, grid=grid, model=model,
                       param=param, method=method, register=register,
                       gridtriple=gridtriple, distribution="Gauss") <= 0) {
        return(DoSimulateRF(n=n, reg=register))
    }
    else {
        return(NULL)
    }
}
"SimulateRF" <-
function (x, y = NULL, z = NULL, grid, model, param, method = NULL, 
    n = 1, register = 0, gridtriple = FALSE, distribution=NA) 
{
    if (InitSimulateRF(x=x, y=y, z=z, grid=grid, model=model,
                       param=param, method=method, 
                       reg=register, gridtriple=gridtriple,
                       distr=distribution) <= 0) {
        return(DoSimulateRF(n=n, reg=register))
    }
    else {
        return(NULL)
    }
}
"Variogram" <-
function (x, model, param, dim = 1) 
{
  if (!is.finite(param[1])) param[1] <- 0
  stopifnot(length(x) != 0,
            all(is.finite(x)),
            is.character(model),
            all(is.finite(param)),
            length(dim) == 1,
            is.finite(dim)) 
    covnr <- as.integer(.C("GetModelNr", as.character(model), nr = integer(1))$nr)
    if (covnr < 0) {
        .C("PrintModelList")
        stop("given model cannot be (uniquely) identified from the above list")
    }
    kappas <- as.integer(.C("GetNrParameters", covnr, k = integer(1))$k)
    if (length(param) < 4 + kappas) 
        stop("not enough parameters for the variogram model")
    storage.mode(x) <- "double"
    storage.mode(param) <- "double"
    storage.mode(dim) <- "integer"
    return(
           .C("Variogram", x, as.integer(length(x)), covnr, param, as.integer(length(param)), 
              dim, res=double(length(x)), DUP = FALSE)$res
           )
}

GetModelNames <- function() {
  p <- .C("GetrfParameters",covmaxchar=integer(1),methodmaxchar=integer(1),
          covnr=integer(1),methodnr=integer(1))
  l <- character(p$covnr)
  for (i in 1:p$covnr) {
    l[i] <- .C("GetModelName",as.integer(i-1),
               n=paste(rep(" ",p$covmaxchar),collapse=""))$n
  }
  return (l)
}


GetMethodNames <- function() {
  p <- .C("GetrfParameters",covmaxchar=integer(1),methodmaxchar=integer(1),
          covnr=integer(1),methodnr=integer(1))
  l <- character(p$methodnr)
  for (i in 1:p$methodnr) {
    l[i] <- .C("GetMethodName",as.integer(i-1),
               n=paste(rep(" ",p$methodmaxchar),collapse=""))$n
  }
  return(l)
}


GetPracticalRange <- function(model,kappas=NULL) {
  covnr <- as.integer(.C("GetModelNr", as.character(model),
                         nr = integer(1))$nr)
  if (covnr < 0) {
    .C("PrintModelList")
    stop("given model cannot be (uniquely) identified from the above list")
  }
  if (length(kappas)!=.C("GetNrParameters", covnr, k = integer(1),DUP=FALSE)$k)
    stop("incorrect number of parameters!")
  nat.scl <- double(1)
  error <- integer(1)
  .C("GetNaturalScaling",
     as.integer(covnr),
     as.double(c(0,1,0,1,kappas)),         ## not stable w.r.t. to changings !!
     as.integer(4+ length(kappas)),        ## dito
     as.integer(11),
     nat.scl,
     error,DUP=FALSE)
  if (error) stop("natural scaling could not be obtained")
  return(1.0 / nat.scl)
}


