q()
# source("RFtest.precision.R")

## file.remove("xx.dat")

if (FALSE) {
  n <- .C("getCov", integer(1), double(0))[[1]]
  C <- matrix(.C("getCov", as.integer(n), double(n*n))[[2]], ncol=n) 
}


if (FALSE) {
  z <- NULL
  a <- 1.1
  repeat {
   z <- c(z, (1-runif(100000))^{-1/a} - a/(a-1))
   print(c(length(z),mean(z)))
  }
}



if (FALSE) {
  library(RandomFields, lib="~/TMP")
  model <- "exp"
  param <- c(0,1,0,10)
  x <- seq(0,5,0.1) 
#plot(x, Variogram(x, model=model, param=param))
  z <-  GaussRF(x,x, model=model, param=param, grid=TRUE, me="TBM2", n=1000)
  unix.time(e <- EmpiricalVariogram(x,x,data=z, bin=c(0.05,seq(0.1,3,0.1)), 
					      grid=TRUE))
  plot (e$c, e$e, pch=16)
  lines(x, Variogram(x, model=model, param=param))
}
 
#runif(1);save(file="randomseed", .Random.seed)
#load("randomseed")


.path <- "/home/schlather/article/R/NEW.RF/RandomFields/R/"
if (!file.exists(paste(.path, "rf.R", sep=""))) q()

if (file.exists("source.R")) source("source.R") else
if (file.exists("../source.R")) source("../source.R")


x <- seq(0, 10, 0.1)
y <- .C("I0ML0", y=as.double(x), as.integer(length(x)))$y
print(cbind(x,y))

models <- NULL ## if NULL then all models are considered
#model <- NULL

#models <- c("exp", "sp")
#models <- "nnst"
#models <- "exp"
jump.models <- c("cone", "bessel", "cauchy",
                 "cauchytbm", "circular", #nochmal!
                 "cubic", "dampedcosine",
                 # "exponential", #: method=TBM2 T=TRUE; funktioniert nicht !!
                 "nsst", "nsst2", ## beide nicht ueberprueft!!

                 #"cauchy" ## ci-emb funktioniert nicht (gut);
                 ##           run again with grids <- TRUE

                 "2d", "3d"
                 )

models <- "spherical"

jump.models <- c("cone", "bessel", "cauchy", "cauchytbm")

#jump.models <- c("cone")
                     
jump.methods <- c("none", "circ", "local", "TBM2", "TBM3", "spectral", "direct",
                  "nugget", "add.MPP", "hyper", "other")
jump.methods <- "none"
jump.methods <- c("none", "circ", "local", "TBM3", "spectral", "direct",
                  "nugget", "add.MPP", "hyper", "other")


if (exists("time.list")) rm("time.list")
#time.list <- list(NULL)  ## is set (and then to list(NULL), no temporal
##                          components are considered
#
#
T <- runif(1,1,5); T <- c(T/2, 3*T/2, T); time.list <- list(T)
#time.list <- list(NULL)


p <- 0 # percent of randomly skipped tests
# p <- 90

dimensions <- 1:3
grids <- c(FALSE, TRUE)
#grids <- TRUE
max.unix.time <- 1 ## 0.4 ungefaehr 1 Minute
##                     bei 1 Durchlauf und repetitions=3000

max.kappas <- 20   ## maximum number of considered combination of kappa
##                     if the value is passed a random subset is taken
repetitions <- 3000 ## number of simulation the estimation is based on
rep.factor <- 10   ## if deviation then simulation is redone with increased
##                   number of simulations
#repetitions <- 10
tol = 0.05 ## tolerate deviation (in %)

#repetitions <- 15; tol <- 100



setparameters <- function(n.sf, pch="!") {
  sf <- c(5, 20) # linesimufactor
  prnt <- c(1, 6) # >6 gives too much output
  trials <- c(3, 6)
  RFparameters(TBM2.linesimufactor=sf[n.sf],
               TBM3D2.linesimufactor=sf[n.sf],
               TBMCE.force=TRUE,
               TBMCE.trials=trials[n.sf],
               TBM2.lines=120,#120,
               CE.trials=4, #trials[n.sf],
               CE.force=FALSE,
               TBM3D3.linesimufactor=sf[n.sf], pch=pch, 
               Print=prnt[n.sf], Storing=TRUE, #TBM.method="di",
               direct.method=0)
}

ENVIR = environment()
refined.simulation <- function(tol, xx, tt, grid, model, repetitions, method,
                               v, rep.factor, pch="!") {              
  DeleteAllRegisters()
  assign("z", NULL, envir=ENVIR)
  zaehler <- 0
  dev <- Inf
  setparameters(2, pch)  
  while (abs(dev) > tol) {
    zaehler <- zaehler + 1
    cat("\nlarge deviance ")
    assign("z", cbind(matrix(GaussRF(xx, T=tt, grid=grid, gridtriple=TRUE,
                              model=model, n=repetitions, method=method)
                      , ncol=repetitions), z), envir=ENVIR)
    RFparameters(Print=1)
    e <- 0.5 * mean((z[1,] - z[nrow(z), ])^2)
    dev <- if (abs(v)>1e-10) (e-v)/v else e * 100000
    cat("; simu=", e, " (", formatC(dev * 100, dig=2),
        "%)",
        #formatC((0.5*(diff(z[c(1,4), 1]))^2-v)/v * 100, dig=2),
        #"; ", ncol(z),
        #paste(formatC(z[c(1,4), 1], dig=3), collapse=", "),
        sep="")
    if (zaehler > rep.factor && abs(dev) > tol) {
      str(model)
      cat("too large deviance:", abs(e-v), ">", tol * v,"\n")
      DeleteAllRegisters()
      return(FALSE)
    }
  } # while dev > tol
  DeleteAllRegisters()
  return(TRUE)
}

allmodels <- dimnames(GetModelList())[[1]]
allmethods <- dimnames(GetModelList(FALSE))[[2]]
if (is.null(models)) models <- allmodels
models <- pmatch(models, allmodels, dup=TRUE)
methods <- GetModelList()[models, , drop=FALSE]
models <- allmodels[models]

variance <- runif(1,1,10)
scale <- runif(1,0.1,1)
T <- runif(1,1,5)
T <- c(T/2, 3*T/2, T)
if (!exists("time.list")) time.list <- list(NULL, T)

if (file.exists("xx.dat")) {
  e <- NULL
  load("xx.dat") 
#repetitions <- 1 ## number of simulation the estimation is based on
#rep.factor <- 100000   ## if deviation then simulation is redone with increased
#tol <- 0
  vv <- v
  ee <- e
  v <- Variogram(cbind(xx[2,,drop=FALSE]- xx[1,,drop=FALSE], tt[3]),
                 model=model)
  cat("\nRECHECK: \n", model$model, "grid=",grid, "; d=", d,
      "; T=", !is.null(tt),
      "; ani=", anisotropy, "; k=", paste(nk, collapse=","),
      "; vario=", v, " ", method, sep="")
  str(model); str(xx); str(tt)
  print(ncol(xx))
  if (FALSE) {
    print(as.matrix(if (grid) {
      switch(ncol(xx),
             expand.grid(xx[1:2], tt[1:2]),
             expand.grid(xx[1:2,1], xx[1:2,2], tt[1:2]),
             expand.grid(xx[1:2,1], xx[1:2,2], xx[1:2,3], tt[1:2]),
             )
    } else {
      switch(ncol(xx),
             expand.grid(xx, tt[1:2]),
             expand.grid(xx[,1], xx[,2], tt[1:2]),
             expand.grid(xx[,1], xx[,2], xx[,3], tt[1:2]),
             )
    }) %*% model[[1]]$aniso)
  }
  v <- Variogram(cbind(xx[2,,drop=FALSE]- xx[1,,drop=FALSE], tt[3]),
                 model=model)
   
  if (refined.simulation(tol, xx, tt, grid, model, repetitions, method,
                         v, rep.factor, pch="!")) {
   file.remove("xx.dat")
    stop("OK")
  }
  else {
    n <- .C("getCov", integer(1), double(0))[[1]]
    C <- .C("getCov", integer(n), double(n*n))[[2]]
    stop("recheck failed")
  }
}

for (mi in 1:length(models)) {
  m <- models[mi]
  if (!all(is.na(pmatch(jump.models, m, dup=TRUE)))) next;
  cat("\n\n\n\n ************   ", m, "var=", variance, "scale=", scale,
      "   ************")
  gridx <- runif(3, 0.1, 1)
  gridx <- rbind(gridx/2, 3*gridx/2, gridx)
  x <- list(rbind(runif(3), runif(3)), gridx)
  for (grid in grids) {
    for (d in dimensions) {
      xx <- x[[grid+1]][, 1:d, drop=FALSE]
      for (tt in time.list) {
        spt <- d + !is.null(tt)
        kappa.range <- parameter.range(m, spt)
        if (!is.null(kappa.range)) {
          if (any(is.na(unlist(kappa.range)))) next;
          kappas <- NULL
          f <- 2
          variable.kappas <- kappa.range$th[[1]][2,] != kappa.range$th[[1]][1,]
          for (k in 1:length(kappa.range$th)) {
            txt <- paste("expand.grid(",
                         paste("seq(",kappa.range$pr[[k]][1,],
                               ",", kappa.range$pr[[k]][2,],
                               ",len=", 1 + f * variable.kappas, ")", sep="",
                               collapse=", "), ")")
            kappas <- rbind(kappas, as.matrix(eval(parse(text=txt))))
          }
          dimnames(kappas) <- NULL
          if (nrow(kappas)> max.kappas)
            kappas <- kappas[ order(runif(nrow(kappas)))[1:max.kappas], ]
          kappas <- lapply(apply(kappas, 1, function(x) list(x)),
                           function(x) x[[1]])
        } else kappas <- list(NULL)
        aniso <- matrix(c(-0.5,0.5,1, 0,1.2,0.7, 0.3,0.1,0.9), ncol=3) *
          runif(1,1,10)
        ## print(kappas); yyyy
       for (anisotropy in (!is.null(tt)):1) {
          if (anisotropy) {
            ani <- aniso[1:d, 1:d]
            if (!is.null(tt))
              ani <- rbind(cbind(ani,0), runif(d + 1, -1.5, 1.5))
            model <- list(list(model=m, var=variance, aniso=ani))
          } else
            model <- list(list(model=m, var=variance, scale=scale))
          ## print(kappas); xxxx
          for (nk in kappas) {
            model[[1]]$kappa <- nk
            RFparameters(Print=2)    
            v <- Variogram(cbind(xx[2,,drop=FALSE]- xx[1,,drop=FALSE], tt[3]),
                           model=model)
            stopifnot(!is.na(v))
            cat("\n\n", m, ": grid=",grid, "; d=", d, "; T=", !is.null(tt),
                "; ani=", anisotropy, "; k=", paste(nk, collapse=","),
                "; vario=", v, sep="")
            last.failed <- FALSE
            for (me in 1:length(allmethods)) {
              if (!all(is.na(pmatch(jump.methods, allmethods[me], dup=TRUE))))
                next;
              if (methods[mi, me] && (grid || allmethods[me]!="circu")) {
                setparameters(1)
                method <- allmethods[me]
                if (runif(1,0,100) < p) {cat("%"); next}
                if (!last.failed) cat("\n")
                cat("method=", allmethods[me], sep="");
                DeleteRegister(0)
                
                save(file="xx.dat", v, xx, tt, model, grid, repetitions,
                       method, d, anisotropy, nk)
                
#            RFparameters(Print=5)    
                ut <- unix.time(z <- GaussRF(xx, T=tt, grid=grid,gridtriple=TRUE,
                                             model=model, n=10,method=method))[1]
                if (last.failed <- is.null(z)) {
                                        # cat(" failed")
                  if (FALSE)
                    {
                      if ((grid || me!=1) && (is.null(tt) || me!=5)
                          && (me!=5 || d<=2) && (me!=3 || d<=2)
                          )
                        readline(" -- what now?")
                    }
                  ## RFparameters(Print=20)
                  next
                }
                if (ut>max.unix.time) {
                  cat(" too time consuming", ut)
                  next
                } else cat(" ")
                z <- matrix(c(as.double(z),
                              as.double(GaussRF(xx, T=tt,grid=grid,gridtriple=TRUE,
                                                model=model, n=repetitions-10,
                                                method=method))),
                            ncol=repetitions)
                e <- 0.5 * mean((z[1,] - z[nrow(z), ])^2)
                dev <- if (abs(v)>1e-10) (e-v)/v else e * 100000
                cat("; simu=", e, " (", formatC(dev * 100, dig=2),"%)", sep="")
                if (abs(dev) > tol &&
                    !refined.simulation(tol, xx, tt, grid, model, repetitions,
                                        method, v, rep.factor)) {
                  save(file="xx.dat", v, xx, tt, model, grid, repetitions,
                       method, e, d, anisotropy, nk)
                  stop("check failed")
                } else file.remove("xx.dat")
              }
            } # for method
          } # for kappas
        } # for anisotropy
      } # for time
    } # for dimensions
  } # for grid yes/no
} # for models

if (FALSE) {

  .C("printKEY", as.integer(0))


}
