

# options(warn=0);library(RandomFields);source("~/R/RF/RandomFields/R/getNset.R");source("~/R/RF/RandomFields/R/rf.R");source("RandomFields/R/ShowModels.R"); ShowModels(x=1:100, y=1:100, model=list(list(model="exp", var=1, aniso=c(1,0,0,1))), method="ci");



ShowModels <- function(x, y=NULL,
                       covx=ifelse(is.null(empirical),diff(range(x))/5,
                         max(empirical$c)),
                       fixed.rs=TRUE,
                       method=NULL,
                       empirical=NULL,
                       model=NULL,
                       param=NULL,
                       all.param=NULL,## var, nugg, scale
                       legends = TRUE,
                       register=0,
                       Mean=NULL,
                       erase=TRUE,
                       x.fraction=0.60,
                       cex.names=1,
                       covx.default = 100,
                       link.fct=NULL,
                       Zlim=NULL,
                       Col.rect="red", Col.bg="blue", Col.sep="grey",
                       Col.line="red", Col.txt="black", Col.flash="red",
                       Col.vario="blue", Col.main="black",
                       Col.model=c("red", "black"), ## rf, transformed,
                       ##                                   vario
                       vario.lty=c(1,2), ## main axis, sec. axis
                       cex.leg =0.7 * cex.names,
                       cex.eval=0.8 * cex.names,
                       update=TRUE,
                       screen.new=TRUE,
                       use.outer.RFparameters=FALSE,
                       debug=FALSE,
                       ...){

  print(link.fct)
  
  stopifnot(!missing(x))
  if (any(diff(x) <= 0)) stop("x should be a sequence of increasing numbers")
  if (!is.null(y) && any(diff(y) <= 0))
    stop("y should be a sequence of increasing numbers")    
  parent.screen <- screen()
  ENVIR <- environment()
  linkfctlist <- c("MaxStable")
  warn.orig <- options()$warn
  options(warn=1)
  bg.save <- par()$bg
  par(bg="white")
  oldRFparameters <- RFparameters(no.readonly=TRUE)

  if (debug) {
    print1 <- print0 <- 5
    RFparameters(Print=print1)
  } else {
    print1 <- 1
    print0 <- 0
  }

  if (!use.outer.RFparameters)
  RFparameters(PracticalRange=FALSE, Print=print1, maxstable.maxGauss=2,
               CE.force=TRUE, CE.trials=1, CE.mmin=-4, CE.useprimes=TRUE
               )
    
  runif(1)
  assign("save.seed", get(".Random.seed", envir=.GlobalEnv, inherits = FALSE),
         envir=ENVIR)

  dim <- 2 - is.null(y)
  maxangle <- 180
  varioangle <- NA ## if NA, than the varioangle follows the angle of the model

  if (!(missing.zlim <- is.null(Zlim)) && !is.list(Zlim))
    Zlim <- list(range(Zlim), if (!is.null(link.fct)) range(link.fct(Zlim)))
  ## if (!(missing.zlim <- (sum(index <- (opts==2))!=1))) zlim <- ll[index][[1]]

  
  #######################################################################
  ## diverse
  #######################################################################
  ## set limits for the steps in parameter choice
  ## check of empirical variogram
  if (variogram <- !is.null(empirical)) {
    if (is.matrix(empirical$e) && ncol(empirical$e)==1)
      empirical$e <- is.double(empirical$e)
    if (is.matrix(empirical$e)) {
      if (dim==1)
        stop("anisotropic variogram but 1d simulations indicated by missing y")
      if (is.null(empirical$T)) {
        if (ncol(empirical$e) %% 2 !=0)
          stop("no orthogonal direction in empirical variogram (ncol is odd)")
        empirical$e <- cbind(empirical$e[ , 1 + c(0, ncol( empirical$e)/2)])
        if (!is.null(empirical$phi)) varioangle <- empirical$phi[1]
        else if (is.null(empirical$theta)) varioangle <- empirical$theta[1]
        else stop("puzzled about both empirical$phi and $theta being null")
     } else {
        maxangle <- 0
        empirical$e <- cbind(empirical$e[,1], empirical$e[1,])
      }
    } else if (is.array(empirical$e))
      stop("anisotropic variograms of higher are not allowed")    
    m.vario <- mean(empirical$e, na.rm=TRUE)
    m.lag <- diff(range(empirical$c))    
    ##                 mean        var      nugg    scale    add. par.
    maxstep <- c( c(sqrt(m.vario), m.vario, m.vario, m.lag)/5, 1)
    minscale <- min(diff(empirical$c)) / 10000
  } else {
    m.lag <- diff(range(x))/10
    maxstep <- c(1, 1, 1, m.lag, 1)
    minscale <- min(diff(x)) / 10000
  }
  if (!is.null(Mean) && !is.na(Mean) && (Mean!=0)) maxstep[1] <- abs(Mean) / 2


  ## prepare initial model if there is any
  anisotropy <- FALSE
  if (is.null(model)) {
    if (is.null(Mean)) Mean <- 0
  } else {
    Model <- PrepareModel(model=model, param=param, timespacedim=dim,
                          method=method)
    covnr <- Model$covnr + 1
    anisotropy <- ((model.aniso <- Model$anisotropy) && (dim!=1) ||
                   is.matrix(empirical$e))
    Model <- convert.to.readable(Model, allowed="list")
    if (is.null(Mean)) Mean <- Model$mean

    err.mess <- "complicated models are not allowed yet"
    if (length(Model$model)>3) stop(err.mess)
    if (length(Model$model)==3) {
      if (Model$model[[2]] != "+") stop(err.mess)
      if (Model$model[[1]]$model=="nugget") {
        swap <- Model$model[[1]]
        Model$model[[1]] <- Model$model[[3]]
        Model$model[[3]] <- swap
        if (!anisotropy && model.aniso) {
          Model$model[[1]]$scale <- abs(1 / Model$model[[1]]$aniso)
          Model$model[[1]]$aniso <- NULL
        }
        covnr <- covnr[2]
      } else {
        if (Model$model[[3]]$model!="nugget") stop(err.mess)
        covnr <- covnr[1]
      }
    } else {
      stopifnot(length(Model$model)==1)
      if (!anisotropy && model.aniso) {
        Model$model[[1]]$scale <- abs(1 / Model$model[[1]]$aniso)
        Model$model[[1]]$aniso <- NULL
      } 
      Model$model[[2]] <- "+"
      Model$model[[3]] <- list(model="nugget", var=0)
    }
    model <- Model
  } 
  if (anisotropy) {
    model$model[[3]]$scale <- NULL
    model$model[[3]]$aniso <- diag(2) # model[[3]] : nugget !!
    model$model[[1]]$aniso <-
              if (model.aniso) matrix(model$model[[1]]$aniso, ncol=dim) else
                     diag(model$model[[1]]$scale, dim)
  } else {
    model$model[[3]]$scale <- 1
    model$model[[3]]$aniso <- NULL
  }
  
  #######################################################################
  ## preparation of the models
  #######################################################################
  namen <- GetModelNames()
  n <- length(namen)

  ################### kappa -- parameter values ###########################
  cur.par <-list()
  cur.par[[n + 1]] <- NA
  if (unknown.par <- is.null(all.param)) {
    if (is.null(empirical)) all.param <- c(1, 0, diff(range(x))/5)
    else all.param <- c(Mean, 0.75 * m.vario, m.lag/5)
  } else stopifnot(length(all.param)==3) ## var, nugg, scale

  npar <- as.list(.C("GetNrParameters", as.integer(0:(n-1)), as.integer(n),
             k=integer(n), PACKAGE="RandomFields")$k)
  npar <- lapply(npar, function(x) x <- rep(1,x))
  type <- lapply(npar, function(x) rep("real", length(x)))

  ## models: models where vaslid starting values of the parameters cannont be
  ## easily be guessed
  ## (goal is that at starting point, user always gets a set of valid
  #   parameter values)
  m.par.except <- list("cone", "cauchy", "cauchytbm", "power", "gengneiting",
                       "lgd1", "FD", )
  npar.except <- list(c(0.5,0.5,0.5), 2, c(1,1,2), 2, c(1,3),
                      c(0.45, 1), 0)

  ## models with discrete parameters: 
  m.type.except <- list("gengneiting", "nsst", "nsst2")
  type.except <- ## which parameters take only discrete values? 
    list(c("discrete","real"),
         c("real", "discrete", "real", "real", "discrete", "real"),
         c("real", "real", "discrete", "real", "real", "discrete", "real"),
         )
  
  scale <- rep(all.param[3], n)
  m.scale.except <- list("bessel", "wave", "gengneiting")
  scale.except <-  c(1/3, 1/3, 1/3)

  ## user readable function names turned into covariance numbers
  ## for the exceptions
  for (i in 1:length(m.par.except))
    npar[[.C("GetModelNr", m.par.except[[i]], as.integer(1),
            nr=integer(1), PACKAGE="RandomFields")$nr + 1]] <-  npar.except[[i]]
  for (i in 1:length(m.type.except))
    type[[.C("GetModelNr", m.type.except[[i]], as.integer(1),
            nr=integer(1), PACKAGE="RandomFields")$nr + 1]] <-  type.except[[i]]
  for (i in 1:length(m.scale.except)) {
    j <- .C("GetModelNr", m.scale.except[[i]], as.integer(1),
            nr=integer(1), PACKAGE="RandomFields")$nr + 1
    scale[j] <- scale[j] * scale.except[i]
  }

  #######################################################################
  ## auxiliary functions
  #######################################################################
  ## interactive menue to choose the name of a class of variogram models  
  col.choose <- rep(Col.txt, n)  
  choose.model <- function(nr=NA) {
    screen(name.dev)
    par(mar=name.mar)
    col.choose[nr] <- Col.flash
    plot(-Inf, -Inf, xlim=c(0,ncol), ylim=c(1,top+1), axes=FALSE,
         xlab="", ylab="")
    title(main=maintitle, col.main="black")
    text(as.integer(((1:n)-1) / maxrow), top - ((1:n)-1) %% maxrow,
         labels=paste(namen), adj=c(0,0), cex=cex.names,
         col=col.choose)
    repeat {
      if (length(loc <- Locator(1))==0) {
        return(NA)
      }
      loc <- floor(unlist(loc))
      covnr <- (1+top-loc[2]) + maxrow * loc[1]
      if ((covnr>0) && (covnr<=n) && (loc[2]>=0) && (loc[2]<=top)) {
        col.choose[nr] <- Col.txt
        return(covnr)
      }
    }    
  }

  
  ## used in menue plot to define non-linear changings by user
  quadratic <- function(d, v, a, mini=0, maxi=Inf) {
    d <- pmin(1, pmax(0, d)) - 0.5
    d <- ((d>0) * 2 - 1) * d^2 * a * 4
    if (missing(v)) d else pmax(mini, pmin(maxi, v + d))
  }

  
  ## simulates Gaussian random field according to the users choice of
  ## parameters (cp) and calculates the variogram
  simulate <- function(cp, param=c("vario", "field", "simu", "rs")) {
    ## "rs" only active if fixed.rs, and new random seed required
    ## "vario" : refresh of variogram plot
    ## "field" : refresh of simulation plot
    ## "simu" : new simulation and recalculation of variogram model
   
    if (cp$PracticalRange != RFparameters()$PracticalRange) {
      DeleteRegister(register)
      RFparameters(PracticalRange=cp$PracticalRange)                   
    }
    
    if (any(param %in% c("field", "simu"))) {
      screen(simu.dev, new=screen.new)
      if (!is.null(y)) {
        plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
        text(0.5, 0.5, lab="calculating...", adj=c(0.5, 0.5), cex=2)
      }
    }

    if (any("vario" == param) || !is.null(link.fct)) { 
      screen(model.dev)
      par(mar=model.mar)
      if (anisotropy) {
        f <- pi / 180 * cp$angle
        u <- matrix(c(cos(f), sin(f), -sin(f), cos(f)), ncol=2)
        cp$model[[1]]$aniso <- u %*% (1/cp$diag * t(u))
        if (is.finite(cp$varioangle)) {
          f <- pi / 180 * cp$varioangle
          uu <- matrix(c(cos(f), -sin(f), +sin(f), cos(f)), ncol=2)
        } else uu <- t(u)
        covxx <-
          ## works ! try out by:
          ## matrix(outer(c(10,20,30), matrix(1:4, ncol=2),  "+"), ncol=2)
          matrix(outer(covx, uu, "*"), ncol=2) # $v : dir in each row
      } else {
        covxx <- covx
      }

      RFparameters(Print=print0)
      covvalue <-
        eval(parse(text=paste(c("CovarianceFct", "Variogram")[cp$variogram + 1],
                     "(covxx, model=cp, dim=dim)")))
      link.cov <- NULL
      linkname <- "transf."
      if (!is.na(covvalue[1]) && !is.null(link.fct)) {
        if (cp$variogram) {
          corvalue <- CovarianceFct(covxx, model=cp, dim=dim)
          if (is.na(corvalue[1])) {
            linkname <- "* trafo"
            corvalue <- Variogram(covxx, model=cp, dim=dim)
            corvalue <- max(corvalue) - corvalue
          }
        } else corvalue <- covvalue
        if (!is.na(sd <- corvalue[1])) {
          corvalue <- if (sd > 0) corvalue / sd else rep(1, length(corvalue))
          sd <- sqrt(sd)
          xx <- rnorm(10000, 0, sd)
          ## create 10000 bivariate norms with mean mn, sd, and corvalue
          ## and calculate covariance empirically for link.fct
          if (is.function(link.fct)) {
            link.cov <-
              apply(link.fct(cp$mean + outer(xx, corvalue, "*") +
                             outer(rnorm(10000, 0, sd), sqrt(1-corvalue^2),"*")),
                    2, cov, y=link.fct(cp$mean + xx), use="complete.obs")
            if (cp$variogram) link.cov <- link.cov[1] - link.cov
          } else {
            switch(pmatch(link.fct, linkfctlist),
                   {
                     link.cov <-
                        if (cp$mean==0) sqrt((1 - corvalue)/2)
                        else rep(NA, length(corvalue))
                     linkname <- "theta-1"
                   },
                   stop("unknown link function")
                   )
          }
          link.cov <- matrix(link.cov, ncol=1 + anisotropy)
        }         
      } # !is.na(covvalue[1]) && !is.null(link.fct)
      if (!any(z.f <- is.finite(z <- c(covvalue, if (cp$variogram) empirical$e,
                                       link.cov)))){
        plot(0,0,col=0,axes=FALSE)
        text(0,0,label="plot not available", adj=c(0.5,0.5), col=Col.flash)
      } else {      
        Ylim <- range(z[z.f])
        covvalue <- matrix(covvalue, ncol=1 + anisotropy)
        plot(covx[1], covvalue[1,1], type="p", xlim=range(covx), ylim=Ylim,
             col=Col.model)
        if (!is.null(link.cov)) points(covx[1], link.cov[1,1], col=Col.model[2])
        if (any(is.finite(covvalue[-1, ])))
          matplot(covx[-1], covvalue[-1, , drop=FALSE],
                  col=Col.model[1], lty=vario.lty, type="l", add=TRUE)
        if (!is.null(link.cov) && any(is.finite(link.cov[-1, ])))
          matplot(covx[-1], link.cov[-1, , drop=FALSE],
                  col=Col.model[2], lty=vario.lty, type="l", add=TRUE)
        if (!is.null(empirical) && (cp$variogram)) {
          vario.pch <- if (is.matrix(empirical$e)) c(16, 1) else 16
          matplot(empirical$c, empirical$e, col=Col.vario, pch=vario.pch,
                  add=TRUE)
        } else vario.pch <- -1
        if (!is.null(link.fct) || anisotropy) {
          if(is.null(link.fct)) {
            leg <- c("1st", "2nd")
            llty <- vario.lty
            lcol <- Col.model[1]
          } else {
            lcol <- Col.model
            leg <- c("Gauss", linkname)
            llty <- rep(vario.lty, each= 1 + anisotropy)
            if (anisotropy) 
              leg <- as.vector(outer(leg, c("1st", "2nd"), paste, sep=", "))
          }
          legend(x=covx[length(covx)], y=Ylim[2 - cp$variogram], leg=leg,
                 pch=vario.pch,
                 lty=llty, col=lcol, xj=1, yj=1 - cp$variogram, cex=cex.leg)
        }
      }
    }

    RFparameters(Print=print1)
    
    if (any(param %in% c("simu", "field"))) {
      if ("simu" %in% param) {
        ## always calculated, even if not used in the plot
        ## (in case of max-stable "link"-function)
        ## in that case it does not cost too much more ...
        if (!fixed.rs || ("rs" %in% param))
          assign("save.seed", get(".Random.seed", envir=.GlobalEnv,
                                  inherits = FALSE), envir=ENVIR)
        else assign(".Random.seed", save.seed, envir=.GlobalEnv)

        z <- GaussRF(x, y, model=cp, grid=TRUE, register=register)
        
        assign("zsimu", z, envir=ENVIR)
        assign("ztrafo", NULL, envir=ENVIR)
      } else z <- zsimu
      if (cp$link) { ## if TRUE then also !is.null(link.fct)
        if (is.null(ztrafo)) { ## either simuparam, or user has swapped
          ##                      presentation
          if (is.function(link.fct)) z <- link.fct(z) else
          switch(pmatch(link.fct, linkfctlist),
                 {
                   assign(".Random.seed", save.seed, envir=.GlobalEnv)
                   z <- MaxStableRF(x, y, model=cp, maxstable='extremalGauss',
                                    grid=TRUE, register=register)
                 },
                 stop("unknown link function")
                 )
          assign("ztrafo", z, envir=ENVIR)
        } else z <- ztrafo
      }
      screen(simu.dev)
      par(mar=simu.mar)
      if (is.null(z) || all(is.na(z))) {
        plot(0, 0, col=0, axes=FALSE, xlab="", ylab="")
        text(0, 0, label="image not available", adj=c(0.5, 0.5))
      } else {
        if (any(is.na(z))) warning("some z-values are NA")
        if (is.null(y)) {
          plot(x, z, ylim=Zlim[[cp$link + 1]], ...)
          if (legends) {
            text(txt.x, max(z, na.rm=TRUE), adj=c(0.5,1), cex=1.5,
                 if (cp$link) if (is.character(link.fct)) link.fct
                 else "transformed" else "Gaussian", col=Col.main);
          }
        } else {
          zlim <- if (missing.zlim) range(z) else Zlim[[cp$link + 1]]
          image(x, y,  z, zlim=zlim, ...)
          if (legends) {
            if (big.legend) ml <- c(format(mean(zlim), dig=2),filler)
            else ml <- NULL
            text(txt.x, txt.y, adj=c(0.5,1), cex=1.5,
                 if (cp$link) if (is.character(link.fct)) link.fct
                 else "transformed" else "Gaussian", col=Col.main);
            legend(lu.x, lu.y, y.i=y.i, x.i=0.1, yj=0,
                   legend=c(format(zlim[2],dig=2),filler,
                     ml, format(zlim[1],dig=2)),
                   lty=1, lwd=lwd, col=col[length(col):1], cex=cex.leg)
          }
        }
      } # is.null(z)
    } # c("simu", "field") 
    return(cp) ## return required by function eval.param
  }


 
  #######################################################################
  ## recognized formulae for covariance functions
  #######################################################################
  covlist <- c("bessel","cauchy","cauchytbm","circular",
               "cone","cubic","dampedcosine","exponential",
               "FD", "fractgauss", "gauss",
               "gencauchy","gengneiting","gneiting",      
               "hyperbolic", "lgd1", "nsst", "nsst2", "nugget",
               "penta","power",
               "qexponential","spherical","stable","wave",
               "whittlematern", "fractalB"
               )
  see.manual <- "formula see help page of CovarianceFct"
  exprlist <- c(expression(2^a *Gamma(a+1)*x^{-a}* J[a](x)),
                expression((1+x^2)^{-a}),
                expression((1+(1-b/c)*x^a)*(1+x^a)^{-b/a-1}),
                expression((1-2/pi*(x *sqrt(1-x^2)+asin(x))) * (x<1)),
                # "cone",...
                expression("C(x) not available"),
                expression((1-7*x^2 + 8.75 *x^3 - 3.5* x^4+0.75 *x^6)*(x<1)),
                expression(e^{-a* x}* cos(x)),
                expression(e^{-x}),
                # "FD",...
                expression((-1)^x * Gamma(1-a)^2/(Gamma(1-a+x) * Gamma(1-a-x))
                    * " x=integer"),
                expression(0.5 * (abs(x+1)^a - 2 * abs(x)^a + abs(x-1)^a)),
                expression(exp(-x^2)),
                # "gencauchy",...
                expression((1+x^a)^{-b/a}),
                expression(see.manual),
                expression((1 + 8 *s *x + 25* s^2* x^2 + 32*
                    s^3 *x^3)*(1-s* x)^8 * (x<1)),
                # "hyperbolic",...
                expression(const[a][b][c] * (c^2 +x^2)^{b/2} *
                    K[b](sqrt(a(c^2 + x^2)))),
                expression({1- b * (a+b)^{-1} * x^a} * (x <= 1) +
                    a * (a+b)^{-1} *x ^b  * (x>1)),
                expression(see.manual),
                expression(see.manual),
                expression((x==0)),
                # "penta",...
                expression(( 1 - 22*x^2/3  +33 *x^4  - 77*x^5/2  + 33* x^7/2
                    - 11* x^9/2 + 5 *x^11 /6)*(x<1)),
                expression((1-x)^a  * (x<1)),
                # "qexponential",...
                expression((2*e^{-x}-a*e^{-2*x})/(2-a)),
                expression((1 - 1.5* x + 0.5* x^3)*(x<1)),
                expression(exp(-x^a)),
                expression(sin(x)/x),
                # "whittlematern",...
                expression(2^{1-a}* Gamma(a)^{-1}* x^a * K[a](x)),
                expression(x^a)
                ) 
  expr <- rep(expression("C(x) unknown"), n)
  idx <- pmatch(covlist, namen)
  if (any(is.na(idx))) stop("programming error -- inform maintainer")
  expr[idx] <- exprlist
  DeleteRegister(register)

 
  #######################################################################
  ## graphical preparations
  #######################################################################

  ## check x coordinates for the plot of variogram models
  stopifnot(all(covx>=0))
  if (length(covx)>1) {
    covx <- covx[covx>0]
    covx <- c(0, sort(c(max(covx)/100000, covx)))
  } else {
    covx <- c(0,seq(covx/100000, covx, l=covx.default))
  }
  covx <- covx[is.finite(covx)]


  ## preparation of legend for 2-dim image plot
  if (dim==2) {
    ll <-  list(...)
    opts <- lapply(names(ll), pmatch, c("col" # ,"zlim"
                                        ), no=0)
    if (length(opts)==0) opts <- 0 else opts <- as.integer(opts)

    if (sum(index <- (opts==1))==1) col <- ll[index][[1]]
    else {
      eval(parse(text = paste("col <-",
                   paste(as.character(as.list(args(image.default))$col),
                         collapse="("),")")))
    }
    cn <- length(col)
    if (cn==1) stop ("more than one colour must be given!")
    lwd <- 1
    y.i <-  0.03
    if (big.legend <- (cn>50)) {
      filler <- vector("character",length=(cn-3)/2)
      if (cn>80) y.i <-  2.4 / cn
    } else {
      filler <- vector("character",length=cn-2)
      if (cn==2) {
        lwd <- 5
        y.i <- 1
      }
      else if (cn<30) {
        y.i <- 0.9 / (cn-2)
        lwd <- 1 + 10/(cn-2)
      }
    }
    lu.x <- min(x)
    lu.y <- min(y)
  }
  txt.x <- sum(range(x))/2
  if (!is.null(y)) txt.y <- max(y)


  ## preparation of screens
  open.screen <- function() {
    if (is.numeric(parent.screen)) screen(parent.screen)
    assign("scr", split.screen(figs=rbind(
                                 c(0.01,x.fraction,0.01,0.49),
                                 c(x.fraction+0.02,0.99,0.01,1),
                                 c(0.01,x.fraction,0.51,0.99),
                                 c(x.fraction+0.02,0.99,0.01,0.99)),
                               erase=erase), envir=ENVIR)
    assign("name.dev", scr[4], envir=ENVIR)
    assign("simu.dev", scr[1], envir=ENVIR)
    assign("model.dev", scr[3], envir=ENVIR)
    assign("par.dev", scr[2], envir=ENVIR)
    assign("model.mar", c(2,2,0,0), envir=ENVIR)
    assign("name.mar", c(0,0,2,0), envir=ENVIR)
    assign("simu.mar", c(2,2,0,0), envir=ENVIR)
  }

  open.screen()
  on.exit({
    close.screen(scr);
    options(warn=warn.orig);
    RFparameters(oldRFparameters); par(bg=bg.save)
  })
  

  screen(name.dev)
  par(mar=name.mar)
  maxrow <- 20
  ncol <- 1 + as.integer( (n-1) / maxrow)
  top <- min(maxrow+1,n,na.rm=TRUE)
  maintitle <- "Models"

  ## initialization of model, choosing first model  
  if (is.null(model$model[[1]])) {
    if (!is.null(empirical)) {
      screen(model.dev)
      par(mar=model.mar)
      Ylim <- range(empirical$e, na.rm=TRUE)
      matplot(empirical$c, empirical$e,xlim=range(covx, na.rm=TRUE),
           col=Col.vario, pch=c(16,1))       
    }
    if (is.na(covnr <- choose.model())) return(NULL)
    cur.par[[covnr]] <-
      list(model=list(model=list(model=namen[covnr], var=all.param[2],
                        kappas=npar[[covnr]]),
             "+",
             model=list(model="nugget", var=all.param[2])))  
    if (anisotropy) {
      cur.par[[covnr]]$model[[1]]$aniso <- diag(1/scale[covnr], 2)
      cur.par[[covnr]]$model[[3]]$aniso <- diag(2)
    } else {
      cur.par[[covnr]]$model[[1]]$scale <- scale[covnr]
      cur.par[[covnr]]$model[[3]]$scale <- 1
    }
  } else {
    cur.par[[covnr]] <- model

    if (anisotropy) {
      SVD <- svd(cur.par[[covnr]]$model[[1]]$aniso)
      u <- SVD$u[,1]
      if (u[1] < 0) u <- -u
      cur.par[[covnr]]$diag <- 1 / SVD$d
      if (is.na(varioangle)) { # no empirical variogram with 2 directions 
        cur.par[[covnr]]$angle <- acos(u[1]) / pi * 180
        if (u[2]<0) cur.par[[covnr]]$angle <- 180 - cur.par[[covnr]]$angle
      } else  cur.par[[covnr]]$angle <- (varioangle + 360) %% 180
    }
  } # else, is.null(model$model[[1]] 

  
  cur.par[[covnr]]$varioangle <- NA  ## default: varioangle follows angle
  cur.par[[covnr]]$PracticalRange <- RFparameters()$PracticalRange
  cur.par[[covnr]]$variogram <-  variogram
  cur.par[[covnr]]$mean <- Mean
  cur.par[[covnr]]$link <- FALSE
  oldnr <- covnr

  ## preparation of menues
  if (anisotropy) {
    model.entry <-
      list(
           list(name="Anisotropy Parameters", var=NULL, val="simulate",
                param=c("simu", "rs")),       
           list(name="first axis scale", var="diag[1]", delta=TRUE,
                val=function(d, v) {quadratic(d=d, v=v, a=maxstep[4],
                  mini=minscale) }),
           list(name="second axis scale", var="diag[2]", delta=TRUE,
                val=function(d, v){quadratic(d=d, v=v, a=maxstep[4],
                  mini=minscale) }),
           list(name="angle (degrees)", var="angle",
                delta=FALSE, val=function(d, v) pmin(maxangle,
                               pmax(0,d * maxangle))),
           )
  } else model.entry <- NULL

  model.entry <-
    c(list(
           if (!anisotropy)
           list(name="scale", var="model[[1]]$scale", delta=TRUE,
                val=function(d, v) {quadratic(d=d, v=v,
                  a=maxstep[4],mini=minscale)}),
           list(name="variance", var="model[[1]]$var", delta=TRUE,
                val=function(d, v) {quadratic(d=d, v=v, a=maxstep[2]) }),
           list(name="nugget", var="model[[3]]$var", delta=TRUE,
                val=function(d, v) {quadratic(d=d, v=v, a=maxstep[3]) }),
           ),
      model.entry,
      list(
           list(name="Global Parameters", var=NULL, val="simulate",
                param=c("simu", "rs")),       
           list(name="mean", var="mean", delta=TRUE,
                val=function(d, v) {quadratic(d=d, v=v, a=maxstep[1],
                  mini=-Inf)}, param="simu"),
           list(name="practical range", var="PracticalRange", val=TRUE),
           list(name="variogram", var="variogram", val=TRUE, param="vario"),
           if (anisotropy) list(name="variogram angle (degrees)",
                                 var="varioangle",
                                 delta=FALSE,
                                 val=function(d, v) pmin(180, pmax(0,d * 180)),
                                 param="vario"
                                 ),
           if (!is.null(link.fct)) list(name="show transformed field",
                                        var="link", val=TRUE, param="field"),
           ),
      )
 
  #######################################################################
  ##  the graphical interface itself
  #######################################################################
  repeat {
    ## is the chosen model class chosen for the first time?
    ## if so copy "standard" values from former choices
    if (is.null(cur.par[[covnr]])) {
      cur.par[[covnr]] <- cur.par[[oldnr]] ## var, nugget, (scale/aniso,) mean
      cur.par[[covnr]]$model[[1]]$model <- namen[covnr]     
      ## event fall unterscheidung npar[[covr]] == / != 0
      cur.par[[covnr]]$model[[1]]$kappas <- npar[[covnr]]
      if (!cur.par[[covnr]]$PracticalRange) {
        if (anisotropy) {
          cur.par[[covnr]]$model[[1]]$aniso <- diag(1/scale[covnr], 2)
          cur.par[[covnr]]$angle <- 0
          cur.par[[covnr]]$diag <- rep(scale[covnr], 2)
        }
        else cur.par[[covnr]]$model[[1]]$scale <- scale[covnr]
      }
    } else {
      cur.par[[covnr]]$PracticalRange <- cur.par[[oldnr]]$PracticalRange
      cur.par[[covnr]]$variogram <- cur.par[[oldnr]]$variogram
      cur.par[[covnr]]$mean <- cur.par[[oldnr]]$mean
    }

    ## define menue points that depend on the choice of the model
    kappa.entry <- NULL
    if ((l <- length(cur.par[[covnr]]$model[[1]]$kappas)) > 0) { 
      for (i in 1:l) {
        rd <- pmatch(type[[covnr]][i], c("real","discrete"))          
        kappa.entry[[l + 1 -i]] <-
          list(name=letters[i], var=paste("model[[1]]$kappas[",i,"]"),
               delta=TRUE,                                
               val=switch(rd,
                 function(d, v) {quadratic(d=d, v=v, a=maxstep[5], mini=-Inf)},
                 function(d, v) {if (missing(v)) v<-rep(0, length(d));
                                 v[d>0.5] <- v[d>0.5] + 1;
                                 v[d<0.5] <- v[d<0.5] - 1;
                                 v},
                 ),
               col=if (rd==1) "blue" else "lightblue"
               )
        }
    }
    
    entry <- c(list(
                    if (!update) list(name="simulate", val="simulate",
                                      param=c("simu", "rs")),
                    list(name=cur.par[[covnr]]$model[[1]]$model, var=NULL,
                         col=Col.flash, cex=0.8, val="simulate",
                         param=c("simu", "rs")),
                    list(name=expr[covnr], var=NULL, col=Col.txt, cex=0.8,
                         val="simulate", param=c("simu", "rs"))
                    ),
               kappa.entry,
               model.entry)
    ## the bad side of having simpler expression for the entries, is
    ## that NULLs might have been introduced...
    entry <- entry[!sapply(entry, is.null)]
    
    if (!is.null(y)) {
      screen(simu.dev)
      par(mar=simu.mar)
      plot(Inf, Inf, xlim=range(x), ylim=range(y), axes=FALSE)
    }
    simulate(cp=cur.par[[covnr]])

    options(warn=-1)
    cur.par[[covnr]] <- ## menu call
      eval.parameters("cp", entry,
                      update=update,  simulate=simulate,
                      dev = par.dev, cp=cur.par[[covnr]],
                      col.rect=Col.rect, col.bg=Col.bg, col.sep=Col.sep,
                      col.line = Col.line, col.txt=Col.txt, sep=NULL,
                      cex=cex.eval, cex.i=cex.eval,
                      param=c("field", "simu", "vario")
                      )
 
    oldnr <- covnr
    if (is.na(covnr <- choose.model(oldnr))) break
    close.screen(scr)
    open.screen()
  }
  if (anisotropy && !update) {
    f <- cur.par[[oldnr]]$angle * pi / 180
    u <- matrix(c(cos(f), sin(f), -sin(f), cos(f)), ncol=2)
    cur.par[[oldnr]]$model[[1]]$aniso <- u %*% (1/cur.par[[oldnr]]$diag * t(u))
  }
  return(convert.to.readable(PrepareModel(cur.par[[oldnr]], timespace=dim)))
}

 
