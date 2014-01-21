## Methods for classes 'RFempVario' and 'RFfit'  #######################

## accessing slots
accessByNameOrNumber <- function(x, i, j, drop=FALSE) {
  stopifnot(length(i)==1)
  if (is.numeric(i))
    i <- slotNames(x)[i]
  return(accessSlotsByName(x=x, i=i, j=j, drop=drop))
}

setMethod("[", signature = "RFfit", def=accessByNameOrNumber)


## print and show
print.RFfit <- function(x, ..., max.level = 3) {
  cat("*** object of Class 'RFfit' ***\n  ** str(object):\n  ")#
  str(x, max.level = max.level)#
  cat("*** end of 'RFfit' ***\n")
} 

setMethod(f="show", signature="RFfit",
          definition=function(object) print(object))


## coercion methods
setAs("RFempVariog", "list",
      def=function(from) {
        li <- list()
        for (name in slotNames(from)) {
          li[[name]] <- eval(parse(text=paste("from@", name, sep="")))
        }
        return(li)
      })

setAs("RFfit", "RFempVariog", def=function(from) from@ev)

list2RFempVariog <- function(li)
  return(new("RFempVariog",
             centers=li$centers,
             emp.vario=li$emp.vario,
             var=var,
             sd=li$sd,
             n.bin=li$n.bin,
             phi.centers=li$phi.centers,
             theta.centers=li$theta.centers,
             T=li$Tbins,
             coord.units = RFoptions()$general$coord_units,
             variab.units = RFoptions()$general$variab_units,
             call=li$call ))



## plot method

setMethod(f="plot", signature(x="RFempVariog", y="missing"),
	  definition=function(
            x, y, nmax.phi=6, nmax.theta=3, nmax.T=3,
            plot.nbin=TRUE, plot.sd=FALSE, model = NULL, variogram=TRUE, ...)
          
          plotRFempVariog(x, nmax.phi=nmax.phi, nmax.theta=nmax.theta,
                          nmax.T=nmax.T, plot.nbin=plot.nbin, plot.sd=plot.sd,
                          model = model, variogram=variogram,
                          ...))

setMethod(f="plot", signature(x="RFfit", y="missing"),
	  definition=function(
            x, y, fit.method="ml", nmax.phi=6, nmax.theta=3, nmax.T=3,
            plot.nbin=TRUE, plot.sd=FALSE, model = NULL, variogram=TRUE, ...)
          
          plotRFempVariog(x, fit.method=fit.method,
                          nmax.phi=nmax.phi, nmax.theta=nmax.theta,
                          nmax.T=nmax.T, plot.nbin=plot.nbin, plot.sd=plot.sd,
                          model = model, variogram=variogram,
                          ...))



plotRFempVariog <- function(x, y, nmax.phi, nmax.theta, nmax.T,
                            plot.nbin, plot.sd, fit.method = NULL, model,
                            variogram=variogram,
                            coord.units=c(""),
                            variab.units=c(""),
                            ...) {
  
  if (!variogram)
    stop("plot of estimated covariance functions not programmed yet.")
  
  dots = list(...)

  dotnames <- names(dots)
  if (!("type" %in% dotnames)) dots$type <- "b"
  
  if (!is(x, "RFempVariog")&&!is(x, "RFfit"))
    stop("method only for objects of class 'RFempVariog' or 'RFfit'")

  model <- c(model)
  if(!is.null(model))
    if (!all(unlist(lapply(model, FUN=function(x) is(x, ZF_MODEL)))))
      stop("model must be (a list of elements) of class ZF_MODEL")
  
  RFfit <- is(x, "RFfit")
  newx <- list()
  methodnames <- double()
  
  if (RFfit){    
    newx$autostart <- x@autostart
    newx$self <- x@self
    newx$plain <- x@plain
    newx$sqrt.nr <- x@sqrt.nr
    newx$sd.inv <- x@sd.inv
    newx$internal1 <- x@internal1
    newx$internal2 <- x@internal2
    newx$internal3 <- x@internal3
    newx$ml <- x@ml
    if(length(x@ev)!=0)
      x <- x@ev[[1]]
    else
      stop("The fit does not contain an empirical variogram.")
    if( is.null(fit.method) )
      fit.method <- if (length(newx$ml@name) > 0) "ml" else "plain"

    methodidx <- fit.method %in% names(newx[unlist(lapply(newx, FUN = function(x) is(x, ZF_MODELEXT)))])
    methodnames <- fit.method[methodidx]
    nomethodnames <- fit.method[!methodidx]

    if( !all(methodidx) )
      warning( paste("The following method does not exist: ", nomethodnames) )
  } 

  variab.name <-
    if (is.matrix(x@emp.vario)) dimnames(x@emp.vario)[[2]][1]
    else names(x@emp.vario)[1]

  if(!is.null(model)){
    len.modelnames <- length(names(model))
    modelnames <- if(len.modelnames > 0) names(model) else paste("model", 1:length(model))
    methodnames <- c(methodnames, modelnames)
    names(model) <- modelnames
    newx <- c(newx, model)
  }
  
  n.methods <- length(methodnames)
  

  ## case without binning
  if (is.list(x@centers)) {

    coords <- GridTopology2gridVectors(cbind(x@centers$x, x@centers$T))
    if (length(coords)>1) {
      coords[[1]] <- sort(unique((coords[[1]] - min(coords[[1]])) *
                                 rep(c(-1, 1), each=length(coords[[1]]))))
      lab.names <- dimnames(x@centers$x)[[2]]
      
      if (!(x@centers$spacedim==1 && x@centers$Time) &&
          !(x@centers$x[3,2]==dim(x@emp.vario)[2]))
        coords[[2]] <- sort(unique((coords[[2]] - min(coords[[2]])) *
                                   rep(c(-1, 1), each=length(coords[[2]]))))
      else
        lab.names[2] <- "T"
    } else {
      x@centers <- coords[[1]]-min(coords[[1]])
      do.call(plot, args=c(dots, list(x=x, plot.nbin=FALSE)))
      return()
    }
   
    if (!("main" %in% dotnames)) {
      main <- "Variogram image plot"
      if (length(variab.name)>0) main <- paste(main, "for", variab.name)
      dots$main <- main
    }
    lab.names <- paste(lab.names, "-distance", sep="")

    idx <- lab.names != "T-distance"
    if (any(idx) && all(coord.units[idx] != ""))
      lab.names[idx] <-
        paste(lab.names[idx], " [", coord.units[idx], "]", sep="")
    if (!all(idx) && all(coord.units[!idx] != ""))
      lab.names[!idx] <-
        paste(lab.names[!idx], " [", coord.units[!idx], "]", sep="")
    
    if (!("xlab" %in% dotnames)) dots$xlab <- lab.names[1]
    if (!("ylab" %in% dotnames)) dots$ylab <- lab.names[2]
    if (!("xlim" %in% dotnames)) dots$xlim <- range(coords[[1]])
    if (!("ylim" %in% dotnames)) dots$ylim <- range(coords[[2]])

    idx1 <- coords[[1]] >= dots$xlim[1] & coords[[1]] <= dots$xlim[2]
    idx2 <- coords[[2]] >= dots$ylim[1] & coords[[2]] <= dots$ylim[2]
    coords[[1]] <- coords[[1]][idx1]
    coords[[2]] <- coords[[2]][idx2]

    dims <- dim(x@emp.vario)
    ev.plot <- matrix(x@emp.vario, nrow=dims[1], ncol=dims[2])
    ev.plot <- ev.plot[idx1, idx2]
    
    range.ev <- range(ev.plot)#x@emp.vario)
    col <- if ("col" %in% dotnames) dots$col else
              default.image.par(NULL, NULL)$default.col

    split.screen( rbind(c(0,.85,0,1), c(.85,1,0,1)))
    screen(2)
    par(mar=c(4,0,3,3))
    z.legend <- seq(range.ev[1], range.ev[2], length=length(col))
    image(y=z.legend, x=c(0,1), z=rbind(z.legend, z.legend), axes=F, col=col,
          xlab="")
    axis(4)
    box()
    screen(1)
    par(mar=c(4,4,3,1))
    do.call(image, args=c(dots[-which(names(dots)=="type")],
                     list(x=coords[[1]], y=coords[[2]], z=ev.plot)))
    #      xlab=lab.names[1],
    #      ylab=lab.names[2],
    #      main="Variogram image plot")
    close.screen(all.screens=TRUE)
    return(invisible())
  }

  ## case with binning
  
  #if (is.null(x@phi.centers)) x@phi.centers <- 0
  #if (is.null(x@theta.centers)) x@theta.centers <- 0
  #if (is.null(x@T)) x@T <- 0

#  if (is.null(dim(x@emp.vario))) x@emp.vario <- as.marix(x@emp.vario)
#  if (!is.null(x@sd))
#    if (is.null(dim(x@sd))) x@sd <- as.marix(x@sd)
#  if (is.null(dim(x@n.bin))) x@n.bin <- as.marix(x@n.bin)

  n.phi <- min(nmax.phi, l.phi <- max(1,length(x@phi.centers)))
  n.theta <- min(nmax.theta, l.theta <- max(1,length(x@theta.centers)))
  n.T <- min(nmax.T, l.T <- max(1,length(x@T)))  
  halfphi.sector <- pi/(2*l.phi)
  halftheta.sector <- pi/(2*l.theta)
  
  dim(x@emp.vario) <- newdim <- c(length(x@centers), l.phi, l.theta, l.T)
  

  if (!is.null(x@n.bin)) dim(x@n.bin) <- newdim
  if (has.sd <- !is.null(x@sd))  dim(x@sd) <- newdim

  range.nbin <- range(c(0, x@n.bin), na.rm=TRUE)
  ylim.nbin <- range.nbin * c(1,1.2)
  
  range.ev <- range(x@emp.vario, na.rm=TRUE)
  if (!("ylim" %in% dotnames)) dots$ylim <- range.ev

  if(!("xlim" %in% dotnames)) dots$xlim <- range(x@centers)
  
  col.v <- col <- 
    if ("col" %in% dotnames) rep(dots$col, len=n.phi) else 1:max(n.phi)
  dots$col <- NULL

  cex <- if ("cex" %in% dotnames) dots$cex else .8
  par(cex=cex)
  dots$cex <- NULL

  if (!("pch" %in% dotnames)) dots$pch <- 19

  if ("main" %in% dotnames)  main <- dots$main
  else {
    main <- "Variogram plot"
    if (length(variab.name)>0) main <- paste(main, "for", variab.name)
  }
  dots$main <- NULL
  if (!is.null(main)) oma.top <- 2 else oma.top <- 0

  ylab.ev <- if ("ylab" %in% dotnames) dots$ylab else "semivariance"
  dots$ylab <- NULL
  
  xlab <- if ("xlab" %in% dotnames) dots$xlab else "distance"
  dots$xlab <- NULL
  
  split.screen(c(n.T, n.theta))
  n.all <- n.theta * n.T

  oma.left <- 6

  if (n.methods!=0){
    dotsRFfit <- dots
    dotsRFfit$type <- "l"
    dotsRFfit$lwd <- 2
    ltyRFfit.v <- 1:n.methods
    dotsRFfit$lty <- NULL
  }
  
    
  screen <- 1
  for (iT in 1:n.T)
    for (ith in 1:n.theta) {
      ## plot n.bin
      if (plot.nbin) {
        screen(screen)
        par(oma=c(4,oma.left,oma.top,0))
        split.screen(rbind(c(0,1,.2,1), c(0,1,0,.2)), screen=screen)
        screen(n.all + 2*screen)
        par(mar=c(0,.5,0,.5))
        for (iph in 1:n.phi) {
          if (n.phi > 1) col <- col.v[iph]
          if (iph==1) {
            lab <- xylabs("bin centers", NULL, units=x@coord.units)
            plot(x@centers, x@n.bin[ ,iph, ith, iT],
                 xlim=dots$xlim, ylim=ylim.nbin,
                 type=if (n.phi>1) "p" else "h",
                 col =if (n.phi>1) col else "darkgray", lwd=4,
                 pch=16, axes=FALSE, ylab="", xlab = lab$x
                 )
            box()
            at <- seq(range.nbin[1], range.nbin[2], len=3)
            if (ith==1)
              axis(2, las=1, at=at, labels=format(at, scient=-1, dig=2),
                   outer=TRUE)
            else  axis(2, las=1, at=at, labels=FALSE)
            
            if (iT==n.T) axis(1)
            if (ith==1) title(ylab="n.bin", line=5, outer=TRUE, adj=0)

          } else {
            points(x@centers, x@n.bin[ ,iph, ith, iT],
                   type="p", col=col, 
                   pch=16)
          }
        }
        screen(n.all + 2*screen-1)
      } else {
        screen(screen)
        par(oma=c(4,oma.left,oma.top,0))
      }

      ## plot emp.vario
      #if (ith==1) par(mar=c(0,6,1,1)) else
      par(mar=c(0,.5,1,.5))

      plotted.meth <- NULL  # needed for legend
      
      for (iph in 1:n.phi) {
        if (n.phi>1) col <- col.v[iph]
        if (iph==1) {
          do.call(plot, args=c(dots, list(
                          x=x@centers, y=x@emp.vario[ ,iph, ith, iT],
                          #ylim=ylim.ev, type=type, pch=19,
                          col=col, axes=FALSE, ylab=""))
                  )
          box()
          axis(2, las=1, labels=(ith==1), outer=(ith==1))
          if (!plot.nbin) axis(1)          
          if (l.theta > 1 || l.T > 1)
            legend("topleft",
                   legend=paste(
                     if (l.T>1) paste("T=", signif(x@T[iT], 3), " ", sep=""),
                     if (l.T>1 & l.theta>1) ",",
                     if (l.theta>1) paste("theta=",
                                          signif(x@theta.centers[ith], 3),
                                          sep="")
                     )
                   )
          if (ith == 1) title(ylab=ylab.ev, line=5, outer=TRUE)
          if (has.sd && plot.sd)
            legend("topright", bty="n", #col=NA, lwd=NA,
                   legend="arrows represent +-1 sd intervals")
         
        } else {
          do.call(points, args=c(dots, list(
                            x=x@centers,
                            y=x@emp.vario[ ,iph, ith, iT], col=col))
                  )   #type=type,  pch=19)
        }
        
        k <- 1  
        if(n.methods!=0){
          for( i in 1:n.methods )
            {        
              method.model <- newx[[methodnames[i]]]
              if(length(method.model@name) == 0){
                warning( paste("The method ", methodnames[i] ,
                               " was not fitted.") )
                next
              }
              phi.angles <- c(halfphi.sector, 0, -halfphi.sector)
              theta.angles <- c(halftheta.sector, 0, -halftheta.sector)
              if(!is.null(x@phi.centers)){
                x.radial <- cbind(cos(x@phi.centers[iph]+phi.angles),
                                  sin(x@phi.centers[iph]+phi.angles))
                if(!is.null(x@theta.centers))
                  x.radial <- cbind(x.radial,
                                    cos(x@theta.centers[ith] +
                                        rep(theta.angles,
                                            each= length(phi.angles))))
                
              }
              else
                x.radial <- matrix(1, nrow=1, ncol=1)              
              x.time <- NULL
              if(!is.null(x@T))
                x.time <- x@T[iT]
              x.eval <- seq(from = max(dotsRFfit$xlim[1],1e-3),
                            to = dotsRFfit$xlim[2], len = 150)

              dir2vario <- function(dir.vec){
                x.space <- as.matrix(x.eval) %*% t(dir.vec)
                for(j in (1:20)){
                  if(j>1)
                    x.space <- cbind(x.space,0)
                  vario.vals <- try(RFvariogram(x = cbind(x.space, x.time),
                                                model = method.model,
                                                grid = FALSE),
                                    silent = FALSE)
                  #xxx
                  if(!is(vario.vals, "try-error"))
                    return(vario.vals)
                }
                return(NULL)
              }              

              dummy.vals <- apply(x.radial, 1, FUN = dir2vario)
              if(!is.null(dummy.vals)){
                vario.vals <- rowMeans(as.matrix(dummy.vals))                
                do.call(points, args=c(dotsRFfit, list(
                                  x=x.eval, y=vario.vals, col=col,
                                  lty = ltyRFfit.v[k])))
                k <- k+1
                if(iph == 1)
                  plotted.meth <- c(plotted.meth, methodnames[i])
              }              
              
            }
        }      
        
        
        if (has.sd && plot.sd) {
          sdnot0 <-  x@sd[ ,iph, ith, iT] != 0
          arrows(x@centers[sdnot0],
                 x@emp.vario[sdnot0 ,iph, ith, iT] - x@sd[sdnot0 ,iph, ith, iT],
                 x@centers[sdnot0],
                 x@emp.vario[sdnot0 ,iph, ith, iT] + x@sd[sdnot0 ,iph, ith, iT],
                 code=2, angle=90, length=0.05, col=col)
        }
        
        
      }

      pifactor <- signif((x@phi.centers[1:n.phi]) / pi, 2)
      
      len.mnames <- length(plotted.meth)
      string.emp <- "empirical"
      if(len.mnames > 0)
        labels <-
          eval(parse(text=paste(
                       "c(",
                       paste('expression(paste("phi=",',
                             rep(pifactor, each = len.mnames+1),
                             ', pi, ',
                             '", ', rep(c(string.emp, plotted.meth), l.phi),
                             '"',
                             '))', collapse=","),
                       ")")))
      else
        labels <-
          eval(parse(text=paste(
                       "c(",
                       paste('expression(paste("phi=",', pifactor, ', pi))',
                             collapse=","),
                       ")")))
               
      
      if (l.phi > 1 || len.mnames > 0)# && iT==1 && ith==1) # auskommentiert auf Sebs wunsch        
        legend("bottomright", col=rep(col.v, each = len.mnames+1),
               lwd=1, pch=rep(c(19,rep(NA, len.mnames)), l.phi), bty="n",
               legend=labels, lty =
               rep(c(1, if(len.mnames==0) NULL else ltyRFfit.v[1:len.mnames]),
                   l.phi))
            
      screen <- screen+1
    }
  dots$type <- NULL
  
  if (!is.null(main)) do.call(title, args=c(dots, main=main, outer=TRUE))
  if (!is.null(xlab)) do.call(title, args=c(dots, xlab=xlab, outer=TRUE))

  close.screen(all.screens=TRUE)
}

if (FALSE){
  
  ### plots empirical variogram
  # ev - an object produced by RFempvario()
  # fit - add fitted variogram if available
  # ... - further argument for plot function
  
  # define layout
  #nf <- layout(matrix(data=c(1,2),ncol=1,nrow=2,byrow=TRUE), c(7,0), c(5,2), TRUE)
  # find dimension of plot
  xm <- c(min(ev$centers), max(ev$centers))
  ym <- c(min(ev$emp.vario), max(ev$emp.vario)*1.1)
  # plot emp.vario
  par(mar=c(1,4,2,2), las=1)
   lab <- xylabs("bin centers", NULL, units=x@coord.units)
  plot(ev$centers, ev$emp.vario, xlim=xm, ylim=ym, xlab=lab$x,
       ylab='emprical variogram', pch=19, xaxt='n',...)
  # add fitted variogram if available
  if(!is.null(fit)){
    gx <- seq(0.001, max(xm), l=100)
    if(!is.null(fit$ml$model))
      lines(gx, RFvariogram(gx, model=fit$ml$model), col=2, lwd=2)
    if(!is.null(fit$plain$model))
      lines(gx, RFvariogram(gx, model=fit$plain$model), col=3, lwd=2)
    if(!is.null(fit$sqrt.nr$model))
      lines(gx, RFvariogram(gx, model=fit$sqrt.nr$model), col=4, lwd=2)
    if(!is.null(fitsd.inv$model))
      lines(gx, RFvariogram(gx, model=fit$sd.inv$model), col=6, lty=2, lwd=2)
    legend('bottomright', c("empirical", "ML", "least Squares", "sqrt(n) least Squares", 'sd^(-1) least Squares'), lty=c(-1, rep(1, 4)), lwd=c(-1,rep(2,4)), pch=c(19, rep(-1, 4)), col=c(1, 2, 3, 4, 6), bty='n')
  }
  # add info with bins
  par(mar=c(3,4,0,2))
  lab <- xylabs("bin centers", NULL, units=x@coord.units)
  plot(x=ev$centers, y=ev$n.bin, ylab='n.bin', xlab=lab$x,
       type='h', lwd=5, col='darkgrey')
  mtext(line=-10,'distance')
  
 }

