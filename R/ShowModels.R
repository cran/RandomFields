
## to.do: total screen as screen 0 then clicking on one of the windows will let you enter into the respective
##        window and then change things
##
##        starting: according to the values given set the plots, and get into one of the plots
##        show sign in which window you are , e.g. red bullet
##        better use Koordinates to jump directly into the other screen!

ShowModels <- function(covx=ifelse(is.null(empirical),4,max(empirical$c)),
                       x=NULL, y=NULL,
                       fixed.rs=FALSE,
                       method=NULL,
                       empirical=NULL,
                       model=NULL,
                       param=NULL,
                       all.param=NULL,
                       PracticalRange=FALSE,   
                       legend = TRUE,
                       register=0,
                       ...){

  stopifnot(all(covx>=0))
  if (length(covx)>1) {
    covx <- covx[covx>0]
    covx <- c(0, sort(c(max(covx)/10000, covx)))
  } else {
    covx <- c(0,seq(covx/10000, covx, l=100))
  }
  
  scr <- split.screen(figs=rbind(
                        c(0.01,0.60,0.01,0.49),
                        c(0.62,0.99,0.01,0.49),
                        c(0.01,0.60,0.51,0.99),
                        c(0.62,0.99,0.51,0.99)))
  if (isnullX <- is.null(x)) {
    if (!is.null(y)) stop("x is null, but not y")
  }

  if (!is.null(y)) {
    ll <-  list(...)
    opts <- lapply(names(ll), pmatch,c("col","zlim"),no=0)
    if (length(opts)==0) opts <- 0 else opts <- as.integer(opts)

    if (sum(index <- (opts==1))==1) col <- ll[index][[1]]
    else {
      if (exists("image.default")) ## R-1.3.0 onwards
        text <- paste("col <-",
                        paste(as.character(as.list(args(image.default))$col),
                              collapse="("),")")
     else  ## R-1.2.3
       text <- paste("col <-",
                    paste(as.character(as.list(args(image))$col),
                          collapse="("),")")         
      eval(parse(text=text))
  }
    if (!(missing.zlim <- (sum(index <- (opts==2))!=1))) zlim <- ll[index][[1]]
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
    lu.y <- max(y)
  }
  runif(1)
  bg.save <- par()$bg
  par(bg="white")
  
  rs <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  screen(scr[4])
  maxrow<- 14
  namen <- GetModelNames()
  n <- length(namen)
  par(mar=c(0,0,2,0))
  ncol <- 1 + as.integer( (n-1) / maxrow)
  top <- min(maxrow+1,n,na.rm=TRUE)
  plot(-1,-1,xlim=c(0,ncol),ylim=c(1,top+1),axes=FALSE)
  maintitle <- "Models"
  title(main=maintitle)
  for (i in 1:n) {
    text(as.integer((i-1) / maxrow), top - (i-1) %% maxrow,
         labels=paste(namen[i]),adj=c(0,0))
  }
  p <- .C("GetParameterIndices",integer(1),integer(1),integer(1),
          scale=integer(1),
          kappa=integer(1),lastkappa=integer(1),integer(1),integer(1))
  p <- lapply(p,function(k) k+1)
  nkappa <- p$lastkappa - p$kappa + 1
  
  globals <- rbind(c("Variogram","CovarianceFct"),
                   c("practical range","math. def")) ##"" to note numericals
  shortglobals <- rbind(c("Vario","CovFct"),
                   c("pract","math")) ##"" to note numericals
  globmin <- c(1,1)
  globmax <- c(2,2)
  globminsteps <- c(1,1)
  globmaxsteps <- c(1,1)
  if (is.null(empirical)) {
    glob <- c(2,2-PracticalRange)
    maxsteps <- c(1,1,1,1,rep(1,nkappa))
  } else {
    glob <- c(1,2-PracticalRange)
    m.vario <- mean(empirical$e)
    m.lag <- mean(empirical$c)
    maxsteps <- c( c(sqrt(m.vario), m.vario, m.vario, m.lag)/5, rep(1,nkappa))
    screen(scr[3])
    par(mar=c(2,2,0,0))
    ylim <- range(empirical$e, na.rm=TRUE)
    plot(empirical$c, empirical$e, pch="*",xlim=range(covx, na.rm=TRUE),
         ylim=ylim) 
  } 
  nglob <- length(glob)
  mins <- c(-Inf,0,0,0.001,rep(-Inf,nkappa))
  maxs <- c(Inf,Inf,Inf,Inf,rep(Inf,nkappa))
  minsteps <- maxsteps / 10000
  pnames <- c("mean","variance","nugget","scale","a","b","c","d","e")
  
  currentparam <- list()
  if (unknown.par <- is.null(all.param)) {
    if (is.null(empirical)) all.param <- c(0,1,0,1)
    else all.param <- c(0, m.vario, 0, m.lag/3)
  } else stopifnot(length(all.param)==4)

  
  for (i in 1:n) 
    currentparam[[i]] <-  c(all.param,rep(1,.C("GetNrParameters",
                                             as.integer(i-1),
                                             k=integer(1))$k))
  i <- .C("GetModelNr","cone",nr=integer(1))$nr + 1
  currentparam[[i]] <- c(all.param,rep(0.5,.C("GetNrParameters",
                                          as.integer(i-1),
                                          k=integer(1))$k))
  currentparam[[.C("GetModelNr","cauchy",nr=integer(1))$nr+1]] <-
    c(all.param,2)
  currentparam[[.C("GetModelNr","cauchytbm",nr=integer(1))$nr+1]] <-
    c(all.param,1,1,2)
  currentparam[[.C("GetModelNr","power",nr=integer(1))$nr+1]] <-
    c(all.param,2)
  
  if (unknown.par) { all.param[p$scale] <- 0.3 }                     
  currentparam[[.C("GetModelNr","bessel",nr=integer(1))$nr+1]] <-
    c(all.param,1)
  currentparam[[.C("GetModelNr","wave",nr=integer(1))$nr+1]] <- c(all.param)
  currentparam[[.C("GetModelNr","gengneiting",nr=integer(1))$nr+1]] <-
    c(all.param,1,3)

  special <- rep(0,n)
  paramAdiscrete1 <- 1
  special[.C("GetModelNr","gengneiting",nr=integer(1))$nr+1] <- paramAdiscrete1
  
  covlist <- c("bessel","cauchy","cauchytbm","circular",
               "cone","cubic","exponential","gaussian",
               "gencauchy","gengneiting","gneiting","gneitingdiff",
               "holeeffect",
               "hyperbolic","nugget","pentamodel","power",
               "qexponential","spherical","stable","wave",
               "whittlematern")
  exprlist <- c(expression(2^a *Gamma(a+1)*x^{-a}* J[a](x)),
                expression((1+x^2)^{-a}),
                expression((1+(1-b/c)*x^a)*(1+x^a)^{-b/a-1}),
                expression((1-2/pi*(x *sqrt(1-x^2)+asin(x))) * (x<1)),
                
                expression("C(x) not available"),
                expression((1-7*x^2 + 8.75 *x^3 - 3.5* x^4+0.75 *x^6)*(x<1)),
                expression(e^{-x}),
                expression(exp(-x^2)),
                
                expression((1+x^a)^{-b/a}),
                expression("formula see help page for CovarianceFct"),
                expression((1 + 8 *s *x + 25* s^2* x^2 + 32*
                    s^3 *x^3)*(1-s* x)^8 * (x<1)),
                expression((1 + 8 *x/b + 25 *(x/b)^2 + 32* (x/b)^3)*(1-x/b)^8
                    * 2^{1-a} *Gamma(a)^{-1}* x^a* K[a](x)*(x<1) ),

                expression(e^{-a* x}* cos(x)),
                
                expression(const[a][b][c] * (c^2 +x^2)^{b/2} *
                    K[b](sqrt(a(c^2 + x^2)))),
                expression((x==0)),
                expression(( 1 - 22*x^2/3  +33 *x^4  - 77*x^5/2  + 33* x^7/2
                    - 11* x^9/2 + 5 *x^11 /6)*(x<1)),
                expression((1-x)^a  * (x<1)),
                
                expression((2*e^{-x}-a*e^{-2*x})/(2-a)),
                expression((1 - 1.5* x + 0.5* x^3)*(x<1)),
                expression(exp(-x^a)),
                expression(sin(x)/x ),
                    
                expression(2^{1-a}* Gamma(a)^{-1}* x^a * K[a](x))
                )
  expr <- NULL
  for (i in 1:n) {
    j <- pmatch(namen[i],covlist)
    if (!is.na(j)) expr[i] <- exprlist[j] else expr[j] <- "C(x) unknown"
  }
  DeleteRegister(register)
  
  oldPracticalRange <-  RFparameters()$PracticalRange
  RFparameters(PracticalRange=2-glob[2])

  first <- TRUE
  new.text <- ""
  repeat {
    #screen(scr[4]);title(main=maintitle,col="red")
    screen(scr[4],new=FALSE)

    if (first && (!is.null(model))) {
      first <- FALSE
      covnr <- as.integer(.C("GetModelNr", as.character(model),
                             nr = integer(1))$nr) + 1
     
      if ((covnr<=0) || (covnr>n)) {
        cat("warning: unknown covariance model")
        next;
      } else {
        if (!is.null(param)) currentparam[[covnr]] <- param
      }
    } else {      
      if (length(loc <-locator(1))==0) break;
      loc <- lapply(loc,floor)
      covnr <- (1+top-loc$y) + maxrow * loc$x
      if ((covnr<=0) || (covnr>n) || (loc$y<0) || (loc$y>top))  next;   
    }
    cov <- namen[covnr];

    cp <- length(currentparam[[covnr]])
    currentp <- cp  + nglob
    ##screen(scr[4]);title(main=maintitle,col="black")
                   
    textcol2 <- "red"
    locy <- nglob
    repeat {
      ##erase.screen(scr[2])
      screen(scr[2])
      par(mar=c(0,0,2,0))
      plot(-1,-1,xlim=c(0,1),ylim=c(0,currentp+1.5),axes=FALSE)
      if (textcol2=="black")  title(main=cov,col=textcol2)
      else title(paste("*",main=cov,"*"),col=textcol2)
      lines(c(0.5,0.5),c(-0.1,currentp+0.5),col="lightgrey")
     
      if (locy<nglob) {
        locy <- locy+1
        glob[locy] <- min(globmax[locy],max(globmin[locy],
              glob[locy] + 2 * (loc$x-0.5) * globmaxsteps[locy] 
                         + sign(loc$x-0.5) * globminsteps[locy]))
        if (locy==2) {
          DeleteRegister(register)
          RFparameters(PracticalRange=2-glob[locy])
        }
      } else if (locy>nglob) {
        locy <- locy - nglob
        
        if (spec <- (locy>=p$kappa) && (special[covnr]!=0)) {
          spec <- FALSE
          switch(special[covnr],
                 { ##paramAdiscrete1
                   if (locy==p$kappa) {
                     spec <- TRUE
                     currentparam[[covnr]][locy] <-
                       min(maxs[locy], max(mins[locy],
                                           currentparam[[covnr]][locy] +
                                           sign(loc$x-0.5)))
                   }
                 }
                 )
        }
        if (!spec)
          currentparam[[covnr]][locy] <- min(maxs[locy],max(mins[locy],
               currentparam[[covnr]][locy] + 2 * (loc$x-0.5) * maxsteps[locy]
                                     + sign(loc$x-0.5) * minsteps[locy]))
      }
      text(0,currentp+1,lab=expr[covnr],adj=c(0,0))
      positions <- ( 1:(currentp+1) )[-c(nglob+1)] - 1
      labels <- rep("-----",length(positions))
      labels[(1:nglob)[globals[,2]!=""]] <- shortglobals[globals[,2]!="",1]
      text(0,positions,labels=labels,adj=c(0,0),col="blue")
      labels <- rep("+++++",length(positions))
      labels[(1:nglob)[globals[,2]!=""]] <- shortglobals[globals[,2]!="",2]
      text(1,positions,labels=labels,adj=c(1,0),col="red")
      txt <- globals[cbind(c(1,2),glob)]
      txt[globals[,2]==""] <- globals[globals[,2]=="",1]
      text(0.5, (1:nglob)-1, labels=txt, adj=c(0.5,0))
      text(0.5, nglob+(1:cp), labels=new.text, adj=c(0.5,0),col="white")
      new.text <- paste(pnames[1:cp]," (", format(currentparam[[covnr]],dig=3),
                        ")",sep="")
      text(0.5, nglob+(1:cp), labels=new.text, adj=c(0.5,0))

      if (textcol2=="black") break;

      erase.screen(scr[3])
      screen(scr[3])
      par(mar=c(2,2,0,0))

      eval(parse(text=paste("covvalue <- ",globals[1,glob[1]],
                   "(covx,model=cov,param=currentparam[[covnr]],dim=2-is.null(y))")))
     if (is.na(covvalue[1])) {
        plot(0,0,col=0,axes=FALSE)
        text(0,0,label="plot not available",adj=c(0.5,0.5))
      }
      else {
        if (!is.null(empirical) && (glob[1]==1))
          ylim <- range(covvalue,empirical$e, na.rm=TRUE)
        else ylim <- range(covvalue, na.rm=TRUE)
        plot(covx[1], covvalue[1], type="p", xlim=range(covx, na.rm=TRUE),
             ylim=ylim)
        lines(covx[-1],covvalue[-1])
        if (!is.null(empirical) && (glob[1]==1)) {
          points(empirical$c, empirical$e, pch="*")
        }
      }
      
      if (!isnullX) {
        screen(scr[1])
        if (!is.null(y)) lines(range(x),range(y),col="grey")
        if (fixed.rs) assign(".Random.seed", rs, envir=.GlobalEnv)
        z <- GaussRF(x,y,model=cov,param=currentparam[[covnr]],grid=TRUE,
                     register=register,method=method)

        ##erase.screen(scr[1])
        screen(scr[1])
        par(mar=c(2,2,0,0))
        if (is.null(z)) {
          plot(0,0,col=0,axes=FALSE)
          text(0,0,label="image not available",adj=c(0.5,0.5))
        } else {          
          if (is.null(y)) { plot(x,z,...) }
          else {
            image(x,y,z,...)
            if (legend) {
              if (missing.zlim) zlim <- range(z)
              if (big.legend) ml <- c(format(mean(zlim),dig=2),filler)
              else ml <- NULL
              legend(lu.x, lu.y, y.i=y.i, x.i=0.1, 
                     legend=c(format(zlim[2],dig=2),filler,
                       ml, format(zlim[1],dig=2)),
                     lty=1, lwd=lwd, col=col[length(col):1])
            }
          }
        } # is.null(z)
      }

      repeat {
        screen(scr[2],FALSE)
        if (length(loc <-locator(1))==0) {
          locy <- nglob; textcol2 <- "black"; break; }
        locy <- floor(loc$y)
        if (locy>currentp+1)  {locy <- nglob; textcol2 <- "black"; break;}
        if ((locy>=0) && (locy<=currentp)&& (locy!=nglob) && (loc$x>=0)) break;
      }
    }
  }
 
  RFparameters(PracticalRange=oldPracticalRange)
  close.screen(scr)
  par(bg=bg.save)
  if (exists("covnr"))
    return(list(model=cov,param=currentparam[[covnr]],
                PracticalRange=as.logical(2-glob[2])))
  else return(NULL)
}
  
