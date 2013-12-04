RFempiricalvariogram <- function(
 x, y = NULL, z = NULL, T = NULL, data, data2, grid, bin, 
 phi,  ## phi, number of anglular segments per PI
 theta, ## aehnlich
 deltaT, ##  deltaT[1] max abstand, deltaT[2] : gitterabstand
 distances, ...
 ) {
  
  ## bin centers will be a vector of scalar distances  (in cylinder coord, e.g.)
  ## for the angles: start always with the first non negative angle, continue
  ##                 counter clockwise [0,2pi]
  ## in 3 d the third angle is zero if vector in the (x,y) plane, positive
  ##                 angle starts with points above plane
  
  ## make sure that exactly one negative value appears, and that zero is
  ## added if bin starts with a positive value 

  
  call <- match.call()

  if( missing(phi) ) phi <- NULL
  if( missing(theta) ) theta <- NULL
  if( missing(T) ) T <- NULL
  if( missing(bin) ) bin <- NULL
  if( missing(deltaT) ) deltaT <- NULL
  
  if (missing.data2 <- missing(data2))
    data2 <- data

  if (is(data, "RFsp")) { 
    stopifnot(is(data2, "RFsp"))
    data@data <- data@data[0:(data@.RFparams$n-1) * data@.RFparams$vdim + 1]
    data2@data <- data2@data[0:(data2@.RFparams$n-1) * data2@.RFparams$vdim +
                             min(2, data2@.RFparams$vdim)]

    stopifnot(all.equal(coordinates(data), coordinates(data2)))
    if (!(missing(x) && is.null(y) && is.null(z) && is.null(T)))
      stop("x, y, z, T may not be given if 'data' is of class 'RFsp'")
    
    gridtmp <- isGridded(data)
    compareGridBooleans(grid, gridtmp)
    grid <- gridtmp
    tmp <- RFspDataFrame2conventional(data)
    tmp2 <- RFspDataFrame2conventional(data2)
    x <- tmp$x
    y <- NULL
    z <- NULL
    T <- tmp$T
    data <- tmp$data
    data2 <- tmp2$data
    rm(tmp)
  } else stopifnot(!is(data2, "RFsp"))
   
 
  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  ## bitte code ueberpruefen. Hiermit soll Nutzer fast nie
  ## "grid" eingeben muessen
 
  new <- CheckXT(x, y, z, T, grid, distances,
                 if (!missing(distances)) spdim=1, length.data=length(data))
  repetitions <- length(data)/new$restotal
  
  variance <- var(data)
  if (!missing.data2) {
    variance <- matrix(ncol=2, c(variance, rep(var(data, data2), 2),
                         var(data2)))
  }
  
 
  if(is.null(bin)){
    ## automatic bin depending on coords
    
    if(new$grid)
      bin <- seq(0, max(new$x[2,] * new$x[3,]) / 2, len = 20) 
    else {
      bin <- seq(0, sqrt(sum((apply(new$x, 2, max)-apply(new$x, 2, min))^2)) /2,
                len = 20)      
    }
    message("Bins in RFempiricalvariogram are chosen automatically:\n",
            paste(signif(bin, 2), collapse=" "))
  }

  


#  Print(RFopt$empvario, new)
    
  fft <- RFopt$empvario$fft && new$grid
  pseudo <- RFopt$empvario$pseudo
  phi0 <- RFopt$empvario$phi0 # 0 if automatic
  theta0 <- RFopt$empvario$theta0 # 0 if automatic
  time <- !is.null(new$T)

  if(new$spacedim <= 2) theta <- NULL
  if(new$spacedim <= 1) phi <- NULL
 
  
  if(time && pseudo)
    stop("Time component is not compatible with Pseudo variogram")
  if(!fft && pseudo)
    stop("Pseudo variogram is only implemented for grids with the FFT method")
  ## IS THE FFT FLAG SET

   if(fft)
    {
      data<-as.array(data)
      data2<-as.array(data2)

      ## missing(bin) not implemented properly yet, instead choose bins automatically (below)
      if( (is.null(bin)) && FALSE) 
	{
          stop("Bin vector")
                    
          if(!(is.null(phi))    ||
             !(is.null(theta)))
            warning("Cannot handle angles if bin is missing!")

          data<-.incrDim3(as.array(data), time)
          data2<-.incrDim3(as.array(data2), time)
          
	  crossvar <- doVario(data, data2, asVector = FALSE, pseudo=pseudo, time = time)
          sumvals <- doMerge(crossvar[[1]])
          nbvals <- doMerge(crossvar[[2]])
          
          idx <- nbvals > 0
	  emp.vario <- array(0, dim=dim(as.array(sumvals)))
          emp.vario[idx] <- sumvals[idx] / nbvals[idx]
          emp.vario[!idx] <- NaN

         if (RFopt$general$spConform)
           l <-  new("RFempVariog",
                    centers=new,
                     emp.vario=emp.vario,
                     sd=NULL,
                     n.bin=NULL,
                     phi.centers=NULL,
                     theta.centers=NULL,
                     T=NULL,
                     call=call )
          else
           l <- list("RFempVariog",
                    centers=new,
                     emp.vario=emp.vario,
                     sd=NULL,
                     n.bin=NULL,
                     phi.centers=NULL,
                     theta.centers=NULL,
                     T=NULL
                     )

#          Print(l)
          
          return(l)
          #return(emp.vario)
	}
      
      bin<- prepareBin(bin)
      centers <- pmax(0,(bin[-1] + bin[-length(bin)])/2)
      n.bins <- length(bin) - 1 
      if(((is.null(phi)) && (is.null(theta)) && !is.null(bin) && FALSE) ||
         (new$spacedim > 3)) {
        cubes <- doVario(data, data2, asVector = FALSE, pseudo=pseudo)
        warning("Not implemented yet!")
        if(!is.null(T)) {
        } else {
        }
      } else {   #3-dimensional case with angles and bins
        stopifnot(ncol(new$x) <= 3)
          if (ncol(new$x)<3)  # not matrix(0,...) here! since x could be a triple
            new$x <- cbind(new$x, matrix(1, nrow=nrow(new$x),
                                         ncol=3-ncol(new$x)))
                       
          data<-.incrDim3(as.array(data), time)
          data2<-.incrDim3(as.array(data2), time)
          
          ## to achieve a reflection in x and z instead of y we transpose the
          ## array
          
          if (!time) {
            T <- c(1, 1, 1) 
            #data<- aperm(data, c(1,3,2))
            #data2<- aperm(data2, c(1,3,2))
          } else {
            T <- new$T            
            #data<- aperm(data, c(1,3,2,4))
            #data2<- aperm(data2, c(1,3,2,4))
          }

          NotimeComponent <- (T[3]==1 ||
                              is.null(deltaT) || any(deltaT==0))
                    
          crossvar <- doVario(data, data2, asVector=TRUE, pseudo, time)
          sumvals <- crossvar[[1]]
          nbvals <- crossvar[[2]]

          ## Prepare data for C call
       
          if (is.null(phi)) phi <- c(0, 0)
          else phi <- c(phi0, phi)         
          if (is.null(theta)) theta <- c(0, 0)
          else theta <- c(theta0, theta)
          if (is.null(deltaT)) deltaT <- c(0,0)
          stepT <- deltaT[2] / T[2]
          if (stepT != as.integer(stepT))
            stop("deltaT not multiple of distance of temporal grid")
          if (any(diff(bin)<=0))
            stop("bin must be a strictly increasing sequence")
          stopifnot(0 <= phi[1], 2 * pi > phi[1],
                    0 <= theta[1], 2 * pi > theta[1],
                    phi[2] >= 0,  phi[2] == as.integer(phi[2]), 
                    theta[2] >= 0, theta[2] == as.integer(theta[2]))
          stopifnot(all(is.finite(deltaT)), all(deltaT >= 0))
          realdelta <- deltaT[2]
        
          stepT <- max(1, stepT)
          nstepT <- as.integer(min(deltaT[1], T[2] * (T[3]-1)) / max(T[2], deltaT[2]))
          
          n.phi <- max(1, phi[2])
          if((!pseudo)&&NotimeComponent)            
            n.phibin <- max(1,n.phi)
          else
            n.phibin <- if(phi[2]==0) 1 else 2 * n.phi           
          n.theta <- max(1, theta[2])
          n.delta <- 1 + nstepT
          totalbins <- n.bins * n.phibin * n.theta * n.delta
          emp.vario <- double(totalbins)
          n.bin <- double(totalbins)
          phibins <- thetabins <- NULL
          Tbins <- 0
          if (length(data) != new$restotal * repetitions)
            stop("number of data does not match coordinates")
          ##print( list(new$x, sumvals,nbvals,repetitions,phi,pseudo))
  
        .C("fftVario3D", as.double(new$x),
             as.double(sumvals), as.double(nbvals),
             as.double(bin), as.integer(n.bins),
             as.integer(T[3]),
             as.integer(stepT), as.integer(nstepT),       
             as.double(c(phi[1], phi[2])),
             as.double(c(theta[1], theta[2])),
             as.integer(repetitions),
             emp.vario,
             n.bin,
             pseudo,
             PACKAGE="RandomFields", DUP = FALSE)
          ## the results are now reformatted into arrays
          ## the angles are given in clear text
          
          if (NotimeComponent) {
            ## no genuine time component
            if (n.theta>1) {             
              if (phi[2]>0) {                                          
                n.bin <- array(n.bin, dim=c(n.bins, n.phibin, n.theta))
                emp.vario <- array(emp.vario, dim=c(n.bins, n.phibin, n.theta))
              } else {                
                n.bin <- matrix(n.bin, nrow=n.bins, ncol=n.theta)
                emp.vario <- matrix(emp.vario, nrow=n.bins, ncol=n.theta)
              }
            } else { # n.theta==1
              if (phi[2]>0) {
                ## phi[2]!=0, theta[2]==0, no time component                  
                n.bin <- matrix(n.bin, ncol=n.phibin)
                emp.vario <-  matrix(emp.vario, ncol=n.phibin)                
              } else {
                ## nichts zu tun
              }
            }
          } else {
            ## genuine time component
            dims <- c(n.bins, n.phibin, n.theta, n.delta)
            dims <- dims[dims!=1] 
            n.bin <-  array(n.bin, dim=dims)
            emp.vario <- array(emp.vario, dim=dims)             
            Tbins <- (0:nstepT) * realdelta            
          }

          thetabins <- theta[1] +  0 : (n.theta-1) * pi / n.theta +  pi / (2 * n.theta)
          phibins <- phi[1] + 0 : (n.phibin - 1) * pi / n.phi          
          
          idx <- n.bin > 0
          emp.vario[idx] <- emp.vario[idx] / n.bin[idx]
          emp.vario[!idx] <- NaN   
          
          
#          return(list(centers=centers, emp.vario=emp.vario,
#                      n.bin=n.bin,
#                      phi.centers=phibins,
#                      theta.centers=thetabins,
#                      T=Tbins
#                      ))
        } # else        
    
    }  # if(fft)
  
  else
    {
#      Print(fft); xxx
  
     
      ## #################################################################################################
      ##
      ## MARPINS CODE WENN FFT == FALSE
      ##
      ## #################################################################################################
      
      stopifnot( #all(is.finite(data)), fuer nicht grid jetzt NA erlaubt
                length(bin)>=2, all(is.finite(bin)))
      
      bin<- prepareBin(bin)
      centers <- pmax(0,(bin[-1] + bin[-length(bin)])/2)
      n.bins <- length(bin) - 1 
      
      if (!missing(distances)) stop("option distances not programmed yet.")

      repetitions <- as.integer(length(data)/new$restotal)
      
      if (length(data) != new$restotal * repetitions) 
        stop("number of data does not match coordinates")
      if (repetitions==0) stop("no data given")
      
      if (is.null(T) && is.null(phi) && is.null(theta)) {
         
        emp.vario <- double(n.bins)
        emp.vario.sd <- double(n.bins)
        n.bin <- integer(n.bins)
##        new$x <- rfConvertToOldGrid(new$x)
        
       .C("empiricalvariogram",
           as.double(new$x), ## new definition
           as.integer(new$spacedim), as.integer(new$l), 
           as.double(data), as.integer(repetitions), as.integer(new$grid), 
           as.double(bin), as.integer(n.bins), as.integer(0), 
           emp.vario, emp.vario.sd,
           n.bin, PACKAGE="RandomFields", NAOK=TRUE, DUP = FALSE)
        emp.vario[is.na(emp.vario) & (centers==0)] <- 0
        phibins <- 0
        thetabins <- 0
        Tbins <- 0       
        
      } else { ## anisotropic space-time
        ## always transform to full 3 dimensional space-time coordinates
        ## with all angles given. Otherwise there would be too many special
        ## cases to treat in the c program. However, there is some lost
        ## of speed in the calculations...
        stopifnot(is.matrix(new$x))
        if (ncol(new$x)<3)  # not matrix(0,...) here! since x could be a triple
          new$x <- cbind(new$x, matrix(1, nrow=nrow(new$x), ncol=3-ncol(new$x)))
    ##    new$x <- rfConvertToOldGrid(new$x)
        phi <- if (is.null(phi)) c(0, 0)  else c(phi0, phi)
        theta <- if (is.null(theta)) c(0, 0)
                 else c(theta0, theta)        
        T <- if (is.null(T)) c(1, 1, 1) else new$T
        if (is.null(deltaT)) deltaT <- c(0,0)
        if (any(diff(bin)<=0)) stop("bin must be a strictly increasing sequence")
        stopifnot(0 <= phi[1], 2 * pi > phi[1], 0 <= theta[1], 2 * pi > theta[1],
                  phi[2] >= 0,  phi[2] == as.integer(phi[2]), 
                  theta[2] >= 0, theta[2] == as.integer(theta[2]))
        stopifnot(all(is.finite(deltaT)), all(deltaT >= 0))
        realdelta <- deltaT[2]
        stepT <-  deltaT[2] / T[2]
        if (stepT != as.integer(stepT))
          stop("deltaT not multiple of distance of temporal grid")
        deltaT <- c(max(1, stepT),
                    as.integer(min(deltaT[1], T[2] * (T[3]-1)) / max(T[2], deltaT[2])))
        n.phi <- max(1, 2 * phi[2])
        n.theta <- max(1, theta[2])
        n.delta <- 1 + deltaT[2]
        totalbins <- n.bins * n.phi * n.theta * n.delta
        emp.vario <- double(totalbins)
        emp.vario.sd <- double(totalbins)
        n.bin <- integer(totalbins)
        phibins <- thetabins <- NULL
        Tbins <- 0
       .C("empvarioXT", as.double(new$x), as.double(T), as.integer(new$l),
           as.double(data), as.integer(repetitions), as.integer(new$grid), 
           as.double(bin), as.integer(n.bins),
           as.double(c(phi[1], phi[2])),
           as.double(c(theta[1], theta[2])),
           as.integer(deltaT),
           ##  input : deltaT[1] max abstand, deltaT[2] : echter gitterabstand,
           ##    c   : delta[1] : index gitterabstand, deltaT[2] : number of bins -1
           ##                   (zero is the additional distance)
           emp.vario, emp.vario.sd,
           n.bin,  PACKAGE="RandomFields", NAOK=TRUE,  DUP = FALSE)
        
        ## the results are now reformatted into arrays
        ## the angles are given in clear text    
        if (T[3]==1 || is.null(deltaT) || any(deltaT==0)) {
          ## no genuine time component
          if (n.theta>1) { 
            
            ## thetabins <- theta[1] +  0 : (n.theta-1) * pi / n.theta
            if (n.phi>1) {
              n.phi <- n.phi / 2              
              n.bin <- array(n.bin, dim=c(n.bins, n.phi, 2, n.theta))
              n.bin <- n.bin[,,1,] + n.bin[,,2,]
              emp.vario <- array(emp.vario, dim=c(n.bins, n.phi, 2, n.theta))
              emp.vario <- emp.vario[,,1,] + emp.vario[,,2,]
              emp.vario.sd <- array(emp.vario.sd, dim=c(n.bins, n.phi, 2, n.theta))
              emp.vario.sd <- emp.vario.sd[,,1,] + emp.vario.sd[,,2,]
            } else {
              n.bin <- matrix(n.bin, ncol=n.theta, nrow=n.bins)
              emp.vario <- matrix(emp.vario, ncol=n.theta, nrow=n.bins)
              emp.vario.sd <- matrix(emp.vario.sd, nrow=n.bins, ncol=n.theta)
            }     
          } else { # n.theta==1, no real 3 dimensional component
            if (n.phi>1) {              
              ## phi[2]!=0, theta[2]==0, no time component
              n.phi <- n.phi /2
              half <- totalbins/2
              n.bin <- matrix(n.bin[1:half] + n.bin[(half+1):totalbins], ncol=n.phi)              
              emp.vario <-  matrix(emp.vario[1:half] + emp.vario[(half+1):totalbins], ncol=n.phi)   
              emp.vario.sd <-
                matrix(emp.vario.sd[1:half] + emp.vario.sd[(half+1):totalbins],
                       ncol=n.phi)              
            } else {
              ## nichts zu tun
            }
          }
          phibins <- phi[1] + 0 : (n.phi - 1) * pi / n.phi
          thetabins <- theta[1] +  0 : (n.theta-1) * pi / n.theta +  pi / (2 * n.theta)          
        } else {
          ## genuine time component
          dims <- c(n.bins, n.phi, n.theta, n.delta)
          dims <- dims[dims!=1]
          n.bin <-  array(n.bin, dim=dims)
          emp.vario <- array(emp.vario, dim=dims) 
          emp.vario.sd <- array(emp.vario.sd, dim=dims)
          Tbins <- (0:deltaT[2]) * realdelta
          thetabins <- theta[1] +  0 : (n.theta-1) * pi / n.theta +  pi / (2 * n.theta)
          phibins <- phi[1] + 0 : (n.phi - 1) * 2 * pi / n.phi
        }
        
        idx <- n.bin > 0
        emp.vario[idx] <- emp.vario[idx] / n.bin[idx]
        emp.vario[!idx] <- NaN
        idx <- n.bin > 1        
        
        evsd <- emp.vario.sd[idx] / (n.bin[idx] - 1) -
          n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2
        if (any(evsd < -1e-14)) {
          Print(idx, n.bin[idx] - 1, emp.vario.sd[idx], #
                emp.vario.sd[idx] / (n.bin[idx] - 1) -
                n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2)
          warning(paste(evsd))
        }
        evsd[evsd < 0] <- 0
        emp.vario.sd[idx] <- sqrt(evsd)   
        emp.vario.sd[!idx] <- NaN
        
      } ## anisotropic space-time
   
      
      ## #################################################################################################
      ##
      ## END OF MARPINS CODE WENN FFT == FALSE
      ##
      ## #################################################################################################
      
    } # if(fft) else {}
  
  sd <- if (fft) NULL else emp.vario.sd

  if(new$spacedim <= 2) thetabins <- NULL
  if(new$spacedim <= 1) phibins <- NULL
  if(!time) Tbins <- NULL  

 if (RFopt$general$spConform)
   l <- new("RFempVariog",
            centers=centers,
            emp.vario=emp.vario,
            var=variance,
            sd=sd,
            n.bin=n.bin,
            phi.centers=phibins,
            theta.centers=thetabins,
            T=Tbins,
            call=call)
  else 
    l <- list("RFempVariog",
              centers=centers,
              emp.vario=emp.vario,
              var=variance,
              sd=sd,
              n.bin=n.bin,
              phi.centers=phibins,
              theta.centers=thetabins,
              T=Tbins
            )
  
  return(l)
 
  
} # function RFempiricalvariogram


## ############################################
## END OF MAIN FUNCTION 
## ############################################

reflection <- function(data, orth, drop=FALSE)
  ##IMPORPANT NOTE! DO NOT CHANGE THE VARIABLE NAMES IN THIS SIGNATURE
  ## why ???
  ## since the variable data is pasted by its name
{
  
  len<- dim(data)[orth]
  n<- length(dim(data))
  text <- paste("data[", reps(orth-1), len, ":1", reps(n - orth),
                ",drop=", drop, "]")
 # Print(text)
  return(eval(parse(text=text)))
    
#  catch(text)

}

#ERWARPET:	
#SICHERP:		
#AUTHOR:		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
stick<- function(X, Y, orth)
{
  len<- dim(X)[orth]
  
  n<- length(dim(X))
  
  prefix <- paste("[", reps(orth-1))
  suffix <- paste(reps(n-orth))
 
  newDim <- dim(X)
  newDim[orth] <- len* 2- 1

  data<- array(NA, newDim)

  text <- paste("data", prefix, "1:", len, suffix, "]",
                "<- Y", prefix, "1:", len, suffix, ",drop=FALSE]")
  eval(parse(text=text))

  text <- paste("data", prefix, len+ 1,":", len* 2- 1, suffix, "]",
                "<- X", prefix, "2:", len, suffix, ",drop=FALSE]")
  eval(parse(text=text))	
  
  return(data)
}


doMerge <- function(cubes)
{
  n<- log2(length(cubes))+ 1

  i<- 1
  
  while (i < n)
  {
      for (j in c(1:(length(cubes)/2)))
        {
          cubes[[j]]<- stick(cubes[[j]], cubes[[j+ 1]], i)
          cubes<- cubes[-(j+ 1)]
        }
          
      i <- i+ 1
  }
  
  return(cubes[[1]])
}


doVario <- function(X, Y, asVector=FALSE, pseudo=FALSE, time=FALSE)
{
  
  n<- length(dim(X)) + pseudo
  len<- 2^(n-1)
  X_list<- list(X)
  Y_list<- list(Y)
 
  ##reflect the data, careful with time reflection
  if(time && !pseudo)
    refl.order <- c(1,3,4)
  else
    refl.order <- c(1,3,2)
  i<- 1
  while (i <= (n-1))
    {
      for (X in X_list)
        {
          X_list<- c(X_list, list(reflection(X, refl.order[i])))
        }
      for (Y in Y_list)
        {
          Y_list<- c(Y_list, list(reflection(Y, refl.order[i])))
        }
      
      i<- i+ 1
    }
  cubes<- numbers <- list()
  ##do the crossvariogram
  
  ## decide which blocks are needed
  blockidx <- array(0,dim=8)
  if(!time && !pseudo){
    if(dim(X)[3] == 1)  ## 2-dim case
      blockidx[1:2] <- 1
    else                ## 3-dim case
      blockidx[1:4] <- 1
  }
  if((time && !pseudo) || (!time && pseudo)){
    if(dim(X)[3] == 1)  ## 2-dim case
      blockidx[c(1:2,5:6)] <- 1
    else                ## 3-dim case
      blockidx[1:8] <- 1
  }
  if(time && pseudo)
    stop("Time component is not compatible with Pseudo variogram")
  
  for (i in c(1:len)){    
      crossvar <- crossvario(X_list[[i]], Y_list[[i]], pseudo, !as.logical(blockidx[i]))      
      cubes <- c(cubes, list(crossvar[[1]]))
      numbers <- c(numbers, list(crossvar[[2]]))
    }
  
  if(asVector) 
    return(list(unlist(lapply(cubes, FUN=as.vector)), unlist(lapply(numbers, FUN=as.vector))))

  ##revert the reflection
  i<- n- 1
  while (i > 0)
    {
      parts<- len/(2^i)
      
      for (j in c(1:parts))
        {
          positions<- 2^(i- 1) 
          
          for (k in c(1:positions))
            {
              cubes[[2* positions* j- positions + k]]<- reflection(cubes[[2* positions* j- positions + k]], i)
              numbers[[2* positions* j- positions + k]]<- reflection(numbers[[2* positions* j- positions + k]], i)
            }
        }
      
      i<-i- 1
    }
  
  return(list(cubes, numbers))
} 

crossvario<-function(f,g, pseudo = FALSE, dummy = FALSE)
{
  d<-dim(f)
  if(dummy) return(list(array(1,d), array(1,d)))
  
  F <- G <- If <- Ig <- array(0,2*d-1)
  eval(parse(text=paste("If[",paste(1,":",d,collapse=","),"]","<-!is.na(f)")))  
  eval(parse(text=paste("Ig[",paste(1,":",d,collapse=","),"]","<-!is.na(g)")))
  
  f[is.na(f)] <- 0
  g[is.na(g)] <- 0 
  eval(parse(text=paste("F[",paste(1,":",d,collapse=","),"]","<-f")))
  eval(parse(text=paste("G[",paste(1,":",d,collapse=","),"]","<-g")))      
  if(!pseudo){    
    fftIfIg <- fft(If*Ig)
    fftFG <- fft(F*G)
    fftIfG <- fft(G*If)
    fftIgF <- fft(F*Ig)   
    z<-fft(Conj(fftFG)*fftIfIg + Conj(fftIfIg)*fftFG - Conj(fftIgF)*fftIfG - Conj(fftIfG)*fftIgF, inverse=TRUE)
    N <- fft( Conj(fftIfIg)*fftIfIg, inverse=TRUE )
  }
  else
    {
      F2 <- F^2
      G2 <- G^2
      fftIf <- fft(If)
      fftIg <- fft(Ig)
      z <- fft( Conj(fft(F2))* fftIg + Conj(fftIf) * fft(G2)- 2* Conj(fft(F)) * fft(G), inverse=TRUE)
      ## N <- 2* fft(Conj(fftIf)*fftIg, inverse=TRUE)
      N <- fft(Conj(fftIf)*fftIg, inverse=TRUE)
    }  
  
  w <- Re(z) / (2 * prod(dim(N))) # sumvals
  Crossvario<-array(0,d)
  eval(parse(text=paste("Crossvario[",paste(1,":",d,collapse=","),"]","<-w[",paste(1,":",d,collapse=","),"]")))
  eval(parse(text=paste("nbvals<-Re(N)[",paste(1,":",d,collapse=","),"]")))     
  nbvals <- nbvals / prod(dim(N))
  return(list(Crossvario,as.array(round(nbvals))))
}

prepareBin <- function(bin)
{
  if(missing(bin))
    {
      return()
    }

  if (bin[1] > 0)
    {
      if (RFoptions()$general$printlevel>1)
	cat("empirical variogram: left bin border 0 added\n")
      bin <- c(0, bin)
    }         
  if (bin[1]==0) bin <- c(-1, bin)
  if (bin[1] < 0) bin <- c(bin[1], bin[bin>=0])
  
  bin
}


.incrDim3<- function(data, time)
{
  newDim<-dim(data)  
  timeidx <- i <- length(dim(data))
  while (i < (3+time))
    {
      newDim<-c(newDim, 1)      
      i<- i+ 1
    }
  if(time)
    newDim[c(timeidx,4)] <- newDim[c(4,timeidx)]
  return(array(data,newDim))
}

hasTimeComponent <- function(coords)
{
	return ("coords.T1" %in% dimnames(coords)[[2]])
}
