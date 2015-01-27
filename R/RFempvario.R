RFempiricalvariogram <- function(
 x, y = NULL, z = NULL, T = NULL, data, grid, bin, 
 phi,  ## phi, number of anglular segments per PI
 theta, ## aehnlich
 deltaT, ##  deltaT[1] max abstand, deltaT[2] : gitterabstand
 distances, vdim, ...
 ) {
  ## repetition is last dimension
  
  ## bin centers will be a vector of scalar distances  (in cylinder coord, e.g.)
  ## for the angles: start always with the first on negative angle, continue
  ##                 counter clockwise [0, 2pi]
  ## in 3 d the third angle is zero if vector in the (x, y) plane, positive
  ##                 angle starts with points above plane
  
  ## make sure that exactly one negative value appears, and that zero is
  ## added if bin starts with a positive value

  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  call <- match.call()

  if( missing(phi) ) phi <- NULL
  if( missing(theta) ) theta <- NULL
  if( missing(T) ) T <- NULL
  if( missing(bin) ) bin <- NULL
  if( missing(deltaT) ) deltaT <- NULL 
  
  variab.units <- RFopt$coords$varunits

  if ((is(data, "RFsp") || isSpObj(data)) && !missing(x))
    stop("x, y, z, T may not be given if 'data' is of class 'RFsp' or an 'sp' object")


  ## to do: distances
  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")

  neu <- CheckData(x=x, y=y, z=z, T=T, distances=distances, grid=grid,
                   data=data, fillall=TRUE, allowFirstCols=FALSE,
                   vdim=if (missing(vdim)) NULL else vdim)

#  options(str=list(strict.width="no", digits.d=3, vec.len=200))
# Print(neu)
  
  if (missing(vdim) || length(vdim) == 0) {
    vdim <- if (!is.null(neu$vdim)) neu$vdim else 1
  } else {
    if (!is.null(neu$vdim) && vdim!=neu$vdim)
      warning("given multivariate dimension 'vdim' does not match multivariate dimension of the data")
  }
  
  
  data <- neu$fulldata
  x <- neu$fullgiven$x
  stopifnot(is.null(neu$fullgiven$y), is.null(neu$fullgiven$z))
  T <- neu$fullgiven$T
  restotal <- neu$fullgiven$restotal
  grid <- neu$fullgiven$grid
  spatialdim <- neu$fullgiven$spatialdim
  repetitions <- as.integer(length(data) / (restotal * vdim))
  
  if (repetitions==0) stop("no data given")
  if (length(data) != restotal * vdim * repetitions)
    stop("number of data does not match coordinates")
  dim.data <- c(restotal, vdim, repetitions)
  dim(data) <- dim.data

#  Print(data); lll

  if (vdim > 1 && repetitions > 1) {
    dataX <- aperm(data, c(1, 3, 2)) ## now: coord, repet, vdim
    dim(dataX) <- c(dim.data[1] * dim.data[3], dim.data[2])
    variance <- cov(dataX)
    rm(dataX)
  } else {
    dim(data) <- if (vdim == 1) prod(dim.data) else dim.data[1:2]    
    variance <- var(data)
    dim(data) <- dim.data
  }
 
  if(is.null(bin) || length(bin)==0) bin <- 20
  if (length(bin) == 1) {
    ## automatic bin depending on coords
    if(grid)
      bin <- seq(0, max(x[2, ] * x[3, ]) / 2, len = bin) 
    else {
      bin <- seq(0, sqrt(sum((apply(x, 2, max)-apply(x, 2, min))^2))/2,
                 len = bin)
    }
    if (RFopt$general$printlevel >= PL_SUBIMPORTANT)
      message("Bins in RFempiricalvariogram are chosen automatically:\n", 
              paste(signif(bin, 2), collapse=" "))
  }


#  Print(RFopt$empvario, neu)
    
  fft <- RFopt$empvario$fft && grid
  pseudo <- RFopt$empvario$pseudovariogram
  phi0 <- RFopt$empvario$phi0 # 0 if automatic
  theta0 <- RFopt$empvario$theta0 # 0 if automatic
  time <- !is.null(T)

  thetagiven <- !is.null(theta) && spatialdim > 2
  phigiven <- !is.null(phi) && spatialdim > 1
  deltagiven <- !is.null(deltaT) && all(deltaT > 0)
  basic <- !(time || phigiven || thetagiven)
 
  if (time && !fft) stop("currently time components are not possible") ## to do
  ## to do multivariate;
    
 if(time && pseudo)
    stop("Time component is not compatible with Pseudo variogram") # to do
  if(!fft && pseudo) ## to do
    stop("Pseudo variogram is only implemented for grids with the FFT method")
  ## IS THE FFT FLAG SET

  #fft <- fft && repetitions == 1 # to do ! fft should allow for repetitions  
  
  bin <- prepareBin(bin)
  stopifnot(length(bin)>=2, all(is.finite(bin)))
  if (any(diff(bin)<=0)) stop("bin must be a strictly increasing sequence")
   ##  is.null(bin) in fft : see version 3.0.12 or earlier ! to do ?! 
  
  centers <- pmax(0, (bin[-1] + bin[-length(bin)])/2)
  n.bins <- length(bin) - 1
  
  
  if (!time) T <-  c(1, 1, 1)
  phi <- if (!phigiven) c(0, 0) else c(phi0, phi)        
  theta <- if (!thetagiven) c(0, 0) else c(theta0, theta)
  if (!deltagiven) deltaT <- c(0, 0)
  stopifnot(0 <= phi[1], 2 * pi > phi[1], 
            0 <= theta[1], 2 * pi > theta[1], 
            phi[2] >= 0,  phi[2] == as.integer(phi[2]), 
            theta[2] >= 0, theta[2] == as.integer(theta[2]),
            all(is.finite(deltaT)), all(deltaT >= 0))
  realdelta <- deltaT[2]

  NotimeComponent <- T[3]==1 || !deltagiven
  stepT <-  deltaT[2] / T[2]
  if (stepT != as.integer(stepT))
    stop("deltaT not multiple of distance of temporal grid")
  stepT <- max(1, stepT)
  nstepT <- as.integer(min(deltaT[1], T[2] * (T[3]-1)) / max(T[2], deltaT[2]))
  n.theta <- max(1, theta[2])
  n.delta <- 1 + nstepT
  n.phi <- max(1, phi[2])
  if (!fft && !basic) {
    n.phibin <- 2 * n.phi
  } else {
    n.phibin <-
      if (!pseudo && NotimeComponent) max(1, n.phi)
      else if(phi[2]==0) 1 else 2 * n.phi
  }

 # Print(n.phi, n.phibin)

  totalbinsOhnevdim <- as.integer(n.bins * n.phibin * n.theta * n.delta)
  totalbins <- totalbinsOhnevdim * vdim^2

  ##emp.vario <- double(totalbins)
  ##  n.bin <- if (fft) double(totalbins) else integer(totalbins)

  phibins <- thetabins <- Tbins <- NULL

  if (!NotimeComponent) Tbins <- (0:nstepT) * realdelta
  if (phi[2] > 0) phibins <- phi[1] + 0 : (n.phibin - 1) * pi / n.phi

  if (n.theta > 1)
    thetabins <- theta[1] + (0 : (n.theta-1) + 0.5) * pi / n.theta
 
  dims <- c(bins=n.bins, phi=n.phibin, theta=n.theta, delta=n.delta,
            vdim=rep(vdim, 2))
  # Print(dims, fft)

  emp.vario.sd <- NULL

  
  if (fft) {
    #Print(neu)

    ## to do: das liest sich alles irgendwie komisch
    maxspatialdim <- 3
        
    if (ncol(x) > maxspatialdim)
      stop("fft does not work yet for spatial dimensions greater than ",
           maxspatialdim)
    if (ncol(x)<maxspatialdim)  # not matrix(0, ...) here!
      ##                              since x is a triple
     x <- cbind(x, matrix(1, nrow=nrow(x), ncol=maxspatialdim-ncol(x)))
     
    neudim <- c(x[3, ], if (time) T[3])
      
    ## last: always repetitions
    ## last but: always vdim
    ## previous ones: coordinate dimensions
    dim(data) <- c(neudim, dim.data[-1])
    
      
    ## to achieve a reflection in x and z instead of y we transpose the
    ## array
    crossvar <- doVario(X=data, asVector=TRUE, pseudo=pseudo, time=time)
    sumvals <- crossvar[[1]]
    nbvals <- crossvar[[2]]
    
    back <- .Call("fftVario3D", as.double(x), 
                  as.double(sumvals), as.double(nbvals), 
                  as.double(bin), as.integer(n.bins), 
                  as.integer(T[3]), 
                  as.integer(stepT), as.integer(nstepT),       
                  as.double(phi), 
                  as.double(theta), 
                  as.integer(repetitions),
                  as.integer(vdim),
                  totalbinsOhnevdim,
                  as.logical(pseudo), 
                  PACKAGE="RandomFields")
   
    ## the results are now reformatted into arrays
    ## the angles are given in clear text
 
#    Print("end fftVario3D");
#    dim(emp.vario) <- c(length(emp.vario) / 4 , 4)
#    dim(n.bin) <- c(length(n.bin) / 4 , 4)
#    print(emp.vario);    print(n.bin)
       
    emp.vario <- back[, 1] / back[, 2] ## might cause 0/0, but OK
    n.bin <- as.integer(round(back[, 2]))
 
 #   Print("Xend fftVario3D")
   
  } else {   
    
    ## #####################################################################
    ##
    ## MARTINS CODE WENN FFT == FALSE
    ##
    ## #####################################################################

    if (vdim > 1) stop("multivariat only progrmmed for fft up to now")

     if (basic) {           
      back <- .C("empiricalvariogram", 
                 as.double(x), ## neu definition
                 as.integer(spatialdim), as.integer(neu$fullgiven$l), 
                 as.double(data),
                 as.integer(repetitions), as.integer(grid), 
                 as.double(bin), as.integer(n.bins), as.integer(0), 
                 emp.vario = double(totalbins), emp.vario.sd=double(totalbins), 
                 n.bin= integer(totalbins), PACKAGE="RandomFields")
      n.bin <- back$n.bin
      emp.vario.sd <- back$emp.vario.sd
      emp.vario <- back$emp.vario      
      emp.vario[is.na(emp.vario) & (centers==0)] <- 0
      rm("back")
    } else { ## anisotropic space-time
      ## always transform to full 3 dimensional space-time coordinates
      ## with all angles given. Otherwise there would be too many special
      ## cases to treat in the c program. However, there is some lost
      ## of speed in the calculations...

      stopifnot(is.matrix(x))
      if (ncol(x)<3)  # not matrix(0, ...) here! since x could be a triple
        x <- cbind(x, matrix(1, nrow=nrow(x), ncol=3-ncol(x)))
       
      back <-
        .C("empvarioXT", as.double(x), as.double(T), as.integer(neu$fullgiven$l), 
           as.double(data), as.integer(repetitions), as.integer(grid), 
           as.double(bin), as.integer(n.bins), 
           as.double(c(phi[1], phi[2])), 
           as.double(c(theta[1], theta[2])), 
           as.integer(c(stepT, nstepT)), 
           ## input : deltaT[1] max abstand, deltaT[2] : echter gitterabstand, 
           ##   c   : delta[1] : index gitterabstand, deltaT[2] : # of bins -1
           ##                   (zero is the additional distance)
           emp.vario=double(totalbins), emp.vario.sd = double(totalbins),  
           n.bin = integer(totalbins),  PACKAGE="RandomFields")
      n.bin <- back$n.bin
      emp.vario.sd <- back$emp.vario.sd
      emp.vario <- back$emp.vario
      rm("back")
      
      if (!time && vdim == 1) {
        ## vario is symmetric in phi;
        ## so the number of phi's can be halfened in this case
        dim(emp.vario) <- dims
        dim(n.bin) <- dims
        dim(emp.vario.sd) <- dims
        
        ##   Print(dims, "here"); print(n.bin[, 1, 1,1,,]); print(emp.vario[, 1, 1,1,,])
        
        if (dims[2] > 1) {
          dims[2] <- as.integer(dims[2] / 2)
          half <- 1 : dims[2]
          n.bin <- n.bin[, half,,,,, drop=FALSE] +n.bin[, -half,,,,, drop=FALSE]
          emp.vario <- emp.vario[, half, , , , , drop=FALSE] +
            emp.vario[, -half, , , , , drop=FALSE]
          emp.vario.sd <- emp.vario.sd[, half, , , , , drop=FALSE] +
            emp.vario.sd[, -half, , , , , drop=FALSE]
          phibins <- phibins[half]
        }
      }

  
    
      emp.vario <- emp.vario / n.bin ## might cause 0/0, but OK
    
      idx <- n.bin > 1 & emp.vario != 0    
      evsd <- emp.vario.sd[idx] / (n.bin[idx] - 1) -
        n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2
      if (any(evsd < -1e-14)) {
        Print(idx, n.bin[idx] - 1, emp.vario.sd[idx], #
              emp.vario.sd[idx] / (n.bin[idx] - 1), #
              emp.vario.sd[idx] / (n.bin[idx] - 1) -
              n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2,
              emp.vario)
        warning(paste(evsd))
      }
      evsd[evsd < 0] <- 0
      emp.vario.sd[idx] <- sqrt(evsd)   
      emp.vario.sd[!idx] <- NaN
    }

    ## ################################################################
    ##
    ## END OF MARPINS CODE WENN FFT == FALSE
    ##
    ## ################################################################


  } # !fft
      
  
  dim(emp.vario) <- dims
  dim(n.bin) <- dims
  if (!is.null(emp.vario.sd)) dim(emp.vario.sd) <- dims
 
#  Print(emp.vario, n.bin, dims, phibins, phibins-pi)
   
#  if (is.array(emp.vario) && length(dims) > 2) {
  name <- list()
  namedim <- names(dims)
  for (i in 1:length(dims)) {
  #  Print(namedim[i], namedim[i] %in% c("vdim1", "vdim2"))
    name[[i]] <-
      if (namedim[i] %in% c("vdim1", "vdim2")) {
          if (length(neu$variab.names) == 0) NULL
          else rep(neu$variab.names, length.out=dims[i])
      } else if (namedim[i] != "bins") paste(namedim[i], 1:dims[i], sep="")  
  }
  dimnames(emp.vario) <- name
#  } else names(emp.vario) <- neu$variab.names[1]

#  Print(emp.vario, fft)

  if (RFopt$general$spConform)
    l <- new("RFempVariog",
            centers=centers,
            emp.vario=emp.vario,
            var=variance,
            sd= emp.vario.sd,
            n.bin=n.bin,
            phi.centers=phibins,
            theta.centers=thetabins,
            T=Tbins,
            vdim = vdim,
            coord.units = neu$fullgiven$coordunits,
            variab.units = variab.units,
            call=call)
  else {
    l <- list(centers=centers,
              emp.vario=emp.vario,
              var=variance,
              sd= emp.vario.sd,
              n.bin=n.bin,
              phi.centers=phibins,
              theta.centers=thetabins,
              T=Tbins,
              vdim = vdim,
              coord.units =  neu$fullgiven$coordunits,
              variab.units = variab.units
           )
    class(l) <- "RF_empVariog"
  }
  
  return(l)
   
} # function RFempiricalvariogram


## ############################################
## END OF MAIN FUNCTION 
## ############################################



doVario <- function(X, asVector=FALSE, pseudo=FALSE, time=FALSE) {
  dimX <- dim(X)
  idx.repet <- length(dimX) 
  idx.vdim <- length(dimX) - 1
  
  d <- length(dimX) - 2## last two dimensions are repet & vdim
  twoD <- dimX[3] == 1
  n <- d + pseudo 
  len<- 2^(n-1)

  
  numbers <- cubes <- array(dim=c(dimX[1:d], len, dimX[idx.repet],
                              rep(dimX[idx.vdim], 2)))
  X_list <- as.list(rep(NA, len))
  X_list[[1]] <- X
  
  ##reflect the data, carefully with time reflection
  refl.order <- if(time && !pseudo) c(1,3,4) else c(1,3,2)

  j <- 2
  for (i in 1:(n-1)) {
    for (k in 1:(2^(i-1))) {
      X_list[[j]] <- reflection(X_list[[k]], refl.order[i])
      j <- j + 1
    }      
  }

#  print(X_list)
 # Print(n, j, X_list); lllll
    ##do the crossvariogram
  
  ## decide which blocks are needed
  blockidx <- rep(FALSE, 8)
  if(!time && !pseudo){
    if(twoD)  ## 2-dim case
      blockidx[1:2] <- TRUE
    else                ## 3-dim case
      blockidx[1:4] <- TRUE
  } else if(time && pseudo) {
    stop("Time component is not compatible with Pseudo variogram")
  } else { # ((time && !pseudo) || (!time && pseudo))
    if(twoD)  ## 2-dim case
      blockidx[c(1:2, 5:6)] <- TRUE
    else                ## 3-dim case
      blockidx[1:8] <- TRUE
  }

 
  
  for (i in c(1:len)){
#    Print(i, len)
    crossvar <- crossvario(X_list[[i]], pseudo=pseudo, dummy=!blockidx[i])
    if (time) {
      cubes[,,,,i ,,,] <- crossvar[[1]]
      numbers[,,,,i ,,,] <- crossvar[[2]]
    } else {
      cubes[,,,i ,,,] <- crossvar[[1]]
      numbers[,,,i ,,,] <- crossvar[[2]]
    }
  }

  if(asVector) return(list(as.vector(cubes), as.vector(numbers)))
 
  ##revert the reflection ## currently not used as asVector
  cubes <- crossvar[[1]]
  numbers <- crossvar[[2]]
  i<- n - 1
  for (i in (n-1):1) {
#    Print(i, n)
    parts<- len / (2^i)      
    positions <- 2^(i - 1)       
    for (j in 1:parts) {
      for (k in 1:positions) {
        idx <- 2* positions * j- positions + k
        if (time) {
          cubes[,,,,idx ,,,] <- reflection(cubes[,,,,idx ,,,], i)
          numbers[,,,,idx ,,,] <- reflection(numbers[,,,,idx ,,,], i)
        } else {
          cubes[,,,idx ,,,] <- reflection(cubes[,,,idx ,,,], i)
          numbers[,,,idx ,,,] <- reflection(numbers[,,,idx ,,,], i)
        }
      }
    }
  }
  return(list(cubes, numbers))
} 

crossvario<-function(f, pseudo = FALSE, dummy = FALSE) {
  d <- dim(f)

#  Print("crossvario", d, f, dummy); 
  idx.repet <- length(d) 
  idx.vdim <- length(d) - 1
  repetvdim <- c(idx.vdim, idx.repet)
  vdim <- d[idx.vdim]
  repet <- d[idx.repet]
  CVd <- c(d[-repetvdim], repet, vdim, vdim)
  if(dummy) return(list(array(1, dim=CVd), array(1, dim=CVd)))

  
  idx <- rep(TRUE, length(d) - 2)
  idx.data <- paste("[", paste(1, ":", d, collapse=", "), "]")
  idx.vario <- paste("[", paste(rep(",", length(d)-2), collapse=""), "r, i, j]")
  idx.w <- paste("[", paste(1, ":", d[-repetvdim], collapse=", "), "]")
 
  dim.coord <- 2 * d[-repetvdim]-1
  F <- If <- array(0, dim=c(dim.coord, d[repetvdim]))

  #  Print(If, f, idx.data, dim.coord)
  
  eval(parse(text=paste("If", idx.data, "<- !is.na(f)")))
  f[is.na(f)] <- 0
  eval(parse(text=paste("F", idx.data,  "<- f")))
  LIf <- list(If)
  LF <- list(F)

  nbvals <- Crossvario <- array(0, CVd)
   
  for (i in 1:vdim) {
    for (j in 1:vdim) {
      for (r in 1:repet) {
        #
        
        
        If <- do.call("[", c(LIf, idx, i, r))
        dim(If) <- dim.coord
        Ig <- do.call("[", c(LIf, idx, j, r))        
        dim(Ig) <- dim.coord
        F <- do.call("[", c(LF, idx, i, r))
        dim(F) <- dim.coord
        G <- do.call("[", c(LF, idx, j, r))
        dim(G) <- dim.coord

        #Print(F, If, i, j, r, d)
        #Print(r, repet)
     
        if (!pseudo) {    
          fftIfIg <- fft(If * Ig)
          fftFG <- fft(F * G)
          fftIfG <- fft(G * If)
          fftIgF <- fft(F * Ig)   
          z <- fft(Conj(fftFG) * fftIfIg
                   + Conj(fftIfIg) * fftFG
                   - Conj(fftIgF) * fftIfG
                   - Conj(fftIfG) * fftIgF, inverse=TRUE)
          N <- fft( Conj(fftIfIg) * fftIfIg, inverse=TRUE )
        } else {
          F2 <- F^2
          G2 <- G^2
          fftIf <- fft(If)
          fftIg <- fft(Ig)
          z <- fft( Conj(fft(F2))* fftIg
                   + Conj(fftIf) * fft(G2)
                   - 2* Conj(fft(F)) * fft(G), inverse=TRUE)
          ## N <- 2* fft(Conj(fftIf)*fftIg, inverse=TRUE)
          N <- fft(Conj(fftIf)*fftIg, inverse=TRUE)
        }

        w <- Re(z) / (2 * prod(dim(N))) # sumvals
  
        
        eval(parse(text=paste("Crossvario", idx.vario, "<- w", idx.w)))
        eval(parse(text=paste("nbvals", idx.vario,
                     "<- Re(N", idx.w, ") / prod(dim(N))")))

#        Print(CVd, w, idx.vario, If, idx.w, r, i,j, dim(Crossvario))
      }
    }
  }  
  return(list(Crossvario, as.array(round(nbvals))))
}


prepareBin <- function(bin)
{
  if(missing(bin)) return(NULL)
  if (bin[1] > 0) {
      if (RFoptions()$general$printlevel>1)
	message("empirical variogram: left bin border 0 added\n")
      bin <- c(0, bin)
    }
  if (bin[1]==0) bin <- c(-1, bin)
  if (bin[1] < 0) bin <- c(bin[1], bin[bin>=0])
  
  bin
}

