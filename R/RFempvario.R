## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  

rfempirical <- function(x, y = NULL, z = NULL, T = NULL, data, grid,
			bin = NULL, 
			phi=NULL,  ## phi, number of anglular segments per PI
			theta = NULL, ## aehnlich
			deltaT = NULL, ##  deltaT[1] max abstand, deltaT[2]
			##                  gitterabstand
			distances, vdim,
                        method=METHOD_VARIOGRAM, ...
			) {
  
  ## repetition is last dimension
  
  ## bin centers will be a vector of scalar distances (in cylinder coord, e.g.)
  ## for the angles: start always with the first on negative angle, continue
  ##                 counter clockwise [0, 2pi]
  ## in 3 d the third angle is zero if vector in the (x, y) plane, positive
  ##                 angle starts with points above plane
  
  ## make sure that exactly one negative value appears, and that zero is
  ## added if bin starts with a positive value
  stopifnot(length(theta) <= 1, length(phi) <= 1)
  if ((is(data, "RFsp") || isSpObj(data)) && !missing(x))
    stop("x, y, z, T may not be given if 'data' is of class 'RFsp' or an 'sp' object")

  ## to do: distances
  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")

  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  varunits <- RFopt$coords$varunits
  call <- match.call()
  
  Z <- UnifyData(x=x, y=y, z=z, T=T, distances=distances, grid=grid,
		 RFopt = RFopt,
		 data=data,  ## allowFirstCols=FALSE, 20.12.17 -- warum war dies
		 ##             gesetzt??
		 vdim=if (missing(vdim)) NULL else vdim)
  grid <- sapply(Z$coord, function(z) z$grid)
  fft <- RFopt$empvario$fft && grid[1] && all(grid == grid[1]) &&
    method %in% c(METHOD_VARIOGRAM, METHOD_PSEUDO)


#  print(fft <- FALSE)
  
  time <- Z$has.time.comp
  
  if (Z$dist.given) stop("option distances not programmed yet.")
  
  if (missing(vdim) || length(vdim) == 0) {
    vdim <- if (!is.na(Z$vdim)) Z$vdim else 1
  } else {
        if (!is.na(Z$vdim) && vdim!=Z$vdim)
	  warning("given multivariate dimension 'vdim' does not match multivariate dimension of the data")
  }


  grid <- sapply(Z$coord, function(z) z$grid)
  data <- RFboxcox(Z$data, vdim=vdim, ignore.na=TRUE)  
  restotal <- sapply(Z$coord, function(z) z$restotal)
  spatialdim <- Z$spatialdim
  repetitions <- Z$repetitions
  sets <- length(Z$data)
  for (i in 1:sets) {
    if (Z$coord[[i]]$dist.given) stop("'distances' not programmed yet.")
    dim.data <- c(restotal[i], vdim, repetitions[i])
    dim(data[[i]]) <- dim.data
    
    if (vdim > 1 && repetitions[i] > 1) {
      dataX <- aperm(data[[i]], c(1, 3, 2)) ## now: coord, repet, vdim
      dim(dataX) <- c(dim.data[1] * dim.data[3], dim.data[2])
      variance <- cov(dataX)
      rm(dataX)
    } else {
      dim(data[[i]]) <- if (vdim == 1) prod(dim.data) else dim.data[1:2]    
      variance <- var(data[[i]])
      dim(data[[i]]) <- dim.data
    }
  }
  
  if(is.null(bin) || length(bin)==0) bin <- 20
  
  if (length(bin) == 1) {
    ## automatic bin depending on coords
    xx <- Z$coord[[1]]$x
    if(grid[1])
      bin <- seq(0, max(xx[2, ] * xx[3, ]) / 2, len = bin) 
    else {
      bin <- seq(0, sqrt(sum((apply(xx, 2, max)-apply(xx, 2, min))^2))/2,
		 len = bin)
    }
    if (RFopt$basic$printlevel >= PL_SUBIMPORTANT)
      message("Bins in RFvariogram are chosen automatically:\n", 
	      paste(signif(bin, 2), collapse=" "))
  }
  
  pseudo <- RFopt$empvario$pseudovariogram
  phi0 <- RFopt$empvario$phi0 # 0 if automatic
  theta0 <- RFopt$empvario$theta0 # 0 if automatic
  has.time.comp <- Z$has.time.comp

  phigiven <-  spatialdim > 1 && (length(phi) > 1 || (length(phi)>0 && phi>1))
  thetagiven <- spatialdim > 2 &&
    (length(theta) > 1 || (length(theta) > 0 && theta > 1))
  deltaTgiven <- length(deltaT)>0 && all(deltaT > 0)
  basic <- !(has.time.comp || phigiven || thetagiven)

  if(pseudo == TRUE) {
    ## change method from cross variogram to pseudo variogram
    if(method == METHOD_VARIOGRAM) method <- METHOD_PSEUDO else
    ## change method from cross madogram to pseudo variogram
    if(method == METHOD_MADOGRAM) method <- METHOD_PSEUDOMADOGRAM
  }
 
  ##fft <- fft && repetitions == 1 # to do ! fft should allow for repetitions  
    
  bin <- prepareBin(bin)
  stopifnot(length(bin)>=2, all(is.finite(bin)))
  if (any(diff(bin)<=0)) stop("bin must be a strictly increasing sequence")
  ##  is.null(bin) in fft : see version 3.0.12 or earlier ! to do ?! 
  
  centers <- pmax(0, (bin[-1] + bin[-length(bin)])/2)
  n.bins <- length(bin) - 1

#  Print(centers, bin)

#Print(phi0, phigiven)
  
  if (!deltaTgiven) deltaT <- c(0, 0)
  if (!phigiven) phi <-c(0, 0) else if (length(phi) == 1) phi <- c(phi0, phi)
  if (!thetagiven) theta <- c(0, 0)
  else if (length(theta) == 1) theta <- c(theta0, theta)  
  stopifnot(0 <= phi[1], 2 * pi > phi[1], 
	    0 <= theta[1], 2 * pi > theta[1], 
	    phi[2] >= 0,  phi[2] == as.integer(phi[2]), 
	    theta[2] >= 0, theta[2] == as.integer(theta[2]),
	    all(is.finite(deltaT)), all(deltaT >= 0))
  
  if (has.time.comp) {
    T.start  <- sapply(Z$coord, function(x) x$T[1])
    T.step <- sapply(Z$coord, function(x) x$T[2])
    T.len  <- sapply(Z$coord, function(x) x$T[3])
    if (sets > 1) {
      if (any(abs(diff(diff(T.step))) > 1e-15))
	stop("only data sets with the same time step allowed") #generalise todo
    }
    T <-  c(0, T.step[1], max(T.len))
  } else {
    T <-  c(1, 1, 1)
  }
  
  if (length(deltaT) == 1) deltaT <- c(deltaT, 1)
  realdelta <- deltaT[2] * T[2]
  
  timeComponent <- T[3] > 1 && deltaTgiven ## T[3] > 1 impliziert time
  stepT <-  deltaT[2] / T[2]
  if (stepT != as.integer(stepT))
    stop("deltaT not multiple of distance of temporal grid")
  stepT <- max(1, stepT)
  nstepT <- as.integer(min(deltaT[1], T[2] * (T[3]-1)) / max(T[2], realdelta))
  ##                                                             , deltaT[2]??
  n.theta <- max(1, theta[2])
  n.delta <- 1 + nstepT
  n.phibin <- n.phi <- max(1, phi[2])
  dplt <- (!fft && !basic && !has.time.comp) || ((pseudo || timeComponent) && phi[2]>0)
  if (dplt) n.phibin <- n.phibin * 2

##  Print(fft, basic, has.time.comp, pseudo, timeComponent, deltaT, deltaTgiven, T, phi, n.phi, n.phibin, dplt); 
  
  totalbinsOhnevdim <- as.integer(n.bins * n.phibin * n.theta * n.delta)
  totalbins <- totalbinsOhnevdim * vdim^2
  
  phibins <- thetabins <- Tbins <- NULL
  
  if (timeComponent) Tbins <- (0:nstepT) * realdelta 
  if (phi[2] > 0) phibins <- phi[1] + 0 : ((n.phibin - 1)) * pi / n.phi

  if (n.theta > 1)
    thetabins <- theta[1] + (0 : (n.theta-1) + 0.5) * pi / n.theta
  
  dims <- c(bins=n.bins, phi=n.phibin, theta=n.theta, delta=n.delta,
	    vdim=rep(vdim, 2))
  
  empirical.sd <- NULL

  if (fft) {
    ## to do: das liest sich alles irgendwie komisch
    maxspatialdim <- 3
    
    if (Z$spatialdim > maxspatialdim)
      stop("fft does not work yet for spatial dimensions greater than ",
	   maxspatialdim)
    
    empirical <- n.bin <- 0
    for (i in 1:sets) {
      xx <- Z$coord[[i]]$x
      if (ncol(xx)<maxspatialdim)  # not matrix(0, ...) here!
        ##                              since x is a triple
        xx <- cbind(xx, matrix(1, nrow=nrow(xx), ncol=maxspatialdim-ncol(xx)))
      T3 <- if (has.time.comp) Z$coord[[i]]$T[3] else 1
      neudim <- c(xx[3, ], if (has.time.comp) T3)
            
      ## last: always repetitions
      ## last but: always vdim
      ## previous ones: coordinate dimensions
      dim(data[[i]]) <- c(neudim,vdim, length(data[[i]]) / vdim / prod(neudim))

      ## to achieve a reflection in x and z instead of y we transpose the
      ## array
      crossvar <- doVario(X=data[[i]], asVector=TRUE, pseudo=pseudo,
                          has.time.comp=has.time.comp)
      sumvals <- crossvar[[1]]
      nbvals <- crossvar[[2]]
      
      back <- .Call(C_fftVario3D, as.double(xx), 
		    as.double(sumvals), as.double(nbvals), 
		    as.double(bin), as.integer(n.bins), 
		    as.integer(T3), 
		    as.integer(stepT), as.integer(nstepT),       
		    as.double(phi), 
		    as.double(theta), 
		    as.integer(repetitions[i]),
		    as.integer(vdim),
		    totalbinsOhnevdim,
		    as.logical(pseudo) )       
      
      ## the results are now reformatted into arrays
      ## the angles are given in clear text
      
      n.bin <- n.bin + back[, EV_FFT_N + 1]
      empirical <- empirical + back[, EV_FFT_EV + 1]# back contains only sums, not averages
    }
    empirical <- empirical / n.bin ## might cause 0/0, but OK
    n.bin <- as.integer(round(n.bin))   
  } else { ## ! fft
    ## #####################################################################
    ##
    ## MARTINS CODE WENN FFT == FALSE
    ##
    ## #####################################################################
    
    if (basic) {
      n.bin <- empirical.sdSq <- empirical <- 0
      
      for (i in 1:sets) {
	back <- .Call(C_empirical, 
		    as.double(Z$coord[[i]]$x), ## Z definition
		    as.integer(spatialdim),
		    as.integer(Z$coord[[i]]$l), 
		    as.double(data[[i]]),
		    as.integer(repetitions[i]), as.integer(grid[i]), 
		    as.double(bin), as.integer(n.bins),
		    as.integer(vdim),
		    as.integer(method) )
	n.new <- back[, EV_N + 1]
	n.bin <- n.bin + n.new
	dummy <- back[, EV_EV + 1]    
        dummy[(is.na(dummy) & (centers==0)) | n.new == 0] <- 0
        empirical <- empirical + dummy * n.new
	
        dummy <- back[, EV_SDSQ + 1]
        dummy[n.new == 0] <- 0
        empirical.sdSq <- empirical.sdSq + dummy * n.new	
      }
      dummy <- n.bin != 0
      empirical[dummy] <- empirical[dummy] / n.bin[dummy]
      empirical.sd <- sqrt(empirical.sdSq / n.bin)      
      rm("back")
    } else { ## anisotropic space-time
      ## always transform to full 3 dimensional space-time coordinates
      ## with all angles given. Otherwise there would be too many special
      ## cases to treat in the c program. However, there is some lost
      ## of speed in the calculations...
      
      for (i in 1:sets) {
	ll <- Z$coord[[i]]$l
	coord <-  Z$coord[[i]]
	
	xx <- coord$x
	stopifnot(is.matrix(xx))
	if (ncol(xx)<3)  # not matrix(0, ...) here! since x could be a triple
	xx <- cbind(xx, matrix(1, nrow=nrow(xx), ncol=3-ncol(xx)))   
	
	## x fuer grid und nicht-grid: spalte x, y, bzw z
	n.bin <- empirical.sdSq <- empirical <- 0
	back <-
          .Call(C_empvarioXT,
		as.double(xx),
		as.double(if (length(coord$T)>0) coord$T else rep(1,3)),
		as.integer(Z$coord[[i]]$l), 
		as.double(data[[i]]),
		as.integer(repetitions[i]),
		as.integer(grid[i]), 
		as.double(bin), as.integer(n.bins), 
		as.double(phi[1:2]), 
		as.double(theta[1:2]), 
		as.integer(c(stepT, nstepT)), 
		## input : deltaT[1] max abstand, deltaT[2]: echter gitterabst.
		##   c   : delta[1]: index gitterabstand, deltaT[2]:#of bins -1
		##                   (zero is the additional distance)
		as.integer(vdim),
		as.integer(method)
		)
	
	n.bin <- n.bin + back[, EV_N + 1]
	dummy <- back[, EV_EV + 1]    
        dummy[(is.na(dummy) & (centers==0)) | back[, EV_N + 1] == 0] <- 0
        empirical <- empirical + dummy * back[, EV_N + 1]
	
        dummy <- back[, EV_SDSQ + 1]
        dummy[back[, EV_N + 1] == 0] <- 0
        empirical.sdSq <- empirical.sdSq + dummy^2 * back[, EV_N + 1]
	
	rm("back")

	if (FALSE) {
	  if (!has.time.comp && vdim == 1) {
	    ## vario is symmetric in phi;
	    ## so the number of phi's can be halfened in this case
	    dim(empirical) <- dims
	    dim(n.bin) <- dims
	    dim(empirical.sdSq) <- dims
	    
	    if (dims[2] > 1) {
	      dims[2] <- as.integer(dims[2] / 2)
	      half <- 1 : dims[2]
	      n.bin <- n.bin[, half,,,,, drop=FALSE] +n.bin[, -half,,,,,drop=FALSE]
	      empirical <- empirical[, half, , , , , drop=FALSE] +
	      empirical[, -half, , , , , drop=FALSE]
	      empirical.sdSq <- empirical.sdSq[, half, , , , , drop=FALSE] +
	      empirical.sdSq[, -half, , , , , drop=FALSE]
	      phibins <- phibins[half]
	    }
	  }
	} ## end false

	
      } ## sets
	
      idx <- n.bin > 1 & !is.nan(empirical) & empirical != 0 
      evsdSq <- empirical.sdSq[idx] / n.bin[idx]
      
      if (any(evsdSq < -1e-14)) {
	Print(idx, n.bin[idx] - 1, empirical.sdSq[idx], #
	      empirical.sdSq[idx] / (n.bin[idx] - 1), empirical)
	warning(paste(evsdSq))
      }
      evsdSq[evsdSq < 0] <- 0
      empirical.sd[idx] <- sqrt(evsdSq)   
      empirical.sd[!idx] <- NaN
    }
 
      
    ## ################################################################
    ##
    ## END OF MARPINS CODE WENN FFT == FALSE
    ##
    ## ################################################################

  } # !fft

  dim(empirical) <- dims
  dim(n.bin) <- dims
  if (!is.null(empirical.sd)) dim(empirical.sd) <- dims

  name <- list()
  namedim <- names(dims)
    for (i in 1:length(dims)) {
      name[[i]] <-
      if (namedim[i] %in% c("vdim1", "vdim2")) {
	if (length(Z$varnames) == 0) NULL
	else rep(Z$varnames, length.out=dims[i])
      } else if (namedim[i] != "bins") paste(namedim[i], 1:dims[i], sep="")  
    }
  dimnames(empirical) <- name
  ##  {} else names(empirical) <- Z$varnames[1]


  if (RFopt$general$spConform) {
    l <- new("RFempVariog",
	     centers=centers,
	     empirical=empirical,
	     var=variance,
	     sd= empirical.sd,
	     n.bin=n.bin,
	     phi.centers=phibins,
	     theta.centers=thetabins,
	     T=Tbins,
	     vdim = vdim,
	     coordunits = Z$coordunits,
	     varunits = varunits,
	     call=call,
	     method=method)
  } else {
    l <- list(centers=centers,
	      empirical=empirical,
	      var=variance,
	      sd= empirical.sd,
	      n.bin=n.bin,
	      phi.centers=phibins,
	      theta.centers=thetabins,
	      T=Tbins,
	      vdim = vdim,
	      coordunits =  Z$coordunits,
	      varunits = varunits,
	      call=call,
	      method=method
	      )
    class(l) <- "RF_empVariog"
  }
  
  return(l)
  
} # function rfempirical


## ############################################
## END OF MAIN FUNCTION 
## ############################################



doVario <- function(X, asVector=FALSE, pseudo=FALSE, has.time.comp=FALSE) {
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
  refl.order <- if(has.time.comp && !pseudo) c(1,3,4) else c(1,3,2)
  
  j <- 2
  for (i in 1:(n-1)) {
    for (k in 1:(2^(i-1))) {
      X_list[[j]] <- reflection(X_list[[k]], refl.order[i])
      j <- j + 1
    }      
  }
  
  ## to do the crossvariogram
  
  ## decide which blocks are needed
  blockidx <- rep(FALSE, 8)
  if(!has.time.comp && !pseudo){
    blockidx[1:(if (twoD) 2 else 4)] <- TRUE ## else 3 D
  } else if(has.time.comp && pseudo) {
    stop("Time component is not compatible with Pseudo variogram")
  } else { # ((has.time.comp && !pseudo) || (!has.time.comp && pseudo))
    blockidx[if (twoD) c(1:2, 5:6) else 1:8] <- TRUE
  }  
    
  for (i in c(1:len)){
    crossvar <- crossvario(X_list[[i]], pseudo=pseudo, dummy=!blockidx[i])
    if (has.time.comp) {
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
    parts<- len / (2^i)      
    positions <- 2^(i - 1)       
    for (j in 1:parts) {
      for (k in 1:positions) {
	idx <- 2* positions * j- positions + k
	if (has.time.comp) {
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

crossvario <- function(f, pseudo = FALSE, dummy = FALSE) {
  d <- dim(f)
  idx.repet <- length(d) 
  idx.vdim <- length(d) - 1
  repetvdim <- c(idx.vdim, idx.repet)
  vdim <- d[idx.vdim]
  repet <- d[idx.repet]
  CVd <- c(d[-repetvdim], repet, vdim, vdim)
  if(dummy) return(list(array(1, dim=CVd), array(1, dim=CVd)))
  
  idx <- rep(TRUE, length(d) - 2)
  idx.data <- paste("[", paste(1, ":", d, collapse=", "), "]")
  idx.vario <- paste("[", paste(rep(",", length(d)-2), collapse=""), "r,i,j]")
  idx.w <- paste("[", paste(1, ":", d[-repetvdim], collapse=", "), "]")
  
  dim.coord <- 2 * d[-repetvdim]-1
  F <- If <- array(0, dim=c(dim.coord, d[repetvdim]))
  eval(parse(text=paste("If", idx.data, "<- !is.na(f)")))
  f[is.na(f)] <- 0
  eval(parse(text=paste("F", idx.data,  "<- f")))
  LIf <- list(If)
  LF <- list(F)
  
  nbvals <- Crossvario <- array(0, CVd)
  
  for (i in 1:vdim) {
    for (j in 1:vdim) {
      for (r in 1:repet) {
	If <- do.call("[", c(LIf, idx, i, r))
	dim(If) <- dim.coord
	Ig <- do.call("[", c(LIf, idx, j, r))        
	dim(Ig) <- dim.coord
	F <- do.call("[", c(LF, idx, i, r))
	dim(F) <- dim.coord
	G <- do.call("[", c(LF, idx, j, r))
	dim(G) <- dim.coord
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
      }
    }
  }  
  return(list(Crossvario, as.array(round(nbvals))))
}


prepareBin <- function(bin) {
  if(missing(bin)) return(NULL)
  if (bin[1] > 0) {
    if (RFoptions()$basic$printlevel>1)
      message("empirical variogram: left bin border 0 added\n")
    bin <- c(0, bin)
  }
  if (bin[1]==0) bin <- c(-1, bin)
  if (bin[1] < 0) bin <- c(bin[1], bin[bin>=0])
  
  bin
}

