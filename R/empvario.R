"EmpiricalVariogram" <-
function (x, y = NULL, z = NULL, T=NULL, data, grid, bin, gridtriple = FALSE,
          phi,  ## phi[1] erste richtung, phi[2] : anzahl der richtungen
          theta, ## aehnlich
          deltaT ##  deltaT[1] max abstand, deltaT[2] : gitterabstand
          ) 
{
  ## bin centers will be a vector of scalar distances  (in cylinder coord, e.g.)
  ## for the angles: start always with the first non negative angle, continue
  ##                 counter clockwise [0,2pi]
  ## in 3 d the third angle is zero if vector in the (x,y) plane, positive
  ##                 angle starts with points above plane
  stopifnot(all(is.finite(data)),
            length(bin)>=2,
            all(is.finite(bin)))
  
  ## make sure that exactly one negative value appears, and that zero is
  ## added if bin starts with a positive value 
  if (bin[1] > 0) {
    if (RFparameters()$Print>1)
      cat("empirical variogram: left bin border 0 added\n")
    bin <- c(0, bin)
  }
  if (bin[1]==0) bin <- c(-1, bin)
  if (bin[1] < 0) bin <- c(bin[1], bin[bin>=0])
    
  new <- CheckXT(x, y, z, T, grid, gridtriple)
  repet <- as.integer(length(data)/new$total)
  if (length(data) != new$total * repet) 
    stop("number of data does not match coordinates")
  if (repet==0) stop("no data given")

  centers <- pmax(0,(bin[-1] + bin[-length(bin)])/2)
  n.bins <- length(bin) - 1 
  
  if (is.null(T) && missing(phi) && missing(theta)) {  
    emp.vario <- double(n.bins)
    emp.vario.sd <- double(n.bins)
    n.bin <- integer(n.bins)
    .C("empiricalvariogram",
       as.double(new$x),
       as.integer(new$spacedim), as.integer(new$l), 
       as.double(data), as.integer(repet), as.integer(new$grid), 
       as.double(bin), as.integer(n.bins), as.integer(0), 
       emp.vario, emp.vario.sd,
       n.bin, PACKAGE="RandomFields", DUP = FALSE)
    emp.vario[is.na(emp.vario) & (centers==0)] <- 0
    return(list(centers=centers, emp.vario=emp.vario, sd=emp.vario.sd,
                n.bin=n.bin))
  } else { ## anisotropic space-time
    ## always transform to full 3 dimensional space-time coordinates
    ## with all angles given. Otherwise there would be too many special
    ## cases to treat in the c program. However, there is some lost
    ## of speed in the calculations...
    stopifnot(is.matrix(new$x))
    if (ncol(new$x)<3)  # not matrix(0,...) here! since x could be a triple
      new$x <- cbind(new$x, matrix(1, nrow=nrow(new$x), ncol=3-ncol(new$x)))
    if (missing(phi) || is.null(phi)) phi <- c(0, 0)
    if (missing(theta) || is.null(theta)) theta <- c(0, 0)
    T <- new$T
    if (is.null(T)) T <- c(1, 1, 1)
    if (missing(deltaT) || is.null(deltaT)) deltaT <- c(0,0)
    if (any(as.integer(deltaT[2] / T[3]) * T[3] != deltaT[2]))
      stop("deltaT not multiple of distance of temporal grid")
    if (any(diff(bin)<=0)) stop("bin must be a strictly increasing sequence")
    stopifnot(0 <= phi[1], 2 * pi > phi[1], 0 <= theta[1], 2 * pi > theta[1],
              phi[2] >= 0,  phi[2] == as.integer(phi[2]), 
              theta[2] >= 0, theta[2] == as.integer(theta[2]))
    stopifnot(all(is.finite(deltaT)), all(deltaT >= 0))
    realdelta <- deltaT[2]
    stepT <-  deltaT[2] / T[3]
    if (stepT != as.integer(stepT)) stop("delta[2] must be a multiple of T[3].")
    deltaT <- c(max(1, stepT),
                min(deltaT[1], T[2] - T[1]) / max(T[3], deltaT[2]))
    n.phi <- max(1, 2 * phi[2])
    n.theta <- max(1, theta[2])
    n.delta <- 1 + deltaT[2]
    totalbins <- n.bins * n.phi * n.theta * n.delta
    emp.vario <- double(totalbins)
    emp.vario.sd <- double(totalbins)
    n.bin <- integer(totalbins)
    phibins <- thetabins <- Tbins <- NULL

    .C("empvarioXT", as.double(new$x), as.double(T), as.integer(new$l),
        as.double(data), as.integer(repet), as.integer(new$grid), 
       as.double(bin), as.integer(n.bins),
       as.double(c(phi[1], phi[2])),
       as.double(c(theta[1], theta[2])),
       as.integer(deltaT),
       ##  input : deltaT[1] max abstand, deltaT[2] : echter gitterabstand,
       ##    c   : delta[1] : index gitterabstand, deltaT[2] : number of bins -1
       ##                   (zero is the additional distance)
       emp.vario, emp.vario.sd,
       n.bin,  PACKAGE="RandomFields", DUP = FALSE)


    ## the results are now reformatted into arrays
    ## the angles are given in clear text    
    NotimeComponent <- (length(seq(T[1], T[2], T[3]))==1 ||
                        missing(deltaT) || any(deltaT==0))  
    if (NotimeComponent) {
      ## no genuine time component
      if (n.theta>1) { ## no real 3 dimensional component
        thetabins <- theta[1] +  0 : (n.theta-1) * pi / n.theta
        if (n.phi>1) {
          n.phi <- n.phi / 2
          phibins <- phi[1] + 0 : (n.phi - 1) * pi / n.phi
          n.bin <- array(n.bin, dim=c(n.bins, n.phi, 2, n.theta))
          n.bin <- n.bin[,,1,] + n.bin[,,2,]
          emp.vario <- array(emp.vario, dim=c(n.bins, n.phi, 2, n.theta))
          emp.vario <- emp.vario[,,1,] + emp.vario[,,2,]
          emp.vario.sd <- array(emp.vario.sd, dim=c(n.bins, n.phi, 2, n.theta))
          emp.vario.sd <- emp.vario.sd[,,1,] + emp.vario.sd[,,2,]
          phibins <- phi[1] + 0 : (n.phi - 1) * 2 * pi / n.phi
        } else {
          n.bin <- matrix(n.bin, ncol=n.bins, nrow=thetabins)
          emp.vario <- matrix(emp.vario, ncol=n.bins, nrow=n.theta)
          emp.vario.sd <- matrix(emp.vario.sd, ncol=n.bins, nrow=n.theta)
        }
      } else { # n.theta==1
        if (n.phi>1) {
          ## phi[2]!=0, theta[2]==0, no time component
          n.phi <- n.phi /2
          half <- totalbins/2
          n.bin <- matrix(n.bin[1:half] + n.bin[(half+1):totalbins], ncol=n.phi)
          emp.vario <-  matrix(emp.vario[1:half] + emp.vario[(half+1):totalbins],
                               ncol=n.phi) 
          emp.vario.sd <-
            matrix(emp.vario.sd[1:half] + emp.vario.sd[(half+1):totalbins],
                               ncol=n.phi)
          phibins <- phi[1] + 0 : (n.phi - 1) * pi / n.phi
        } else {
          ## nothing to do
        }
      }
    } else {
      ## genuine time component
      dims <- c(n.bins, n.phi, n.theta, n.delta)
      dims <- dims[dims!=1]
      n.bin <-  array(n.bin, dim=dims)
      emp.vario <- array(emp.vario, dim=dims) 
      emp.vario.sd <- array(emp.vario.sd, dim=dims)
      Tbins <- (0:deltaT[2]) * realdelta
      if (n.theta>1) thetabins <- theta[1] +  0 : (n.theta-1) * pi / n.theta
      if (n.phi>1) phibins <- phi[1] + 0 : (n.phi - 1) * 2 * pi / n.phi
    }

    idx <- n.bin > 0
    emp.vario[idx] <- emp.vario[idx] / n.bin[idx]
    emp.vario[!idx] <- NaN
    idx <- n.bin > 1


    evsd <- emp.vario.sd[idx] / (n.bin[idx] - 1) -
      n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2
    if (any(evsd) < -1e-14) {
      print(idx)
      print(n.bin[idx] - 1)
      print(emp.vario.sd[idx])
      print(emp.vario.sd[idx] / (n.bin[idx] - 1) -
            n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2)
      warning(paste(evsd))
    }
    evsd[evsd < 0] <- 0
    emp.vario.sd[idx] <- sqrt(evsd)   
    emp.vario.sd[!idx] <- NaN
    
    return(list(centers=centers, emp.vario=emp.vario, sd=emp.vario.sd,
                n.bin=n.bin,
                phi=phibins,
                theta=thetabins,
                T=Tbins
                ))
  }
}

