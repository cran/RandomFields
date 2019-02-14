
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2017 -- 2017 Martin Schlather
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



#   stop("big data sets currently not allowed") 
#    printlevel <- mle.methods <- lsq.methods <- recall <- sdvar <- general <- TRUE




## to do: grid
GetNeighbourhoods <- function(model, Z, X,
                              splitfactor, maxn, split_vec, shared=FALSE,
                              consider.correlations = TRUE) {

  model.nr <- MODEL_AUX
  rfInit(model, Z$coord, dosimulate=FALSE, reg=model.nr)

  
  lc <- nrow(Z$data[[1]]) ## length(Z$coords)
  maxn <- as.integer(maxn)        ## max number of points, including neighbours
  minimum <- as.integer(split_vec[1])## min. number of points in a neighbourhood
  splitn <- as.integer(split_vec[2]) ## number of location when split up
  maximum <- as.integer(split_vec[3])## maximum number of points when still
  ##                                    neighbours of neighbours are included.
  ##                         Note that, mostly, an additional box is included.
  locfactor <- as.integer(splitfactor * 2 + 1) ## total number of
  ##    neighbouring boxes in each direction, including the current box itself
  xdimOZ <- Z$xdimOZ
  tsdim <- Z$tsdim
  newZ <- list()
 
  if (Z$dist.given) { ## to do
   stop("space splitting for 'distances' not programmed yet")
   if (is.vector(dist)) {
      j <- 1:lc
      composite <- list()
      li <- 0
      maxi <-  (splitfactor[2] - 1) / Z$vdim 
        # mini <-  splitfactor[3] / Z$vdim
      while (length(j) >= maxi) {
        distpos <- (j[1]-1) * (lc - j[1] / 2) + 1
          distidx <- distpos : (distpos + (lc-j[1]-1))
        
        locs <- (j[1]+1):lc
        locidx <- which(!is.na(pmatch(locs, j)))
        locidx <- c(locs[locidx], j[1])            
        li <- li + 1
        composites[[li]] <- locidx
        j <- j[is.na(pmatch(j, locidx))]
      }
      if (length(j) > 0) {
        li <- li + 1
        composite[[li]] <- j
      }
      
      ## kann jetzt noch optimiert werden hinsichtlich schaetzqualitaet durch
      ## zusammenfassen
    } else {
    }
   
   ## return(result)
   return(newZ)
  }

  newZ <- list()
  restotal <- sapply(Z$coord, function(x) x$restotal) ## passt das ? 2.2.19
  for (set in 1:lc) {
    ## n <- as.integer(dim(Z$coord[[set]])[2])
    oldlen <- length(newZ)
    nOT <- nT <- restotal
    if (Z$has.time.comp) nOT <- nOT / coord$T[3]
    coord <- Z$coord[[set]]
    
    if (consider.correlations) {
      natsc <- .C(C_MultiDimRange, as.integer(model.nr),
                  as.integer(set),
                  natsc = double(xdimOZ)
                  )$natsc
    } else natsc <- rep(1, xdimOZ)  

    data <- Z$data[[i]]    
    if (coord$grid) {
      stop("space splitting for 'grid' not programmed yet")
      X <- cbind(coord$x, coord$T)
      dim(data) <- c(X[3, ], dim(data)[2])
      pts <- natsc / X[2, ]
      total <- prod(pts)
      factor <-(total / splitn) ^ (1/tsdim)
      pts <- pts / factor
      blocks <- round(X[3,] / pts)
      blocks[blocks < 1] <- 1
      old <- FALSE
      while(prod(ceiling(X[3, ] / blocks)) > maximum) {
        stopifnot(any(!old))
        idx <- blocks == max(blocks[!old])        
        old <- old | idx
        blocks[idx] <- blocks[idx] + 1
      }
      minpts <- trunc(X[3, ] / blocks)
      remaining <- X[3, ] - minpts * blocks
      blocknpts <- (blocks - remaining) * minpts
      listseq <- combi <- list()
      for (i in 1:tsdim) {
        combi[[i]] <- c(-1, 1)
        listseq[[i]] <- 1:blocknpts[i]
      }      
      combi <- if (tsdim > 1) do.call(expand.grid, combi)
               else list(as.matrix(combi))
      for (j in 1:nrow(combi)) {
        L <- list(data)
        idx.combi <- combi[j, ] 
        for (i in 1:tsdim) L[[i+1]] <- idx.combi[i] * listseq[[i]]
        L[[length(L) + 1]] <- TRUE ## vdim + repet
     
        if (all(idx.combi > 0 | remaining[idx.combi < 0] > 0)) {
          newZ[[length(newZ) + 1]] <- Z[[i]]
          pts <- minpts + (idx.combi < 0)
          newZ[[length(newZ)]]$coords$x[3, ] <- pts
          newZ[[length(newZ)]]$data <- matrix(nrow=prod(pts), do.call("[", L))
        }
      }
    } else { # !coord$grid
      maxnOT <- maxn
      nDsplitn <- nT / splitn
      x <- t(coord$x)

      u <- numeric(tsdim)
      for (i in 1:xdimOZ) {
        u[i] <- length(unique(coord$x[i,]))
      }
      Range <- apply(coord, 1, range)
      if (Z$has.time.comp) {
        T <- coord$T
        u <- c(u, T[3])
        Range <- cbind(Range, (1:T[3]) * T[2])
        x <- cbind(x, rep(T[2], prod(T[2:3])))
      }
      
      rd <- apply(Range, 2, diff) 
      len <- pmax(1e-10 * max(rd), rd * (1 + 1e-10))
      units <- pmax(1, len * natsc)
            
      ## * "gerechte" Aufteilung in alle Richtungen waere nDsplitn
      ## * die Richtung in die viele units sind, soll eher aufgespalten werden
      ## * ebenso : wo viele Werte sind eher aufspalten
      blockidx <- (nDsplitn / prod(units * u))^{1/tsdim} * units * u > 0.5
      reddim <- sum(blockidx)
      units <- units[blockidx]
      zaehler <- 1
      blocks <- rep(1, tsdim)
      OK <- integer(1)
      
      repeat {
        blocks[blockidx] <- (nDsplitn / prod(units))^{1/reddim} *
          locfactor * zaehler * units * Z$vdim
        blocks <- as.integer(ceiling(blocks))
        cumblocks <- cumprod(blocks)
        Ccumblocks <- as.integer(c(1, cumblocks))
        cumblocks <- Ccumblocks[-length(Ccumblocks)]
        totblocksOT <- as.integer(cumblocks[xdimOZ])

        if (Z$has.time.comp && blocks[tsdim] > 1) {
          maxnOT <- as.integer(maxn / ceiling(T[3] / blocks[tsdim]))
        }
        
        ## zuordnung der coordinaten_Werte zu den jeweiligen "blocks"
        ## given ist liegend
        coord.idx <- floor((x - Range[1,]) / (len / blocks))
        cumidx <- as.integer(colSums(coord.idx * cumblocks))
      
        elms.in.boxes <- .Call(C_countelements, cumidx, nOT, totblocksOT)
        neighbours <- .Call(C_countneighbours, xdimOZ, blocks, locfactor,
                            Ccumblocks, elms.in.boxes, maxnOT)
      
        ## if there too many points within all the neighbours, then split
        ## into smaller boxes
        zaehler <- zaehler * 2
      
        ## image(neighbours, zlim=c(0:(prod(blocks)-1)))
        if (!is.null(neighbours)) break;
      } # repeat
    
      l <- list()
      l[[1]] <- .Call(C_getelements, cumidx, xdimOZ, nOT, Ccumblocks,
                      elms.in.boxes)
      l1len <- sapply(l[[1]], length)

      if (FALSE) {

	l <- .Call(C_getelements, cumidx, xdimOZ, nOT, Ccumblocks,
		   elms.in.boxes)
      for (idx in l) {
        new.x <- Z[[i]]$coords[idx, ]
         
        if (Z$has.time.comp && blocks[tsdim] > 1) {
          minpts <- trunc(T[3] / blocks[tsdim])
          remaining <- T[3] - minpts * blocks
          blocknpts <- (blocks - remaining) * minpts
          combi <- c(-1, 1)
          listseq <- 1:blocknpts
          new.data <- Z[[i]]$data
          dim(new.data) <- c(nrow(Z[[i]]$data) / T[3], T[3], nrow(new.data))
          for (j in 1:length(combi)) {
            if (combi[j] >0 || remaining[combi[j] < 0] > 0) {
              newZ[[length(newZ) + 1]] <- Z[[i]]
              newZ[[length(newZ)]]$coords <- new.x
              pts <- minpts + (idx.combi < 0)
              newZ[[length(newZ)]]$coords$T[3] <- pts
              newZ[[length(newZ)]]$data <-
                as.matrix(new.data[ ,combi[j] * listseq, ],
                          nrow=length(idx) * pts)
            }
          }
        } else {
          newZ[[length(newZ) + 1]] <- Z[[i]]
          newZ[[length(newZ)]]$coords$x <- new.x
          newZ[[length(newZ)]]$data <- Z[[i]]$data[idx, ]
        }
      } # for idx

      }
      
    }# !coord$grid
    
    if (length(X$x) > 0) {
 ##      l1len <- sapply(l[[1]], length)
      if (X$grid) {
        stop("not programmed yet")
      } else { 
        ## now calculate the boxes for the locations where we will interpolate
        i <- pmax(0, pmin(blocks-1,
                          floor((t(coord$x) - Range[1,]) / (len / blocks))))
	#### lenXXXXX
        dim(i) <- rev(dim(coord$x))
        i <- as.integer(colSums(i * cumblocks))
      
        res.in.boxes <- .Call(C_countelements, i, nrow(coord$x), totblocksOT)
        
        notzeros <- res.in.boxes > 0
        l[[3]] <-
          .Call(C_getelements, i, xdimOZ, as.integer(nrow(coord$x)),
                Ccumblocks, res.in.boxes)[notzeros]
        ## TO DO : idx[[3]] passt nicht, da sowohl fuer Daten
        ##          als auch coordinaten verwendet wird. Bei repet > 1
        ##         ist da ein Problem -- ueberpruefen ob repet=1
        
        
        ll <- .Call(C_getneighbours, xdimOZ, blocks, locfactor, Ccumblocks,
                    neighbours)[notzeros]
        less <-
          sapply(ll, function(x) sum(elms.in.boxes[x]) < minimum) | !shared
        ##                  if !shared then all(less)==TRUE
      
        if (any(less)) {
          not.considered.yet <- sapply(l[[1]], length) > 0   
          newll <- ll
          for (i in which(less)) {
            current <- ll[[i]]
            elements <- sum(elms.in.boxes[current] *
                            (shared | not.considered.yet[current]))# number of pts in a neighbourhood
            while (elements < minimum) {
              new <- unique(unlist(ll[current])) # neighbours of neighbours, but not
              new <- new[which(is.na(pmatch(new, current)))]# neighbours themselves
              nn <- elms.in.boxes[new] * (shared | not.considered.yet[new]) # how many pts are in each of these boxes?
              ordr <- order(nn)
              new <- new[ordr]
              nn <- nn[ordr]
              cs <- elements + cumsum(nn)
              smaller <- sum(cs <= maximum) ## now, check which neighbours of
              ## the neigbours can be included in the list of neighbours of i
              ## to increase the number of points in the kriging neighbourhood
              if (smaller == 0) break; ## none
              if (smaller == length(cs) || cs[smaller] >= minimum ||
                  cs[smaller+1] > maxn) {
                if ( (elements <- cs[length(cs)]) <= maxn ) {            
                  current <- c(current, new)            
                } else {
                  current <- c(current, new[1:smaller])
                  elements <- cs[smaller]
                }
                if (smaller != length(cs)) break
              } else {
                ## smaller < length(cs) && cs[smaller] < minimum &&
                ## cs[smaller+1]<=maxn
                ## i.e., include the next one, but there is no chance to include
                ## more within the rules.
                elements <- cs[smaller+1]
                current <- c(current, new[1:(smaller+1)])
                break;
              }
            }
            current <- current[l1len[current] > 0]
            if (!shared) current <- current[not.considered.yet[current]]
            newll[[i]] <- current
            not.considered.yet[current] <- FALSE                            
          }
          newll <- newll[sapply(newll, length) > 0]
          l[[2]] <- newll
        } else l[[2]] <- ll
      } ## locations to be estimated not on grid
    } ## locations to be estimated
  } ## for
  return(if (shared) l else lapply(l[[2]], function(x) unlist(l[[1]][x])))
}

GetComposites <- function(Z, cliquesize) {
  stopifnot(cliquesize == 2)
  return(Z)
}


BigDataSplit <- function(Z, RFopt) {
  fit <- RFopt$fit
  method <- fit$likelihood   
  if (is.na(method <- pmatch(method, RC_LIKELIHOOD_NAMES)))
    stop("unknown value for 'likelihood'.")
  method <- RC_LIKELIHOOD_NAMES[method] # kein + 1 notwendig

  restotal <- sapply(Z$coord, function(x) x$restotal)
  
  if (method == "full" ||
      (method %in% c("auto", "tesselation") && all(restotal<=fit$max_neighbours)
       )) return(Z)
  if (method == "auto") {
    method <- "tesselation"
    if (RFopt$basic$printlevel>=PL_IMPORTANT)
      message("Too many locations to use standard estimation methods.\n",
              "Hence an approximative methods is used. However, it is *not* ",
              "ensured\n",
              "that they will do well in detecting long memory dependencies.")
  }

  stop("not programmed yet")

  model <- list("RFfctn", Z$model)


  stop("Neigbour kriging/selection not programmed yet")

  if (method == "tesselation") {
    if (any(diff(fit$cliquesize) <= 0)) stop("in case of 'tesselation', 'cliquesize' must contain three increasing numbers")
    return(GetNeighbourhoods(model, Z=Z, 
                             splitfactor=fit$splitfactor_neighbours,
                             maxn=fit$max_neighbours,
                             split_vec=fit$cliquesize,
                             shared=FALSE)
           )
  } else if (method == "composite") {
    if (any(diff(fit$cliquesize) != 0)) stop("in case of 'composite', 'cliquesize' must be a single value")
    return(GetComposites(Z=Z, cliquesize = fit$cliquesize))
  } else stop("unknown 'likelihood' value")  
}
