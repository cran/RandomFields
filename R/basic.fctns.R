
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 -- 2017  Martin Schlather
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




## FUNCTION-STARP***********************************************************************************
## NAME		extract VarNames
## REQUIRE	$model is a formula or a RMmodel
## ENSURE		the return value is a vector with the names of the response variables of the formula 
#			or NULL if no left side given
## SEE		
## AUTHOR		Sebastian Gross 
## DATE		26.08.2011; 2014 Martin Schlather modified
## FUNCTION-END*************************************************************************************
extractVarNames <- function(model) {
  if (missing(model) || length(model) <=2 || !is(model, "formula")) return(NULL)
  tmp <- as.character(model)[2]
  tmp <- strsplit(tmp, "c\\(")[[1]]
  tmp <- paste(tmp, sep="", collapse="")
  tmp <- strsplit(tmp, "\\)")[[1]]
  tmp <- paste(tmp, sep="", collapse="")
  varnames <- strsplit(tmp, ", ")[[1]]

  ## ignore numeric formated varnames
  i<- 1
  while (i <= length(varnames)) {
    if (regexpr("^[[:digit:]]", varnames[i]) > 0) varnames <- varnames[-i]
    else i <- i+ 1
  }
  if (length(varnames) == 0) return (NULL)
  return (varnames)
}


earth_coordinate_names<- function(names) {
  ## Earth coordinates + possibly radius
  n <- substr(tolower(names), 1, 6)
  nc <- nchar(n)
  lon <- lat <- logical(length(n))
  for (i in 1:length(n)) {
    lon[i] <- substr(COORD_NAMES_EARTH[1], 1, nc[i]) == n[i]
    lat[i] <- substr(COORD_NAMES_EARTH[2], 1, nc[i]) == n[i]
  }
  lonORlat <- lon | lat  
  earth <- all(nc[lonORlat] >= 2) && sum(lon==1) && sum(lat == 1)
  
  return(if (length(names)==2 | !earth) earth else         
         if ( (lo <- which(lon)) < (la <- which(lat)))
         which(lonORlat[]))
}

cartesian_coordinate_names <- function(names) {
  n <- substr(tolower(names), 1, 1)
  coords <- COORD_NAMES_CART[c(4, 1:3)] # c("T", "x", "y", "z")
  Txyz <- outer(n, coords, "==")
  cs <- colSums(Txyz)
  if (any(cs > 1) || sum(cs[1:2]) == 0 || any(diff(cs[-1]) > 0))
    return (integer(0))
  Txyz <- Txyz[, c(2:4, 1), drop=FALSE]
  ord <- apply(Txyz, 2, function(x) which(x > 0))
  ord <- order(unlist(ord))
  rs <- which(rowSums(Txyz) > 0)
  return(rs[ord])
}



general_coordinate_names <- function(names) {
  n <- substr(tolower(names), 1, 5)
  return(which(n == "coord"))
}


extractFromNames <- function(varidx, varnames, RFopt, cn) {
  ## this function partially interprets RFoptions$coord$varnames/varidx,
  ## namely when varnames is not NA
  if (length(varnames) == 0) {
    is.data <- varidx
    if (is.na(is.data[2])) is.data[2] <- ncol(data)
    return(is.data[1] : is.data[2])
  } 
  
  if (RFopt$general$vdim_close_together)
    stop("'vdim_close_together' must be FALSE")
  l <- list()
  vdim <- length(varnames)
  for (v in 1:vdim)
    l[[v]] <- substring(cn, 1, nchar(varnames[v])) == varnames[v]
  repet <- sapply(l, sum)
  if (repet[1] == 0) return(NULL)
  if (any(repet != repet[1]))
    stop("detected repetitions are not equal for all components")
  m <- matrix(unlist(l), ncol=vdim)
  if (any(rowSums(m) > 1))
    stop("names of multivariate components not unique")
  return(as.vector(t(apply(m, 2, which))))
}


data.columns <- function(data, model=NULL, xdim=0, force=FALSE, halt=TRUE,
			 RFopt=RFoptions(), vdim=0) {
  ## this function interprets RFoptions$coord$varnames/varidx or model
  ## also if varnames is NA
  ## gives the columns in the original data

  PL <-  as.integer(RFopt$basic$printlevel)
  info <- RFopt$coords

  is.data <- is.x <- NULL
  DATA <- if (is(data, "RFsp")) data@data else data    
  if (!missing(model) && !is.null(model )) {
    ## ehemals selectDataAccordingFormula & selectAccordingFormula
    ##
    ## a vector of column indices  or N ULL indicating 'all' in case
    ## of model is given
    
    varnames <- extractVarNames(model) 
    ## are there any variablenames given
    if (length(varnames) > 0) { ## if == 0 continue as if model is N ULL
      if (is(data, "RFsp")) 
	cleanNames <- (if (data@.RFparams$n == 1) colnames(DATA)
		       else sapply(colnames(DATA),
                                   ## AUTHOR: Sebastian Gross, 26.08.2011
                                   FUN= function(x) strsplit(x, "\\.n[[:digit:]]+$")[[1]][1]
                                   ))
      else {
	cleanNames <- try(colnames(data))
	if (is(cleanNames, "try-error"))
	  stop("could not detect colnames of data")
      }
      is.data <- match(varnames, cleanNames)
      if (any(is.na(is.data)))
	stop("response variable names could not be found in colnames ",
	     "of the data")
      
      if (PL > 0)
	message("formula yields the following variables for the data: ",
		paste(cleanNames[is.data], sep=", "), "\n")
      names(is.data) <- varnames
    }
  }

  if (is.null(is.data)) 
    is.data <-
	 if (ncol(DATA) == 1) 1L
	 else if (length(vdim) == 1 && vdim > 0 && ncol(DATA) == vdim)
	   1:vdim

  cn <- colnames(DATA)  
  if (isRFsp <- is(data, "RFsp")) {
    xdim <- if (is(data, "SpatialPointsDataFrame") ||
		is(data, "RFpointsDataFrame")) ncol(data@coords)
            else if (is(data, "SpatialGridDataFrame") ||
		     is(data, "RFgridDataFrame")) length(data@grid@cells.dim)
	    else stop("unknown class for 'data'")  
  } else {
    if (length(xdim) == 0) xdim <- 0
    if (xdim>0 && xdim >= ncol(DATA)) stop("not enough columns in 'data'.")

    if (length(info$coordnames) > 0 || !is.na(info$coordidx[2])) {
      is.x <- extractFromNames(info$coordidx, info$coordnames, RFopt, cn)
      if (xdim > 0 && xdim != length(is.x))
	stop("expected dimension of coordinates does not match the found coordinates")
    }

    ## ehemals StandardizeData:    
    if (is.null(is.x) && !is.null(cn)) {
      if (length(is.x <- earth_coordinate_names(cn)) == 2) {
	txt <- paste("keywords for earth coordinates (", sep="",
		     paste(cn[is.x], collapse=","), ")")
      } else if (length(is.x <- cartesian_coordinate_names(cn)) > 0) {
	txt <- paste("keywords for cartesian coordinates (", sep="",
		     paste(cn[is.x], collapse=","), ")")
      } else if  (length(is.x <- general_coordinate_names(cn)) > 0) {
	txt <- paste("standard keywords for cartesian coordinates (",
		     paste(cn[is.x], collapse=","), ")", sep="")
      }
      if (length(is.x) > 0 && PL > 0)
	message("columns ", paste(is.x, collapse=","),
		" are assumed to define the coordinates as ", txt,
		" are recognized.\n")
    }
  }

  ## ehemals data.columns
  if (length(is.data) == 0 &&
      (length(info$coordnames) > 0 || !is.na(info$coordidx[2]))) {
    is.data <- extractFromNames(info$varidx, info$varnames, RFopt, cn)
  }

  if (length(is.data) == 0  && !is.null(cn)) {
    is.data <- which(tolower(substr(cn, 1, 4)) == "data" |
		     tolower(substr(cn, 1, 4)) == "value" |
		     tolower(substr(cn, 1, 8)) == "variable")
    if (length(is.data) > 0 && is.data[1] > 1)
      is.data <- is.data[1] : ncol(DATA)  #coord am Anfang
    ##     dann wird Rest alles als data angenommen, egal welcher Name
  }

  if (length(is.data)==0 && length(is.x)==0) {
    if (is.null(cn)) {
      if (!force) {
	if (halt)
	  stop(if (is.null(cn)) "colnames of data argument must contain"
	       else 'no colname starts with', ' "data" or "variable"')
	else return(list(is.data=TRUE, is.x=NULL));
      }
      is.data <- (xdim+1):ncol(DATA)
    } else stop("data name(s) could not be found in the column names")
  }

  if (length(is.x) == 0 && !isRFsp) {
    is.x <- (1:ncol(DATA))[-is.data]
    if (xdim > 0) {
      if (length(is.x) < xdim)
        stop("not enough columns for coordinates found ")
      if (xdim < length(is.x) && PL >= PL_SUBIMPORTANT)
        message("column(s) ",
                paste("'", is.x[-1:-xdim], "'", , sep="", collapse=", "),
		" unused.\n")
      is.x <- is.x[1:xdim]
    }
  }
  
  if (length(is.data) == 0) {  #  (all(is.na(info$varnames))) 
    is.data <- (1:ncol(DATA))[-is.x]
    if (length(is.data) == 0) stop("no columns for data found")
  } else {
    if (any(is.x %in% is.data)) stop("column names and data names overlap.")
    if ( (isRFsp && length(is.data) < ncol(DATA)) ||
	 (!isRFsp && length(is.x) + length(is.data) < ncol(DATA)) ) {
      if (PL >= PL_SUBIMPORTANT)
	message("column(s) ",
                paste((1:ncol(DATA))[-c(is.x, is.data)], collapse=", "),
		" unused.\n")
    }
  }

  if (length(cn) > 0) {
    names(is.data) <- cn[is.data]
    if (length(is.x) > 0) names(is.x) <- cn[is.x]
  } 
 
  return(list(is.data=is.data, is.x=is.x))
}


    
SystemCoordNames <- function(locinfo, RFopt = RFopt) {
  ## this function tries to combine all information available on
  ## coordinate names and variable names und returns the names if
  ## ensured that the names are the correct one.
  has.time.comp <- locinfo$has.time.comp
  tsdim <- locinfo$spatialdim + locinfo$has.time.comp
  
  system <- RFopt$coords$coord_system
  if (system == "earth") {
    coordnames <- if (tsdim == 4) COORD_NAMES_EARTH
		  else if (tsdim == 2) COORD_NAMES_EARTH[1:2]
		  else if (tsdim == 3) c(COORD_NAMES_EARTH[1:2], "HeightOrTime")
  } else if (system == "cartesian" && tsdim <= 4) {
    ##coords <- COORD_NAMES_CART[1:tsdim]
    if (has.time.comp) coordnames[tsdim] <- COORD_NAMES_CART[3 + 1]
  } else {
    coordnames <- paste(COORD_NAMES_GENERAL[1], 1:tsdim, sep="")
    if (has.time.comp) coordnames[tsdim] <- COORD_NAMES_GENERAL[2]
  } 
  
  return(coordnames)
}



search.model.name <- function(cov, name, level) {
  if (length(name) == 0 || length(cov) ==0) return(cov);
  if (!is.na(pmatch(name[1], cov))) return(search.model.name(cov, name[-1], 1))

  for (i in 1:length(cov$submodels)) {
    found <- search.model.name(cov$submodels[[i]], name, 1)
    if (!is.null(found)) return(found)      
  }
  found <- search.model.name(cov$internal, name, 1)
  if (!is.null(found)) return(found)
  if (level == 0) stop("model name not found")
  return(NULL)
}






vectordist <- function(x, diag=FALSE) {
  storage.mode(x) <- "double"
  res <- .Call(C_vectordist, t(x), diag)
  dimnames(res) <- list(dimnames(x)[[2]], NULL)
  return(t(res));
}



my.legend <- function(lu.x, lu.y, zlim, col, cex=1, ...) {
  ## uses already the legend code of R-1.3.0
  cn <- length(col)
  if (cn < 43) {
    col <- rep(col, each=ceiling(43 / cn))
    cn <- length(col)
  }
  filler <- vector("character", length=(cn-3)/2)
  legend(lu.x, lu.y, y.intersp=0.03, x.intersp=0.1, 
         legend=c(format(zlim[2], dig=2), filler,
             format(mean(zlim), dig=2), filler,
             format(zlim[1], dig=2)),
         lty=1, col=rev(col),cex=cex, ...)
}

