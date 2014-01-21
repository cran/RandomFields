RMcoord <- function(C0, coord, dist)
{
	cl <- match.call()

	submodels <- list() 
	if(hasArg(C0)) submodels[['C0']] <- C0

	par.model <- list() 
	if(hasArg(coord)) par.model[['coord']] <- coord
	if(hasArg(dist)) par.model[['dist']] <- dist

	par.general <- list()
	par.general[['var']] <- ZF_DEFAULT_STRING
	par.general[['scale']] <- ZF_DEFAULT_STRING
	par.general[['Aniso']] <- ZF_DEFAULT_STRING
	par.general[['proj']]  <- ZF_DEFAULT_STRING

	model <- new(ZF_MODEL, call = cl, name = ZF_COORD,
	submodels = submodels,
	par.model = par.model, par.general = par.general)

	return(model) 
} 

RMcoord <- new('RMmodelgenerator',
               .Data = RMcoord,
               type = ZF_TYPE[OtherType + 1],
               domain = ZF_DOMAIN[PREVMODELD + 1],
               isotropy = ZF_ISOTROPY[CARTESIAN_COORD + 1],
               operator = TRUE,
               normalmix = TRUE,
               finiterange = TRUE,
               maxdim = Inf,
               vdim = -2
               )


internalRMmixed <- function(X, beta, cov, coord, dist, element)
{
	cl <- match.call()

	submodels <- list() 

	par.model <- list() 
	if(hasArg(element)) par.model[['element']] <- element
	if(hasArg(X)) par.model[['X']] <- X
	if(hasArg(beta)) par.model[['beta']] <- beta
	if(hasArg(coord)) par.model[['coord']] <- coord
	if(hasArg(dist)) par.model[['dist']] <- dist
	if(hasArg(cov)) submodels[['cov']] <- cov #par.model[['cov']] <- cov

	par.general <- list()
	par.general[['var']] <- ZF_DEFAULT_STRING
	par.general[['scale']] <- ZF_DEFAULT_STRING
	par.general[['Aniso']] <- ZF_DEFAULT_STRING
	par.general[['proj']]  <- ZF_DEFAULT_STRING

	model <- new(ZF_MODEL, call = cl, name = ZF_MIXED[1],
	submodels = submodels,
	par.model = par.model, par.general = par.general)

	return(model) 
} 

internalRMmixed <- new('RMmodelgenerator',
                       .Data = internalRMmixed,
                       type = ZF_TYPE[OtherType + 1],
                       domain = ZF_DOMAIN[PREVMODELD + 1],
                       isotropy = ZF_ISOTROPY[CARTESIAN_COORD + 1],
                       operator = TRUE,
                       normalmix = TRUE,
                       finiterange = TRUE,
                       maxdim = Inf,
                       vdim = -2
                       )


RRdistr <- function(fct, nrow, ncol, envir) {
  if (!missing(fct)) {
    u <- try(isModel(fct), silent=TRUE)
    if (is.logical(u) && u)
      return(fct)
  }

  cl <- match.call()
  
  par.general <- submodels <- par.model <- list()

 # Print(substitute(fct), rate);  xxxx
  
  ll <- as.list(substitute(fct))
  name <- as.character(ll[[1]])
  ll <- ll[-1]
  if (length(ll) > 0) {
    par.names <- names(ll)
    if (length(par.names) == 0 || any(par.names == ""))
      stop("In distribution families, all parameters must be named.")
   # print(par.names)
    num <- sapply(ll, function(x) {
      is.numeric(x) || is.symbol(x) || {u <- try(eval(x), silent=TRUE);
                                        !(class(u)=="try-error")}
    })
    if (!all(num)) {
      subs <- ll[!num]
      n <- names(subs)
      for (i in 1:length(subs)) {
      #  Print(subs[[i]], n[[i]], deparse(subs[[i]]), substr(deparse(subs[[i]]), 1, 1))
   #     ffff
        if (!is.language(subs[[i]]))
          stop("type of parameter (function, constant) cannot be determined")
        par.model[[n[i]]] <-
          if (substr(deparse(subs[[i]]), 1, 1)=='R') eval(subs[[i]]) else
              do.call(ZF_DISTR[1], list(subs[[i]]))
      }
   # Print(submodels); xxx
    }
    if (any(num)) {
      ll <- ll[num]
      n <- names(ll)
      for (i in 1:length(ll)) par.model[n[i]] <- eval(ll[[i]])
    }
    #Print(submodels, par.model, ll, num, eval(par.model[[1]]))
  }
  
  var <- c('x', 'q', 'p', 'n')
  fctnames <-  c('d', 'p', 'q', 'r')
  for (ii in 1:length(fctnames)) {
    i <- fctnames[ii]
    par.model[[paste(i, "distr", sep="")]] <-
      eval(parse(text=paste("quote(", i, name, "(",
                   var[ii], "=", var[ii],
                   if (length(ll) > 0 && length(par.names)>0)
                      paste(",",
                            paste(par.names, "=", par.names, collapse=", ")),
                   "))", sep="")))
  }
  
  if (hasArg(ncol)) par.model[['ncol']] <- ncol
  if (hasArg(nrow)) par.model[['nrow']] <- nrow
  par.model[['envir']] <-  if (hasArg(envir)) envir else new.env()
    
  model <- new(ZF_MODEL, call = cl, name = ZF_DISTR[1], submodels = submodels, 
               par.model = par.model, par.general = par.general)
  return(model) 
}

RRdistr <- new('RMmodelgenerator',
	.Data = RRdistr,
	type = ZF_TYPE[.RandomType + 1],
        domain = ZF_DOMAIN[PREVMODELD + 1],
	isotropy = ZF_ISOTROPY[PREVMODELI + 1],
	operator = FALSE,
	normalmix = FALSE,
	finiterange = FALSE,
	maxdim = Inf,
	vdim = 1
	)




GetSymbols <- function(ll) {
  idx <- sapply(ll, function(x) is.symbol(x) || !is.language(x))
  symbols <- as.character(ll[idx])
  #Print(ll, idx)
  if (!all(idx)) 
    for (i in which(!idx)) {
      symbols <- c(symbols, GetSymbols(as.list(ll[[i]])))
    }
  return(symbols)
}
       

RMuser <- function(type, domain, isotropy, vdim, beta,
                   fctn, fst, snd, envir, var, scale, Aniso, proj) {
	cl <- match.call()
	submodels <- par.general <- par.model <- list() 
	
	if (hasArg(type)) {
	  if (is.numeric(type)) par.model[['type']] <- type
          else if (is.character(type))
                   par.model[['type']] <- pmatch(type, ZF_TYPE) - 1
	  else stop("undefined value for 'type'")
	}
	if (hasArg(domain)) {
	  if (is.numeric(domain)) par.model[['domain']] <- domain
	  else if (is.character(domain))
                   par.model[['domain']] <- pmatch(domain, ZF_DOMAIN) - 1
	  else stop("wrong value for 'domain'")
	}
	if (hasArg(isotropy)) {
	  if (is.numeric(isotropy)) par.model[['isotropy']] <- isotropy
	  else if (is.character(isotropy))
                   par.model[['isotropy']] <- pmatch(isotropy, ZF_ISOTROPY) - 1
	  else stop("wrong value for 'isotropy'")
	}

        
   #     Print(par.model, type, domain, isotropy, ZF_TYPE, ZF_DOMAIN, ZF_ISOTROPY); xxx        
	if (hasArg(vdim)) {
	  if (is.numeric(vdim)) par.model[['vdim']] <- vdim
	  else stop("wrong value for 'vdim'")
	}
	if (hasArg(beta)) {
	  if (is.numeric(beta) || is.language(beta) || 
	      is.list(beta))
	     par.model[['beta']] <- beta
	  else if (substr(deparse(substitute(beta)), 1, 1)=='R')
	    submodels[['beta']] <- beta
	  else submodels[['beta']] <- RRdistr(beta)
	}
        par.names <- c("x", "y", "z", "T")
	if (hasArg(fctn)) {
          f <- substitute(fctn)
          par <- par.names %in% GetSymbols(as.list(as.list(f)[-1]))
          par.model[['fctn']] <- f
	} else stop("'fctn' not given")
       
	if (hasArg(fst)) {
          f <- substitute(fst)
          if (any(xor(par, par.names %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'fst' do not match the ones in 'fctn'")
          par.model[['fst']] <- f
	} else if (hasArg(snd)) stop("'fst' not given")
        
	if (hasArg(snd)) {         
	  f <- substitute(snd)          
          if (any(xor(par, par.names %in% GetSymbols(as.list(as.list(f)[-1])))))
            stop("the variables in 'snd' do not match the ones in 'fctn'")
	  par.model[['snd']] <- f
	}

        ##      Print(par.names, par);        xxx
        par.model[['variab.names']] <- which(par)        
	par.model[['envir']] <- if (hasArg(envir)) envir else new.env()
	par.general[['var']] <-if (hasArg(var)) var else ZF_DEFAULT_STRING
	par.general[['scale']] <-if (hasArg(scale)) scale else ZF_DEFAULT_STRING
	par.general[['Aniso']] <-if (hasArg(Aniso)) Aniso else ZF_DEFAULT_STRING
	par.general[['proj']] <-if (hasArg(proj)) proj else ZF_DEFAULT_STRING

	model <- new(ZF_MODEL, call = cl, name = ZF_USER[1], 
			submodels = submodels, 
			par.model = par.model, par.general = par.general)
	return(model) 
 } 

RMuser <- new('RMmodelgenerator',
              .Data = RMuser,
              type = ZF_TYPE[PosDefType + 1],
              domain = ZF_DOMAIN[PREVMODELD + 1],
              isotropy = ZF_ISOTROPY[PREVMODELI + 1],
              operator = FALSE,
              normalmix = FALSE,
              finiterange = TRUE,
              maxdim = Inf,
              vdim = -1
              )


f <- substitute(exp(-x-y+z(), zz=a))
