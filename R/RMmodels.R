
## This file is created automatically by 'rfGenerateModels'.


RMtrend <- function(mean, plane, polydeg, polycoeff, arbitraryfct, fctcoeff) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(mean) && !is.null(subst <- substitute(mean))) {
    u <- try(is.numeric(mean) || is.logical(mean) || is.language(mean)
	 || is.list(mean) || is(mean, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mean']] <- mean
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mean']] <- mean
    else  stop('random parameter not allowed')
  }
  if (hasArg(plane) && !is.null(subst <- substitute(plane))) {
    u <- try(is.numeric(plane) || is.logical(plane) || is.language(plane)
	 || is.list(plane) || is(plane, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['plane']] <- plane
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['plane']] <- plane
    else  stop('random parameter not allowed')
  }
  if (hasArg(polydeg) && !is.null(subst <- substitute(polydeg))) {
    u <- try(is.numeric(polydeg) || is.logical(polydeg) || is.language(polydeg)
	 || is.list(polydeg) || is(polydeg, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['polydeg']] <- polydeg
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['polydeg']] <- polydeg
    else  stop('random parameter not allowed')
  }
  if (hasArg(polycoeff) && !is.null(subst <- substitute(polycoeff))) {
    u <- try(is.numeric(polycoeff) || is.logical(polycoeff) || is.language(polycoeff)
	 || is.list(polycoeff) || is(polycoeff, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['polycoeff']] <- polycoeff
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['polycoeff']] <- polycoeff
    else  stop('random parameter not allowed')
  }
  if (hasArg(arbitraryfct) && !is.null(subst <- substitute(arbitraryfct))) {
    u <- try(is.numeric(arbitraryfct) || is.logical(arbitraryfct) || is.language(arbitraryfct)
	 || is.list(arbitraryfct) || is(arbitraryfct, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['arbitraryfct']] <- arbitraryfct
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['arbitraryfct']] <- arbitraryfct
    else  stop('random parameter not allowed')
  }
  if (hasArg(fctcoeff) && !is.null(subst <- substitute(fctcoeff))) {
    u <- try(is.numeric(fctcoeff) || is.logical(fctcoeff) || is.language(fctcoeff)
	 || is.list(fctcoeff) || is(fctcoeff, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['fctcoeff']] <- fctcoeff
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['fctcoeff']] <- fctcoeff
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RMtrend', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtrend <- new('RMmodelgenerator',
	.Data = RMtrend,
	type = 'trend',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = -1
	)



RMplus <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(C0)) submodels[['C0']] <- C0
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMplus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMplus <- new('RMmodelgenerator',
	.Data = RMplus,
	type = 'undefined',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMmult <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(C0)) submodels[['C0']] <- C0
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMmult', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmult <- new('RMmodelgenerator',
	.Data = RMmult,
	type = 'tail correlation function',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMS  <- function(phi, var, scale, Aniso, proj, anisoT) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['var']] <- var
    else par.model[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['scale']] <- scale
    else par.model[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(anisoT) && !is.null(subst <- substitute(anisoT))) {
    u <- try(is.numeric(anisoT) || is.logical(anisoT) || is.language(anisoT)
	 || is.list(anisoT) || is(anisoT, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['anisoT']] <- anisoT
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['anisoT']] <- anisoT
    else par.model[['anisoT']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['Aniso']] <- Aniso
    else par.model[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['proj']] <- proj
    else par.model[['proj']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMS', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMS <- new('RMmodelgenerator',
	.Data = RMS,
	type = 'undefined',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMave <- function(phi, A, z, spacetime, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(A) && !is.null(subst <- substitute(A))) {
    u <- try(is.numeric(A) || is.logical(A) || is.language(A)
	 || is.list(A) || is(A, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['A']] <- A
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['A']] <- A
    else par.model[['A']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(z) && !is.null(subst <- substitute(z))) {
    u <- try(is.numeric(z) || is.logical(z) || is.language(z)
	 || is.list(z) || is(z, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['z']] <- z
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['z']] <- z
    else par.model[['z']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(spacetime) && !is.null(subst <- substitute(spacetime))) {
    u <- try(is.numeric(spacetime) || is.logical(spacetime) || is.language(spacetime)
	 || is.list(spacetime) || is(spacetime, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['spacetime']] <- spacetime
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['spacetime']] <- spacetime
    else par.model[['spacetime']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMave', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMave <- new('RMmodelgenerator',
	.Data = RMave,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 1
	)



RMbcw <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(beta) && !is.null(subst <- substitute(beta))) {
    u <- try(is.numeric(beta) || is.logical(beta) || is.language(beta)
	 || is.list(beta) || is(beta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['beta']] <- beta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['beta']] <- beta
    else par.model[['beta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbcw', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbcw <- new('RMmodelgenerator',
	.Data = RMbcw,
	type = 'undefined',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMbessel <- function(nu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbessel', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbessel <- new('RMmodelgenerator',
	.Data = RMbessel,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMbigneiting <- function(kappa, mu, s, sred12, gamma, cdiag, rhored, c, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(kappa) && !is.null(subst <- substitute(kappa))) {
    u <- try(is.numeric(kappa) || is.logical(kappa) || is.language(kappa)
	 || is.list(kappa) || is(kappa, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['kappa']] <- kappa
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['kappa']] <- kappa
    else par.model[['kappa']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else par.model[['mu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else par.model[['s']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(sred12) && !is.null(subst <- substitute(sred12))) {
    u <- try(is.numeric(sred12) || is.logical(sred12) || is.language(sred12)
	 || is.list(sred12) || is(sred12, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['sred12']] <- sred12
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['sred12']] <- sred12
    else par.model[['sred12']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(gamma) && !is.null(subst <- substitute(gamma))) {
    u <- try(is.numeric(gamma) || is.logical(gamma) || is.language(gamma)
	 || is.list(gamma) || is(gamma, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['gamma']] <- gamma
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['gamma']] <- gamma
    else par.model[['gamma']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(cdiag) && !is.null(subst <- substitute(cdiag))) {
    u <- try(is.numeric(cdiag) || is.logical(cdiag) || is.language(cdiag)
	 || is.list(cdiag) || is(cdiag, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['cdiag']] <- cdiag
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['cdiag']] <- cdiag
    else par.model[['cdiag']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(rhored) && !is.null(subst <- substitute(rhored))) {
    u <- try(is.numeric(rhored) || is.logical(rhored) || is.language(rhored)
	 || is.list(rhored) || is(rhored, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['rhored']] <- rhored
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['rhored']] <- rhored
    else par.model[['rhored']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(c) && !is.null(subst <- substitute(c))) {
    u <- try(is.numeric(c) || is.logical(c) || is.language(c)
	 || is.list(c) || is(c, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['c']] <- c
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['c']] <- c
    else par.model[['c']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbigneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbigneiting <- new('RMmodelgenerator',
	.Data = RMbigneiting,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 2
	)



RMbernoulli <- function(phi, threshold, correlation, centred, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(threshold) && !is.null(subst <- substitute(threshold))) {
    u <- try(is.numeric(threshold) || is.logical(threshold) || is.language(threshold)
	 || is.list(threshold) || is(threshold, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['threshold']] <- threshold
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['threshold']] <- threshold
    else par.model[['threshold']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(correlation) && !is.null(subst <- substitute(correlation))) {
    u <- try(is.numeric(correlation) || is.logical(correlation) || is.language(correlation)
	 || is.list(correlation) || is(correlation, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['correlation']] <- correlation
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['correlation']] <- correlation
    else par.model[['correlation']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(centred) && !is.null(subst <- substitute(centred))) {
    u <- try(is.numeric(centred) || is.logical(centred) || is.language(centred)
	 || is.list(centred) || is(centred, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['centred']] <- centred
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['centred']] <- centred
    else par.model[['centred']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbernoulli', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbernoulli <- new('RMmodelgenerator',
	.Data = RMbernoulli,
	type = 'tail correlation function',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbiwm <- function(nudiag, nured12, nu, s, cdiag, rhored, c, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nudiag) && !is.null(subst <- substitute(nudiag))) {
    u <- try(is.numeric(nudiag) || is.logical(nudiag) || is.language(nudiag)
	 || is.list(nudiag) || is(nudiag, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nudiag']] <- nudiag
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nudiag']] <- nudiag
    else par.model[['nudiag']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(nured12) && !is.null(subst <- substitute(nured12))) {
    u <- try(is.numeric(nured12) || is.logical(nured12) || is.language(nured12)
	 || is.list(nured12) || is(nured12, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nured12']] <- nured12
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nured12']] <- nured12
    else par.model[['nured12']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else par.model[['s']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(cdiag) && !is.null(subst <- substitute(cdiag))) {
    u <- try(is.numeric(cdiag) || is.logical(cdiag) || is.language(cdiag)
	 || is.list(cdiag) || is(cdiag, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['cdiag']] <- cdiag
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['cdiag']] <- cdiag
    else par.model[['cdiag']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(rhored) && !is.null(subst <- substitute(rhored))) {
    u <- try(is.numeric(rhored) || is.logical(rhored) || is.language(rhored)
	 || is.list(rhored) || is(rhored, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['rhored']] <- rhored
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['rhored']] <- rhored
    else par.model[['rhored']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(c) && !is.null(subst <- substitute(c))) {
    u <- try(is.numeric(c) || is.logical(c) || is.language(c)
	 || is.list(c) || is(c, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['c']] <- c
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['c']] <- c
    else par.model[['c']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(notinvnu) && !is.null(subst <- substitute(notinvnu))) {
    u <- try(is.numeric(notinvnu) || is.logical(notinvnu) || is.language(notinvnu)
	 || is.list(notinvnu) || is(notinvnu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['notinvnu']] <- notinvnu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['notinvnu']] <- notinvnu
    else par.model[['notinvnu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbiwm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbiwm <- new('RMmodelgenerator',
	.Data = RMbiwm,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 2
	)



RMbrownresnick <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbrownresnick', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbrownresnick <- new('RMmodelgenerator',
	.Data = RMbrownresnick,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbr2bg <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbr2bg', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbr2bg <- new('RMmodelgenerator',
	.Data = RMbr2bg,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbr2eg <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMbr2eg', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbr2eg <- new('RMmodelgenerator',
	.Data = RMbr2eg,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMcauchy <- function(gamma, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(gamma) && !is.null(subst <- substitute(gamma))) {
    u <- try(is.numeric(gamma) || is.logical(gamma) || is.language(gamma)
	 || is.list(gamma) || is(gamma, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['gamma']] <- gamma
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['gamma']] <- gamma
    else par.model[['gamma']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMcauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcauchy <- new('RMmodelgenerator',
	.Data = RMcauchy,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMcircular <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMcircular', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcircular <- new('RMmodelgenerator',
	.Data = RMcircular,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'Gneiting-Schaback class',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMconstant <- function(M, vdim, element, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(M) && !is.null(subst <- substitute(M))) {
    u <- try(is.numeric(M) || is.logical(M) || is.language(M)
	 || is.list(M) || is(M, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['M']] <- M
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['M']] <- M
    else par.model[['M']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(vdim) && !is.null(subst <- substitute(vdim))) {
    u <- try(is.numeric(vdim) || is.logical(vdim) || is.language(vdim)
	 || is.list(vdim) || is(vdim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['vdim']] <- vdim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['vdim']] <- vdim
    else par.model[['vdim']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(element) && !is.null(subst <- substitute(element))) {
    u <- try(is.numeric(element) || is.logical(element) || is.language(element)
	 || is.list(element) || is(element, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['element']] <- element
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['element']] <- element
    else par.model[['element']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMconstant', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMconstant <- new('RMmodelgenerator',
	.Data = RMconstant,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'completely monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = -1
	)



RMcoxisham <- function(phi, mu, D, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else par.model[['mu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(D) && !is.null(subst <- substitute(D))) {
    u <- try(is.numeric(D) || is.logical(D) || is.language(D)
	 || is.list(D) || is(D, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['D']] <- D
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['D']] <- D
    else par.model[['D']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(beta) && !is.null(subst <- substitute(beta))) {
    u <- try(is.numeric(beta) || is.logical(beta) || is.language(beta)
	 || is.list(beta) || is(beta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['beta']] <- beta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['beta']] <- beta
    else par.model[['beta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMcoxisham', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcoxisham <- new('RMmodelgenerator',
	.Data = RMcoxisham,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'zero-space-isotropic',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMcubic <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMcubic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcubic <- new('RMmodelgenerator',
	.Data = RMcubic,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMcurlfree <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMcurlfree', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcurlfree <- new('RMmodelgenerator',
	.Data = RMcurlfree,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMcutoff <- function(phi, diameter, a, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(diameter) && !is.null(subst <- substitute(diameter))) {
    u <- try(is.numeric(diameter) || is.logical(diameter) || is.language(diameter)
	 || is.list(diameter) || is(diameter, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['diameter']] <- diameter
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['diameter']] <- diameter
    else par.model[['diameter']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(a) && !is.null(subst <- substitute(a))) {
    u <- try(is.numeric(a) || is.logical(a) || is.language(a)
	 || is.list(a) || is(a, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['a']] <- a
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['a']] <- a
    else par.model[['a']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMcutoff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcutoff <- new('RMmodelgenerator',
	.Data = RMcutoff,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = TRUE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RMdagum <- function(beta, gamma, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(beta) && !is.null(subst <- substitute(beta))) {
    u <- try(is.numeric(beta) || is.logical(beta) || is.language(beta)
	 || is.list(beta) || is(beta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['beta']] <- beta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['beta']] <- beta
    else par.model[['beta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(gamma) && !is.null(subst <- substitute(gamma))) {
    u <- try(is.numeric(gamma) || is.logical(gamma) || is.language(gamma)
	 || is.list(gamma) || is(gamma, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['gamma']] <- gamma
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['gamma']] <- gamma
    else par.model[['gamma']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMdagum', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdagum <- new('RMmodelgenerator',
	.Data = RMdagum,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMdampedcos <- function(lambda, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(lambda) && !is.null(subst <- substitute(lambda))) {
    u <- try(is.numeric(lambda) || is.logical(lambda) || is.language(lambda)
	 || is.list(lambda) || is(lambda, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['lambda']] <- lambda
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['lambda']] <- lambda
    else par.model[['lambda']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMdampedcos', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdampedcos <- new('RMmodelgenerator',
	.Data = RMdampedcos,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMdewijsian <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMdewijsian', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdewijsian <- new('RMmodelgenerator',
	.Data = RMdewijsian,
	type = 'negative definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMdivfree <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMdivfree', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdivfree <- new('RMmodelgenerator',
	.Data = RMdivfree,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMepscauchy <- function(alpha, beta, eps, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(beta) && !is.null(subst <- substitute(beta))) {
    u <- try(is.numeric(beta) || is.logical(beta) || is.language(beta)
	 || is.list(beta) || is(beta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['beta']] <- beta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['beta']] <- beta
    else par.model[['beta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(eps) && !is.null(subst <- substitute(eps))) {
    u <- try(is.numeric(eps) || is.logical(eps) || is.language(eps)
	 || is.list(eps) || is(eps, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['eps']] <- eps
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['eps']] <- eps
    else par.model[['eps']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMepscauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMepscauchy <- new('RMmodelgenerator',
	.Data = RMepscauchy,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMexp <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMexp <- new('RMmodelgenerator',
	.Data = RMexp,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'completely monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMexponential <- function(phi, n, standardised, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(n) && !is.null(subst <- substitute(n))) {
    u <- try(is.numeric(n) || is.logical(n) || is.language(n)
	 || is.list(n) || is(n, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['n']] <- n
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['n']] <- n
    else par.model[['n']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(standardised) && !is.null(subst <- substitute(standardised))) {
    u <- try(is.numeric(standardised) || is.logical(standardised) || is.language(standardised)
	 || is.list(standardised) || is(standardised, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['standardised']] <- standardised
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['standardised']] <- standardised
    else par.model[['standardised']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMexponential', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMexponential <- new('RMmodelgenerator',
	.Data = RMexponential,
	type = 'positive definite',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMschlather <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMschlather', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMschlather <- new('RMmodelgenerator',
	.Data = RMschlather,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMfractdiff <- function(a, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(a) && !is.null(subst <- substitute(a))) {
    u <- try(is.numeric(a) || is.logical(a) || is.language(a)
	 || is.list(a) || is(a, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['a']] <- a
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['a']] <- a
    else par.model[['a']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMfractdiff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfractdiff <- new('RMmodelgenerator',
	.Data = RMfractdiff,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



RMfbm <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMfbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfbm <- new('RMmodelgenerator',
	.Data = RMfbm,
	type = 'negative definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'Bernstein',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMfractgauss <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMfractgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfractgauss <- new('RMmodelgenerator',
	.Data = RMfractgauss,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



RMgauss <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgauss <- new('RMmodelgenerator',
	.Data = RMgauss,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgenfbm <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(beta) && !is.null(subst <- substitute(beta))) {
    u <- try(is.numeric(beta) || is.logical(beta) || is.language(beta)
	 || is.list(beta) || is(beta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['beta']] <- beta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['beta']] <- beta
    else par.model[['beta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMgenfbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgenfbm <- new('RMmodelgenerator',
	.Data = RMgenfbm,
	type = 'negative definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgencauchy <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(beta) && !is.null(subst <- substitute(beta))) {
    u <- try(is.numeric(beta) || is.logical(beta) || is.language(beta)
	 || is.list(beta) || is(beta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['beta']] <- beta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['beta']] <- beta
    else par.model[['beta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMgencauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgencauchy <- new('RMmodelgenerator',
	.Data = RMgencauchy,
	type = 'undefined',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgengneiting <- function(kappa, mu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(kappa) && !is.null(subst <- substitute(kappa))) {
    u <- try(is.numeric(kappa) || is.logical(kappa) || is.language(kappa)
	 || is.list(kappa) || is(kappa, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['kappa']] <- kappa
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['kappa']] <- kappa
    else par.model[['kappa']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else par.model[['mu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMgengneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgengneiting <- new('RMmodelgenerator',
	.Data = RMgengneiting,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMgneiting <- function(orig, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(orig) && !is.null(subst <- substitute(orig))) {
    u <- try(is.numeric(orig) || is.logical(orig) || is.language(orig)
	 || is.list(orig) || is(orig, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['orig']] <- orig
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['orig']] <- orig
    else par.model[['orig']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMgneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgneiting <- new('RMmodelgenerator',
	.Data = RMgneiting,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMhyperbolic <- function(nu, lambda, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(lambda) && !is.null(subst <- substitute(lambda))) {
    u <- try(is.numeric(lambda) || is.logical(lambda) || is.language(lambda)
	 || is.list(lambda) || is(lambda, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['lambda']] <- lambda
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['lambda']] <- lambda
    else par.model[['lambda']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(delta) && !is.null(subst <- substitute(delta))) {
    u <- try(is.numeric(delta) || is.logical(delta) || is.language(delta)
	 || is.list(delta) || is(delta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['delta']] <- delta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['delta']] <- delta
    else par.model[['delta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMhyperbolic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMhyperbolic <- new('RMmodelgenerator',
	.Data = RMhyperbolic,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMiaco <- function(nu, lambda, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(lambda) && !is.null(subst <- substitute(lambda))) {
    u <- try(is.numeric(lambda) || is.logical(lambda) || is.language(lambda)
	 || is.list(lambda) || is(lambda, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['lambda']] <- lambda
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['lambda']] <- lambda
    else par.model[['lambda']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(delta) && !is.null(subst <- substitute(delta))) {
    u <- try(is.numeric(delta) || is.logical(delta) || is.language(delta)
	 || is.list(delta) || is(delta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['delta']] <- delta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['delta']] <- delta
    else par.model[['delta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMiaco', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMiaco <- new('RMmodelgenerator',
	.Data = RMiaco,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'space-isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMid <- function(phi, vdim, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(vdim) && !is.null(subst <- substitute(vdim))) {
    u <- try(is.numeric(vdim) || is.logical(vdim) || is.language(vdim)
	 || is.list(vdim) || is(vdim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['vdim']] <- vdim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['vdim']] <- vdim
    else par.model[['vdim']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMid', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMid <- new('RMmodelgenerator',
	.Data = RMid,
	type = 'undefined',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMkolmogorov <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMkolmogorov', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMkolmogorov <- new('RMmodelgenerator',
	.Data = RMkolmogorov,
	type = 'negative definite',
	domain = 'single variable',
	isotropy = 'vector-isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 3
	)



RMlgd <- function(alpha, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(beta) && !is.null(subst <- substitute(beta))) {
    u <- try(is.numeric(beta) || is.logical(beta) || is.language(beta)
	 || is.list(beta) || is(beta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['beta']] <- beta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['beta']] <- beta
    else par.model[['beta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMlgd', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMlgd <- new('RMmodelgenerator',
	.Data = RMlgd,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMmastein <- function(phi, nu, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(delta) && !is.null(subst <- substitute(delta))) {
    u <- try(is.numeric(delta) || is.logical(delta) || is.language(delta)
	 || is.list(delta) || is(delta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['delta']] <- delta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['delta']] <- delta
    else par.model[['delta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMmastein', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmastein <- new('RMmodelgenerator',
	.Data = RMmastein,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'space-isotropic',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMma <- function(phi, alpha, theta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(theta) && !is.null(subst <- substitute(theta))) {
    u <- try(is.numeric(theta) || is.logical(theta) || is.language(theta)
	 || is.list(theta) || is(theta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['theta']] <- theta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['theta']] <- theta
    else par.model[['theta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMma', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMma <- new('RMmodelgenerator',
	.Data = RMma,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMintexp <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMintexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMintexp <- new('RMmodelgenerator',
	.Data = RMintexp,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMmatrix <- function(phi, M, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(M) && !is.null(subst <- substitute(M))) {
    u <- try(is.numeric(M) || is.logical(M) || is.language(M)
	 || is.list(M) || is(M, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['M']] <- M
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['M']] <- M
    else par.model[['M']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMmatrix', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmatrix <- new('RMmodelgenerator',
	.Data = RMmatrix,
	type = 'positive definite',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMmatern <- function(nu, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(notinvnu) && !is.null(subst <- substitute(notinvnu))) {
    u <- try(is.numeric(notinvnu) || is.logical(notinvnu) || is.language(notinvnu)
	 || is.list(notinvnu) || is(notinvnu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['notinvnu']] <- notinvnu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['notinvnu']] <- notinvnu
    else par.model[['notinvnu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMmatern', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmatern <- new('RMmodelgenerator',
	.Data = RMmatern,
	type = 'undefined',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMmqam <- function(phi, C1, C2, C3, C4, C5, C6, C7, C8, C9, theta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg(theta) && !is.null(subst <- substitute(theta))) {
    u <- try(is.numeric(theta) || is.logical(theta) || is.language(theta)
	 || is.list(theta) || is(theta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['theta']] <- theta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['theta']] <- theta
    else par.model[['theta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMmqam', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmqam <- new('RMmodelgenerator',
	.Data = RMmqam,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMnatsc <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMnatsc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnatsc <- new('RMmodelgenerator',
	.Data = RMnatsc,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMnonstwm <- function(nu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMnonstwm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnonstwm <- new('RMmodelgenerator',
	.Data = RMnonstwm,
	type = 'positive definite',
	domain = 'kernel',
	isotropy = 'symmetric',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMnsst <- function(phi, psi, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(psi)) submodels[['psi']] <- psi
  
  if (hasArg(delta) && !is.null(subst <- substitute(delta))) {
    u <- try(is.numeric(delta) || is.logical(delta) || is.language(delta)
	 || is.list(delta) || is(delta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['delta']] <- delta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['delta']] <- delta
    else par.model[['delta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMnsst', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnsst <- new('RMmodelgenerator',
	.Data = RMnsst,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'space-isotropic',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMnugget <- function(tol, vdim, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(tol) && !is.null(subst <- substitute(tol))) {
    u <- try(is.numeric(tol) || is.logical(tol) || is.language(tol)
	 || is.list(tol) || is(tol, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tol']] <- tol
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tol']] <- tol
    else par.model[['tol']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(vdim) && !is.null(subst <- substitute(vdim))) {
    u <- try(is.numeric(vdim) || is.logical(vdim) || is.language(vdim)
	 || is.list(vdim) || is(vdim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['vdim']] <- vdim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['vdim']] <- vdim
    else par.model[['vdim']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMnugget', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnugget <- new('RMmodelgenerator',
	.Data = RMnugget,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -2
	)



RMflatpower <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMflatpower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMflatpower <- new('RMmodelgenerator',
	.Data = RMflatpower,
	type = 'negative definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'Bernstein',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMparswm <- function(nudiag, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nudiag) && !is.null(subst <- substitute(nudiag))) {
    u <- try(is.numeric(nudiag) || is.logical(nudiag) || is.language(nudiag)
	 || is.list(nudiag) || is(nudiag, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nudiag']] <- nudiag
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nudiag']] <- nudiag
    else par.model[['nudiag']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMparswm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMparswm <- new('RMmodelgenerator',
	.Data = RMparswm,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMpenta <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMpenta', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpenta <- new('RMmodelgenerator',
	.Data = RMpenta,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMaskey <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMaskey', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMaskey <- new('RMmodelgenerator',
	.Data = RMaskey,
	type = 'undefined',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMpower <- function(phi, alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMpower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpower <- new('RMmodelgenerator',
	.Data = RMpower,
	type = 'positive definite',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMqam <- function(phi, C1, C2, C3, C4, C5, C6, C7, C8, C9, theta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg(theta) && !is.null(subst <- substitute(theta))) {
    u <- try(is.numeric(theta) || is.logical(theta) || is.language(theta)
	 || is.list(theta) || is(theta, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['theta']] <- theta
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['theta']] <- theta
    else par.model[['theta']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMqam', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMqam <- new('RMmodelgenerator',
	.Data = RMqam,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMqexp <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMqexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMqexp <- new('RMmodelgenerator',
	.Data = RMqexp,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMschur <- function(phi, M, diag, rhored, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(M) && !is.null(subst <- substitute(M))) {
    u <- try(is.numeric(M) || is.logical(M) || is.language(M)
	 || is.list(M) || is(M, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['M']] <- M
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['M']] <- M
    else par.model[['M']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(diag) && !is.null(subst <- substitute(diag))) {
    u <- try(is.numeric(diag) || is.logical(diag) || is.language(diag)
	 || is.list(diag) || is(diag, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['diag']] <- diag
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['diag']] <- diag
    else par.model[['diag']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(rhored) && !is.null(subst <- substitute(rhored))) {
    u <- try(is.numeric(rhored) || is.logical(rhored) || is.language(rhored)
	 || is.list(rhored) || is(rhored, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['rhored']] <- rhored
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['rhored']] <- rhored
    else par.model[['rhored']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMschur', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMschur <- new('RMmodelgenerator',
	.Data = RMschur,
	type = 'positive definite',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMdelay <- function(phi, s, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else par.model[['s']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMdelay', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdelay <- new('RMmodelgenerator',
	.Data = RMdelay,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMspheric <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMspheric', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMspheric <- new('RMmodelgenerator',
	.Data = RMspheric,
	type = 'tail correlation function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'Gneiting-Schaback class',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMstable <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMstable', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstable <- new('RMmodelgenerator',
	.Data = RMstable,
	type = 'undefined',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMintrinsic <- function(phi, diameter, rawR, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(diameter) && !is.null(subst <- substitute(diameter))) {
    u <- try(is.numeric(diameter) || is.logical(diameter) || is.language(diameter)
	 || is.list(diameter) || is(diameter, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['diameter']] <- diameter
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['diameter']] <- diameter
    else par.model[['diameter']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(rawR) && !is.null(subst <- substitute(rawR))) {
    u <- try(is.numeric(rawR) || is.logical(rawR) || is.language(rawR)
	 || is.list(rawR) || is(rawR, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['rawR']] <- rawR
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['rawR']] <- rawR
    else par.model[['rawR']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMintrinsic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMintrinsic <- new('RMmodelgenerator',
	.Data = RMintrinsic,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RMstein <- function(nu, z, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(z) && !is.null(subst <- substitute(z))) {
    u <- try(is.numeric(z) || is.logical(z) || is.language(z)
	 || is.list(z) || is(z, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['z']] <- z
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['z']] <- z
    else par.model[['z']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMstein', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstein <- new('RMmodelgenerator',
	.Data = RMstein,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMstp <- function(xi, phi, S, z, M, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(xi)) submodels[['xi']] <- xi
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(S) && !is.null(subst <- substitute(S))) {
    u <- try(is.numeric(S) || is.logical(S) || is.language(S)
	 || is.list(S) || is(S, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['S']] <- S
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['S']] <- S
    else par.model[['S']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(z) && !is.null(subst <- substitute(z))) {
    u <- try(is.numeric(z) || is.logical(z) || is.language(z)
	 || is.list(z) || is(z, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['z']] <- z
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['z']] <- z
    else par.model[['z']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(M) && !is.null(subst <- substitute(M))) {
    u <- try(is.numeric(M) || is.logical(M) || is.language(M)
	 || is.list(M) || is(M, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['M']] <- M
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['M']] <- M
    else par.model[['M']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMstp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstp <- new('RMmodelgenerator',
	.Data = RMstp,
	type = 'positive definite',
	domain = 'kernel',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 1
	)



RMtbm <- function(phi, fulldim, reduceddim, layers, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(fulldim) && !is.null(subst <- substitute(fulldim))) {
    u <- try(is.numeric(fulldim) || is.logical(fulldim) || is.language(fulldim)
	 || is.list(fulldim) || is(fulldim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['fulldim']] <- fulldim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['fulldim']] <- fulldim
    else par.model[['fulldim']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(reduceddim) && !is.null(subst <- substitute(reduceddim))) {
    u <- try(is.numeric(reduceddim) || is.logical(reduceddim) || is.language(reduceddim)
	 || is.list(reduceddim) || is(reduceddim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['reduceddim']] <- reduceddim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['reduceddim']] <- reduceddim
    else par.model[['reduceddim']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(layers) && !is.null(subst <- substitute(layers))) {
    u <- try(is.numeric(layers) || is.logical(layers) || is.language(layers)
	 || is.list(layers) || is(layers, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['layers']] <- layers
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['layers']] <- layers
    else par.model[['layers']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMtbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtbm <- new('RMmodelgenerator',
	.Data = RMtbm,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = -3
	)



RMvector <- function(phi, a, Dspace, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(a) && !is.null(subst <- substitute(a))) {
    u <- try(is.numeric(a) || is.logical(a) || is.language(a)
	 || is.list(a) || is(a, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['a']] <- a
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['a']] <- a
    else par.model[['a']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Dspace) && !is.null(subst <- substitute(Dspace))) {
    u <- try(is.numeric(Dspace) || is.logical(Dspace) || is.language(Dspace)
	 || is.list(Dspace) || is(Dspace, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['Dspace']] <- Dspace
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['Dspace']] <- Dspace
    else par.model[['Dspace']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMvector', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMvector <- new('RMmodelgenerator',
	.Data = RMvector,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'symmetric',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMwave <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMwave', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMwave <- new('RMmodelgenerator',
	.Data = RMwave,
	type = 'positive definite',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMwhittle <- function(nu, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else par.model[['nu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(notinvnu) && !is.null(subst <- substitute(notinvnu))) {
    u <- try(is.numeric(notinvnu) || is.logical(notinvnu) || is.language(notinvnu)
	 || is.list(notinvnu) || is(notinvnu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['notinvnu']] <- notinvnu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['notinvnu']] <- notinvnu
    else par.model[['notinvnu']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMwhittle', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMwhittle <- new('RMmodelgenerator',
	.Data = RMwhittle,
	type = 'undefined',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMangle <- function(angle, lat.angle, ratio, diag) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(angle) && !is.null(subst <- substitute(angle))) {
    u <- try(is.numeric(angle) || is.logical(angle) || is.language(angle)
	 || is.list(angle) || is(angle, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['angle']] <- angle
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['angle']] <- angle
    else par.model[['angle']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(lat.angle) && !is.null(subst <- substitute(lat.angle))) {
    u <- try(is.numeric(lat.angle) || is.logical(lat.angle) || is.language(lat.angle)
	 || is.list(lat.angle) || is(lat.angle, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['lat.angle']] <- lat.angle
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['lat.angle']] <- lat.angle
    else par.model[['lat.angle']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(ratio) && !is.null(subst <- substitute(ratio))) {
    u <- try(is.numeric(ratio) || is.logical(ratio) || is.language(ratio)
	 || is.list(ratio) || is(ratio, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['ratio']] <- ratio
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['ratio']] <- ratio
    else par.model[['ratio']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(diag) && !is.null(subst <- substitute(diag))) {
    u <- try(is.numeric(diag) || is.logical(diag) || is.language(diag)
	 || is.list(diag) || is(diag, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['diag']] <- diag
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['diag']] <- diag
    else par.model[['diag']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMangle', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMangle <- new('RMmodelgenerator',
	.Data = RMangle,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMball <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(var) && !is.null(subst <- substitute(var))) {
    u <- try(is.numeric(var) || is.logical(var) || is.language(var)
	 || is.list(var) || is(var, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['var']] <- var
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['var']] <- var
    else par.general[['var']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['scale']] <- scale
    else par.general[['scale']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(Aniso) && !is.null(subst <- substitute(Aniso))) {
    u <- try(is.numeric(Aniso) || is.logical(Aniso) || is.language(Aniso)
	 || is.list(Aniso) || is(Aniso, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['Aniso']] <- Aniso
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['Aniso']] <- Aniso
    else par.general[['Aniso']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(proj) && !is.null(subst <- substitute(proj))) {
    u <- try(is.numeric(proj) || is.logical(proj) || is.language(proj)
	 || is.list(proj) || is(proj, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.general[['proj']] <- proj
    else if (substr(deparse(subst), 1, 1)=='R') par.general[['proj']] <- proj
    else par.general[['proj']] <- do.call('RRdistr', list(subst))
  }
  model <- new('RMmodel', call = cl, name = 'RMball', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMball <- new('RMmodelgenerator',
	.Data = RMball,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMeaxxa <- function(E, A) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(E) && !is.null(subst <- substitute(E))) {
    u <- try(is.numeric(E) || is.logical(E) || is.language(E)
	 || is.list(E) || is(E, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['E']] <- E
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['E']] <- E
    else par.model[['E']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(A) && !is.null(subst <- substitute(A))) {
    u <- try(is.numeric(A) || is.logical(A) || is.language(A)
	 || is.list(A) || is(A, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['A']] <- A
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['A']] <- A
    else par.model[['A']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMeaxxa', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMeaxxa <- new('RMmodelgenerator',
	.Data = RMeaxxa,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = -1
	)



RMetaxxa <- function(E, A, alpha) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(E) && !is.null(subst <- substitute(E))) {
    u <- try(is.numeric(E) || is.logical(E) || is.language(E)
	 || is.list(E) || is(E, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['E']] <- E
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['E']] <- E
    else par.model[['E']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(A) && !is.null(subst <- substitute(A))) {
    u <- try(is.numeric(A) || is.logical(A) || is.language(A)
	 || is.list(A) || is(A, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['A']] <- A
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['A']] <- A
    else par.model[['A']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else par.model[['alpha']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMetaxxa', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMetaxxa <- new('RMmodelgenerator',
	.Data = RMetaxxa,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 3
	)



RMtrafo <- function(isotropy) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(isotropy) && !is.null(subst <- substitute(isotropy))) {
    u <- try(is.numeric(isotropy) || is.logical(isotropy) || is.language(isotropy)
	 || is.list(isotropy) || is(isotropy, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['isotropy']] <- isotropy
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['isotropy']] <- isotropy
    else par.model[['isotropy']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMtrafo', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtrafo <- new('RMmodelgenerator',
	.Data = RMtrafo,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMpolygon <- function(lambda) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(lambda) && !is.null(subst <- substitute(lambda))) {
    u <- try(is.numeric(lambda) || is.logical(lambda) || is.language(lambda)
	 || is.list(lambda) || is(lambda, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['lambda']] <- lambda
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['lambda']] <- lambda
    else par.model[['lambda']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMpolygon', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpolygon <- new('RMmodelgenerator',
	.Data = RMpolygon,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMrational <- function(A, a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(A) && !is.null(subst <- substitute(A))) {
    u <- try(is.numeric(A) || is.logical(A) || is.language(A)
	 || is.list(A) || is(A, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['A']] <- A
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['A']] <- A
    else par.model[['A']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(a) && !is.null(subst <- substitute(a))) {
    u <- try(is.numeric(a) || is.logical(a) || is.language(a)
	 || is.list(a) || is(a, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['a']] <- a
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['a']] <- a
    else par.model[['a']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMrational', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrational <- new('RMmodelgenerator',
	.Data = RMrational,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMrotat <- function(speed, phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(speed) && !is.null(subst <- substitute(speed))) {
    u <- try(is.numeric(speed) || is.logical(speed) || is.language(speed)
	 || is.list(speed) || is(speed, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['speed']] <- speed
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['speed']] <- speed
    else par.model[['speed']] <- do.call('RRdistr', list(subst))
  }
  if (hasArg(phi) && !is.null(subst <- substitute(phi))) {
    u <- try(is.numeric(phi) || is.logical(phi) || is.language(phi)
	 || is.list(phi) || is(phi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['phi']] <- phi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['phi']] <- phi
    else par.model[['phi']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMrotat', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrotat <- new('RMmodelgenerator',
	.Data = RMrotat,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMrotation <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(phi) && !is.null(subst <- substitute(phi))) {
    u <- try(is.numeric(phi) || is.logical(phi) || is.language(phi)
	 || is.list(phi) || is(phi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['phi']] <- phi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['phi']] <- phi
    else par.model[['phi']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMrotation', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrotation <- new('RMmodelgenerator',
	.Data = RMrotation,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = -1
	)



RMsign <- function(phi, p) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(p) && !is.null(subst <- substitute(p))) {
    u <- try(is.numeric(p) || is.logical(p) || is.language(p)
	 || is.list(p) || is(p, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['p']] <- p
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['p']] <- p
    else par.model[['p']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMsign', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMsign <- new('RMmodelgenerator',
	.Data = RMsign,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMm2r <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  
  model <- new('RMmodel', call = cl, name = 'RMm2r', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMm2r <- new('RMmodelgenerator',
	.Data = RMm2r,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMm3b <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  
  model <- new('RMmodel', call = cl, name = 'RMm3b', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMm3b <- new('RMmodelgenerator',
	.Data = RMm3b,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'isotropic',
	operator = TRUE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMmps <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  
  model <- new('RMmodel', call = cl, name = 'RMmps', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmps <- new('RMmodelgenerator',
	.Data = RMmps,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = TRUE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMtruncsupport <- function(phi, radius) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(radius) && !is.null(subst <- substitute(radius))) {
    u <- try(is.numeric(radius) || is.logical(radius) || is.language(radius)
	 || is.list(radius) || is(radius, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['radius']] <- radius
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['radius']] <- radius
    else par.model[['radius']] <- do.call('RRdistr', list(subst))
  }
  
  model <- new('RMmodel', call = cl, name = 'RMtruncsupport', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtruncsupport <- new('RMmodelgenerator',
	.Data = RMtruncsupport,
	type = 'shape function',
	domain = 'single variable',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RRdeterm <- function(mean) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(mean) && !is.null(subst <- substitute(mean))) {
    u <- try(is.numeric(mean) || is.logical(mean) || is.language(mean)
	 || is.list(mean) || is(mean, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mean']] <- mean
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mean']] <- mean
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RRdeterm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRdeterm <- new('RMmodelgenerator',
	.Data = RRdeterm,
	type = 'distribution family',
	domain = 'framework dependent',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RRgauss <- function(mu, sd, log) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(sd) && !is.null(subst <- substitute(sd))) {
    u <- try(is.numeric(sd) || is.logical(sd) || is.language(sd)
	 || is.list(sd) || is(sd, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['sd']] <- sd
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['sd']] <- sd
    else  stop('random parameter not allowed')
  }
  if (hasArg(log) && !is.null(subst <- substitute(log))) {
    u <- try(is.numeric(log) || is.logical(log) || is.language(log)
	 || is.list(log) || is(log, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['log']] <- log
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['log']] <- log
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RRgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRgauss <- new('RMmodelgenerator',
	.Data = RRgauss,
	type = 'distribution family',
	domain = 'framework dependent',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RRloc <- function(phi, mu, scale, pow) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(scale) && !is.null(subst <- substitute(scale))) {
    u <- try(is.numeric(scale) || is.logical(scale) || is.language(scale)
	 || is.list(scale) || is(scale, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['scale']] <- scale
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['scale']] <- scale
    else  stop('random parameter not allowed')
  }
  if (hasArg(pow) && !is.null(subst <- substitute(pow))) {
    u <- try(is.numeric(pow) || is.logical(pow) || is.language(pow)
	 || is.list(pow) || is(pow, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['pow']] <- pow
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['pow']] <- pow
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RRloc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRloc <- new('RMmodelgenerator',
	.Data = RRloc,
	type = 'distribution family',
	domain = 'framework dependent',
	isotropy = 'cartesian system',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RRrectangular <- function(phi, safety, minsteplen, maxsteps, parts, maxit, innermin, outermax, mcmc_n, normed, approx, onesided) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(safety) && !is.null(subst <- substitute(safety))) {
    u <- try(is.numeric(safety) || is.logical(safety) || is.language(safety)
	 || is.list(safety) || is(safety, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['safety']] <- safety
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['safety']] <- safety
    else  stop('random parameter not allowed')
  }
  if (hasArg(minsteplen) && !is.null(subst <- substitute(minsteplen))) {
    u <- try(is.numeric(minsteplen) || is.logical(minsteplen) || is.language(minsteplen)
	 || is.list(minsteplen) || is(minsteplen, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['minsteplen']] <- minsteplen
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['minsteplen']] <- minsteplen
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxsteps) && !is.null(subst <- substitute(maxsteps))) {
    u <- try(is.numeric(maxsteps) || is.logical(maxsteps) || is.language(maxsteps)
	 || is.list(maxsteps) || is(maxsteps, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxsteps']] <- maxsteps
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxsteps']] <- maxsteps
    else  stop('random parameter not allowed')
  }
  if (hasArg(parts) && !is.null(subst <- substitute(parts))) {
    u <- try(is.numeric(parts) || is.logical(parts) || is.language(parts)
	 || is.list(parts) || is(parts, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['parts']] <- parts
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['parts']] <- parts
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxit) && !is.null(subst <- substitute(maxit))) {
    u <- try(is.numeric(maxit) || is.logical(maxit) || is.language(maxit)
	 || is.list(maxit) || is(maxit, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxit']] <- maxit
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxit']] <- maxit
    else  stop('random parameter not allowed')
  }
  if (hasArg(innermin) && !is.null(subst <- substitute(innermin))) {
    u <- try(is.numeric(innermin) || is.logical(innermin) || is.language(innermin)
	 || is.list(innermin) || is(innermin, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['innermin']] <- innermin
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['innermin']] <- innermin
    else  stop('random parameter not allowed')
  }
  if (hasArg(outermax) && !is.null(subst <- substitute(outermax))) {
    u <- try(is.numeric(outermax) || is.logical(outermax) || is.language(outermax)
	 || is.list(outermax) || is(outermax, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['outermax']] <- outermax
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['outermax']] <- outermax
    else  stop('random parameter not allowed')
  }
  if (hasArg(mcmc_n) && !is.null(subst <- substitute(mcmc_n))) {
    u <- try(is.numeric(mcmc_n) || is.logical(mcmc_n) || is.language(mcmc_n)
	 || is.list(mcmc_n) || is(mcmc_n, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mcmc_n']] <- mcmc_n
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mcmc_n']] <- mcmc_n
    else  stop('random parameter not allowed')
  }
  if (hasArg(normed) && !is.null(subst <- substitute(normed))) {
    u <- try(is.numeric(normed) || is.logical(normed) || is.language(normed)
	 || is.list(normed) || is(normed, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['normed']] <- normed
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['normed']] <- normed
    else  stop('random parameter not allowed')
  }
  if (hasArg(approx) && !is.null(subst <- substitute(approx))) {
    u <- try(is.numeric(approx) || is.logical(approx) || is.language(approx)
	 || is.list(approx) || is(approx, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['approx']] <- approx
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['approx']] <- approx
    else  stop('random parameter not allowed')
  }
  if (hasArg(onesided) && !is.null(subst <- substitute(onesided))) {
    u <- try(is.numeric(onesided) || is.logical(onesided) || is.language(onesided)
	 || is.list(onesided) || is(onesided, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['onesided']] <- onesided
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['onesided']] <- onesided
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RRrectangular', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRrectangular <- new('RMmodelgenerator',
	.Data = RRrectangular,
	type = 'distribution family',
	domain = 'framework dependent',
	isotropy = 'cartesian system',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RRspheric <- function(spacedim, balldim, R) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(spacedim) && !is.null(subst <- substitute(spacedim))) {
    u <- try(is.numeric(spacedim) || is.logical(spacedim) || is.language(spacedim)
	 || is.list(spacedim) || is(spacedim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['spacedim']] <- spacedim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['spacedim']] <- spacedim
    else  stop('random parameter not allowed')
  }
  if (hasArg(balldim) && !is.null(subst <- substitute(balldim))) {
    u <- try(is.numeric(balldim) || is.logical(balldim) || is.language(balldim)
	 || is.list(balldim) || is(balldim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['balldim']] <- balldim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['balldim']] <- balldim
    else  stop('random parameter not allowed')
  }
  if (hasArg(R) && !is.null(subst <- substitute(R))) {
    u <- try(is.numeric(R) || is.logical(R) || is.language(R)
	 || is.list(R) || is(R, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['R']] <- R
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['R']] <- R
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RRspheric', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRspheric <- new('RMmodelgenerator',
	.Data = RRspheric,
	type = 'distribution family',
	domain = 'single variable',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



RRunif <- function(min, max, normed) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg(min) && !is.null(subst <- substitute(min))) {
    u <- try(is.numeric(min) || is.logical(min) || is.language(min)
	 || is.list(min) || is(min, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['min']] <- min
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['min']] <- min
    else  stop('random parameter not allowed')
  }
  if (hasArg(max) && !is.null(subst <- substitute(max))) {
    u <- try(is.numeric(max) || is.logical(max) || is.language(max)
	 || is.list(max) || is(max, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['max']] <- max
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['max']] <- max
    else  stop('random parameter not allowed')
  }
  if (hasArg(normed) && !is.null(subst <- substitute(normed))) {
    u <- try(is.numeric(normed) || is.logical(normed) || is.language(normed)
	 || is.list(normed) || is(normed, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['normed']] <- normed
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['normed']] <- normed
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RRunif', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRunif <- new('RMmodelgenerator',
	.Data = RRunif,
	type = 'distribution family',
	domain = 'framework dependent',
	isotropy = 'cartesian system',
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMmppplus <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, p) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(C0)) submodels[['C0']] <- C0
  if (hasArg(C1)) submodels[['C1']] <- C1
  if (hasArg(C2)) submodels[['C2']] <- C2
  if (hasArg(C3)) submodels[['C3']] <- C3
  if (hasArg(C4)) submodels[['C4']] <- C4
  if (hasArg(C5)) submodels[['C5']] <- C5
  if (hasArg(C6)) submodels[['C6']] <- C6
  if (hasArg(C7)) submodels[['C7']] <- C7
  if (hasArg(C8)) submodels[['C8']] <- C8
  if (hasArg(C9)) submodels[['C9']] <- C9
  
  if (hasArg(p) && !is.null(subst <- substitute(p))) {
    u <- try(is.numeric(p) || is.logical(p) || is.language(p)
	 || is.list(p) || is(p, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['p']] <- p
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['p']] <- p
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RMmppplus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmppplus <- new('RMmodelgenerator',
	.Data = RMmppplus,
	type = 'shifted shape function',
	domain = 'framework dependent',
	isotropy = 'parameter dependent',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RPaverage <- function(phi, shape, intensity) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(shape)) submodels[['shape']] <- shape
  
  if (hasArg(intensity) && !is.null(subst <- substitute(intensity))) {
    u <- try(is.numeric(intensity) || is.logical(intensity) || is.language(intensity)
	 || is.list(intensity) || is(intensity, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['intensity']] <- intensity
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['intensity']] <- intensity
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPaverage', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPaverage <- new('RMmodelgenerator',
	.Data = RPaverage,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPcirculant <- function(phi, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(force) && !is.null(subst <- substitute(force))) {
    u <- try(is.numeric(force) || is.logical(force) || is.language(force)
	 || is.list(force) || is(force, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['force']] <- force
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['force']] <- force
    else  stop('random parameter not allowed')
  }
  if (hasArg(mmin) && !is.null(subst <- substitute(mmin))) {
    u <- try(is.numeric(mmin) || is.logical(mmin) || is.language(mmin)
	 || is.list(mmin) || is(mmin, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mmin']] <- mmin
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mmin']] <- mmin
    else  stop('random parameter not allowed')
  }
  if (hasArg(strategy) && !is.null(subst <- substitute(strategy))) {
    u <- try(is.numeric(strategy) || is.logical(strategy) || is.language(strategy)
	 || is.list(strategy) || is(strategy, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['strategy']] <- strategy
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['strategy']] <- strategy
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxGB) && !is.null(subst <- substitute(maxGB))) {
    u <- try(is.numeric(maxGB) || is.logical(maxGB) || is.language(maxGB)
	 || is.list(maxGB) || is(maxGB, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxGB']] <- maxGB
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxGB']] <- maxGB
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxmem) && !is.null(subst <- substitute(maxmem))) {
    u <- try(is.numeric(maxmem) || is.logical(maxmem) || is.language(maxmem)
	 || is.list(maxmem) || is(maxmem, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxmem']] <- maxmem
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxmem']] <- maxmem
    else  stop('random parameter not allowed')
  }
  if (hasArg(tolIm) && !is.null(subst <- substitute(tolIm))) {
    u <- try(is.numeric(tolIm) || is.logical(tolIm) || is.language(tolIm)
	 || is.list(tolIm) || is(tolIm, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tolIm']] <- tolIm
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tolIm']] <- tolIm
    else  stop('random parameter not allowed')
  }
  if (hasArg(tolRe) && !is.null(subst <- substitute(tolRe))) {
    u <- try(is.numeric(tolRe) || is.logical(tolRe) || is.language(tolRe)
	 || is.list(tolRe) || is(tolRe, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tolRe']] <- tolRe
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tolRe']] <- tolRe
    else  stop('random parameter not allowed')
  }
  if (hasArg(trials) && !is.null(subst <- substitute(trials))) {
    u <- try(is.numeric(trials) || is.logical(trials) || is.language(trials)
	 || is.list(trials) || is(trials, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['trials']] <- trials
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['trials']] <- trials
    else  stop('random parameter not allowed')
  }
  if (hasArg(useprimes) && !is.null(subst <- substitute(useprimes))) {
    u <- try(is.numeric(useprimes) || is.logical(useprimes) || is.language(useprimes)
	 || is.list(useprimes) || is(useprimes, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['useprimes']] <- useprimes
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['useprimes']] <- useprimes
    else  stop('random parameter not allowed')
  }
  if (hasArg(dependent) && !is.null(subst <- substitute(dependent))) {
    u <- try(is.numeric(dependent) || is.logical(dependent) || is.language(dependent)
	 || is.list(dependent) || is(dependent, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['dependent']] <- dependent
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['dependent']] <- dependent
    else  stop('random parameter not allowed')
  }
  if (hasArg(approx_step) && !is.null(subst <- substitute(approx_step))) {
    u <- try(is.numeric(approx_step) || is.logical(approx_step) || is.language(approx_step)
	 || is.list(approx_step) || is(approx_step, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['approx_step']] <- approx_step
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['approx_step']] <- approx_step
    else  stop('random parameter not allowed')
  }
  if (hasArg(approx_maxgrid) && !is.null(subst <- substitute(approx_maxgrid))) {
    u <- try(is.numeric(approx_maxgrid) || is.logical(approx_maxgrid) || is.language(approx_maxgrid)
	 || is.list(approx_maxgrid) || is(approx_maxgrid, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['approx_maxgrid']] <- approx_maxgrid
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['approx_maxgrid']] <- approx_maxgrid
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPcirculant', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcirculant <- new('RMmodelgenerator',
	.Data = RPcirculant,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = -3
	)



RPcutoff <- function(phi, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid, diameter, a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(force) && !is.null(subst <- substitute(force))) {
    u <- try(is.numeric(force) || is.logical(force) || is.language(force)
	 || is.list(force) || is(force, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['force']] <- force
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['force']] <- force
    else  stop('random parameter not allowed')
  }
  if (hasArg(mmin) && !is.null(subst <- substitute(mmin))) {
    u <- try(is.numeric(mmin) || is.logical(mmin) || is.language(mmin)
	 || is.list(mmin) || is(mmin, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mmin']] <- mmin
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mmin']] <- mmin
    else  stop('random parameter not allowed')
  }
  if (hasArg(strategy) && !is.null(subst <- substitute(strategy))) {
    u <- try(is.numeric(strategy) || is.logical(strategy) || is.language(strategy)
	 || is.list(strategy) || is(strategy, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['strategy']] <- strategy
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['strategy']] <- strategy
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxGB) && !is.null(subst <- substitute(maxGB))) {
    u <- try(is.numeric(maxGB) || is.logical(maxGB) || is.language(maxGB)
	 || is.list(maxGB) || is(maxGB, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxGB']] <- maxGB
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxGB']] <- maxGB
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxmem) && !is.null(subst <- substitute(maxmem))) {
    u <- try(is.numeric(maxmem) || is.logical(maxmem) || is.language(maxmem)
	 || is.list(maxmem) || is(maxmem, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxmem']] <- maxmem
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxmem']] <- maxmem
    else  stop('random parameter not allowed')
  }
  if (hasArg(tolIm) && !is.null(subst <- substitute(tolIm))) {
    u <- try(is.numeric(tolIm) || is.logical(tolIm) || is.language(tolIm)
	 || is.list(tolIm) || is(tolIm, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tolIm']] <- tolIm
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tolIm']] <- tolIm
    else  stop('random parameter not allowed')
  }
  if (hasArg(tolRe) && !is.null(subst <- substitute(tolRe))) {
    u <- try(is.numeric(tolRe) || is.logical(tolRe) || is.language(tolRe)
	 || is.list(tolRe) || is(tolRe, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tolRe']] <- tolRe
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tolRe']] <- tolRe
    else  stop('random parameter not allowed')
  }
  if (hasArg(trials) && !is.null(subst <- substitute(trials))) {
    u <- try(is.numeric(trials) || is.logical(trials) || is.language(trials)
	 || is.list(trials) || is(trials, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['trials']] <- trials
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['trials']] <- trials
    else  stop('random parameter not allowed')
  }
  if (hasArg(useprimes) && !is.null(subst <- substitute(useprimes))) {
    u <- try(is.numeric(useprimes) || is.logical(useprimes) || is.language(useprimes)
	 || is.list(useprimes) || is(useprimes, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['useprimes']] <- useprimes
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['useprimes']] <- useprimes
    else  stop('random parameter not allowed')
  }
  if (hasArg(dependent) && !is.null(subst <- substitute(dependent))) {
    u <- try(is.numeric(dependent) || is.logical(dependent) || is.language(dependent)
	 || is.list(dependent) || is(dependent, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['dependent']] <- dependent
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['dependent']] <- dependent
    else  stop('random parameter not allowed')
  }
  if (hasArg(approx_step) && !is.null(subst <- substitute(approx_step))) {
    u <- try(is.numeric(approx_step) || is.logical(approx_step) || is.language(approx_step)
	 || is.list(approx_step) || is(approx_step, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['approx_step']] <- approx_step
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['approx_step']] <- approx_step
    else  stop('random parameter not allowed')
  }
  if (hasArg(approx_maxgrid) && !is.null(subst <- substitute(approx_maxgrid))) {
    u <- try(is.numeric(approx_maxgrid) || is.logical(approx_maxgrid) || is.language(approx_maxgrid)
	 || is.list(approx_maxgrid) || is(approx_maxgrid, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['approx_maxgrid']] <- approx_maxgrid
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['approx_maxgrid']] <- approx_maxgrid
    else  stop('random parameter not allowed')
  }
  if (hasArg(diameter) && !is.null(subst <- substitute(diameter))) {
    u <- try(is.numeric(diameter) || is.logical(diameter) || is.language(diameter)
	 || is.list(diameter) || is(diameter, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['diameter']] <- diameter
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['diameter']] <- diameter
    else  stop('random parameter not allowed')
  }
  if (hasArg(a) && !is.null(subst <- substitute(a))) {
    u <- try(is.numeric(a) || is.logical(a) || is.language(a)
	 || is.list(a) || is(a, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['a']] <- a
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['a']] <- a
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPcutoff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcutoff <- new('RMmodelgenerator',
	.Data = RPcutoff,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RPintrinsic <- function(phi, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid, diameter, rawR) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(force) && !is.null(subst <- substitute(force))) {
    u <- try(is.numeric(force) || is.logical(force) || is.language(force)
	 || is.list(force) || is(force, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['force']] <- force
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['force']] <- force
    else  stop('random parameter not allowed')
  }
  if (hasArg(mmin) && !is.null(subst <- substitute(mmin))) {
    u <- try(is.numeric(mmin) || is.logical(mmin) || is.language(mmin)
	 || is.list(mmin) || is(mmin, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mmin']] <- mmin
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mmin']] <- mmin
    else  stop('random parameter not allowed')
  }
  if (hasArg(strategy) && !is.null(subst <- substitute(strategy))) {
    u <- try(is.numeric(strategy) || is.logical(strategy) || is.language(strategy)
	 || is.list(strategy) || is(strategy, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['strategy']] <- strategy
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['strategy']] <- strategy
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxGB) && !is.null(subst <- substitute(maxGB))) {
    u <- try(is.numeric(maxGB) || is.logical(maxGB) || is.language(maxGB)
	 || is.list(maxGB) || is(maxGB, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxGB']] <- maxGB
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxGB']] <- maxGB
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxmem) && !is.null(subst <- substitute(maxmem))) {
    u <- try(is.numeric(maxmem) || is.logical(maxmem) || is.language(maxmem)
	 || is.list(maxmem) || is(maxmem, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxmem']] <- maxmem
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxmem']] <- maxmem
    else  stop('random parameter not allowed')
  }
  if (hasArg(tolIm) && !is.null(subst <- substitute(tolIm))) {
    u <- try(is.numeric(tolIm) || is.logical(tolIm) || is.language(tolIm)
	 || is.list(tolIm) || is(tolIm, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tolIm']] <- tolIm
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tolIm']] <- tolIm
    else  stop('random parameter not allowed')
  }
  if (hasArg(tolRe) && !is.null(subst <- substitute(tolRe))) {
    u <- try(is.numeric(tolRe) || is.logical(tolRe) || is.language(tolRe)
	 || is.list(tolRe) || is(tolRe, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tolRe']] <- tolRe
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tolRe']] <- tolRe
    else  stop('random parameter not allowed')
  }
  if (hasArg(trials) && !is.null(subst <- substitute(trials))) {
    u <- try(is.numeric(trials) || is.logical(trials) || is.language(trials)
	 || is.list(trials) || is(trials, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['trials']] <- trials
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['trials']] <- trials
    else  stop('random parameter not allowed')
  }
  if (hasArg(useprimes) && !is.null(subst <- substitute(useprimes))) {
    u <- try(is.numeric(useprimes) || is.logical(useprimes) || is.language(useprimes)
	 || is.list(useprimes) || is(useprimes, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['useprimes']] <- useprimes
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['useprimes']] <- useprimes
    else  stop('random parameter not allowed')
  }
  if (hasArg(dependent) && !is.null(subst <- substitute(dependent))) {
    u <- try(is.numeric(dependent) || is.logical(dependent) || is.language(dependent)
	 || is.list(dependent) || is(dependent, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['dependent']] <- dependent
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['dependent']] <- dependent
    else  stop('random parameter not allowed')
  }
  if (hasArg(approx_step) && !is.null(subst <- substitute(approx_step))) {
    u <- try(is.numeric(approx_step) || is.logical(approx_step) || is.language(approx_step)
	 || is.list(approx_step) || is(approx_step, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['approx_step']] <- approx_step
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['approx_step']] <- approx_step
    else  stop('random parameter not allowed')
  }
  if (hasArg(approx_maxgrid) && !is.null(subst <- substitute(approx_maxgrid))) {
    u <- try(is.numeric(approx_maxgrid) || is.logical(approx_maxgrid) || is.language(approx_maxgrid)
	 || is.list(approx_maxgrid) || is(approx_maxgrid, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['approx_maxgrid']] <- approx_maxgrid
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['approx_maxgrid']] <- approx_maxgrid
    else  stop('random parameter not allowed')
  }
  if (hasArg(diameter) && !is.null(subst <- substitute(diameter))) {
    u <- try(is.numeric(diameter) || is.logical(diameter) || is.language(diameter)
	 || is.list(diameter) || is(diameter, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['diameter']] <- diameter
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['diameter']] <- diameter
    else  stop('random parameter not allowed')
  }
  if (hasArg(rawR) && !is.null(subst <- substitute(rawR))) {
    u <- try(is.numeric(rawR) || is.logical(rawR) || is.language(rawR)
	 || is.list(rawR) || is(rawR, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['rawR']] <- rawR
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['rawR']] <- rawR
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPintrinsic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPintrinsic <- new('RMmodelgenerator',
	.Data = RPintrinsic,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RPdirect <- function(phi, root_method, svdtolerance, max_variab) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(root_method) && !is.null(subst <- substitute(root_method))) {
    u <- try(is.numeric(root_method) || is.logical(root_method) || is.language(root_method)
	 || is.list(root_method) || is(root_method, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['root_method']] <- root_method
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['root_method']] <- root_method
    else  stop('random parameter not allowed')
  }
  if (hasArg(svdtolerance) && !is.null(subst <- substitute(svdtolerance))) {
    u <- try(is.numeric(svdtolerance) || is.logical(svdtolerance) || is.language(svdtolerance)
	 || is.list(svdtolerance) || is(svdtolerance, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['svdtolerance']] <- svdtolerance
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['svdtolerance']] <- svdtolerance
    else  stop('random parameter not allowed')
  }
  if (hasArg(max_variab) && !is.null(subst <- substitute(max_variab))) {
    u <- try(is.numeric(max_variab) || is.logical(max_variab) || is.language(max_variab)
	 || is.list(max_variab) || is(max_variab, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['max_variab']] <- max_variab
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['max_variab']] <- max_variab
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPdirect', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPdirect <- new('RMmodelgenerator',
	.Data = RPdirect,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPhyperplane <- function(phi, superpos, maxlines, mar_distr, mar_param) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(superpos) && !is.null(subst <- substitute(superpos))) {
    u <- try(is.numeric(superpos) || is.logical(superpos) || is.language(superpos)
	 || is.list(superpos) || is(superpos, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['superpos']] <- superpos
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['superpos']] <- superpos
    else  stop('random parameter not allowed')
  }
  if (hasArg(maxlines) && !is.null(subst <- substitute(maxlines))) {
    u <- try(is.numeric(maxlines) || is.logical(maxlines) || is.language(maxlines)
	 || is.list(maxlines) || is(maxlines, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['maxlines']] <- maxlines
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['maxlines']] <- maxlines
    else  stop('random parameter not allowed')
  }
  if (hasArg(mar_distr) && !is.null(subst <- substitute(mar_distr))) {
    u <- try(is.numeric(mar_distr) || is.logical(mar_distr) || is.language(mar_distr)
	 || is.list(mar_distr) || is(mar_distr, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mar_distr']] <- mar_distr
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mar_distr']] <- mar_distr
    else  stop('random parameter not allowed')
  }
  if (hasArg(mar_param) && !is.null(subst <- substitute(mar_param))) {
    u <- try(is.numeric(mar_param) || is.logical(mar_param) || is.language(mar_param)
	 || is.list(mar_param) || is(mar_param, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mar_param']] <- mar_param
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mar_param']] <- mar_param
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPhyperplane', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPhyperplane <- new('RMmodelgenerator',
	.Data = RPhyperplane,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RPnugget <- function(phi, tol, vdim) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(tol) && !is.null(subst <- substitute(tol))) {
    u <- try(is.numeric(tol) || is.logical(tol) || is.language(tol)
	 || is.list(tol) || is(tol, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['tol']] <- tol
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['tol']] <- tol
    else  stop('random parameter not allowed')
  }
  if (hasArg(vdim) && !is.null(subst <- substitute(vdim))) {
    u <- try(is.numeric(vdim) || is.logical(vdim) || is.language(vdim)
	 || is.list(vdim) || is(vdim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['vdim']] <- vdim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['vdim']] <- vdim
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPnugget', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPnugget <- new('RMmodelgenerator',
	.Data = RPnugget,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -2
	)



RPcoins <- function(phi, shape, intensity) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(shape)) submodels[['shape']] <- shape
  
  if (hasArg(intensity) && !is.null(subst <- substitute(intensity))) {
    u <- try(is.numeric(intensity) || is.logical(intensity) || is.language(intensity)
	 || is.list(intensity) || is(intensity, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['intensity']] <- intensity
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['intensity']] <- intensity
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPcoins', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcoins <- new('RMmodelgenerator',
	.Data = RPcoins,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPsequential <- function(phi, max_variables, back_steps, initial) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(max_variables) && !is.null(subst <- substitute(max_variables))) {
    u <- try(is.numeric(max_variables) || is.logical(max_variables) || is.language(max_variables)
	 || is.list(max_variables) || is(max_variables, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['max_variables']] <- max_variables
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['max_variables']] <- max_variables
    else  stop('random parameter not allowed')
  }
  if (hasArg(back_steps) && !is.null(subst <- substitute(back_steps))) {
    u <- try(is.numeric(back_steps) || is.logical(back_steps) || is.language(back_steps)
	 || is.list(back_steps) || is(back_steps, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['back_steps']] <- back_steps
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['back_steps']] <- back_steps
    else  stop('random parameter not allowed')
  }
  if (hasArg(initial) && !is.null(subst <- substitute(initial))) {
    u <- try(is.numeric(initial) || is.logical(initial) || is.language(initial)
	 || is.list(initial) || is(initial, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['initial']] <- initial
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['initial']] <- initial
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPsequential', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPsequential <- new('RMmodelgenerator',
	.Data = RPsequential,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RPspectral <- function(phi, sp_lines, sp_grid, prop_factor, sigma) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(sp_lines) && !is.null(subst <- substitute(sp_lines))) {
    u <- try(is.numeric(sp_lines) || is.logical(sp_lines) || is.language(sp_lines)
	 || is.list(sp_lines) || is(sp_lines, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['sp_lines']] <- sp_lines
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['sp_lines']] <- sp_lines
    else  stop('random parameter not allowed')
  }
  if (hasArg(sp_grid) && !is.null(subst <- substitute(sp_grid))) {
    u <- try(is.numeric(sp_grid) || is.logical(sp_grid) || is.language(sp_grid)
	 || is.list(sp_grid) || is(sp_grid, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['sp_grid']] <- sp_grid
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['sp_grid']] <- sp_grid
    else  stop('random parameter not allowed')
  }
  if (hasArg(prop_factor) && !is.null(subst <- substitute(prop_factor))) {
    u <- try(is.numeric(prop_factor) || is.logical(prop_factor) || is.language(prop_factor)
	 || is.list(prop_factor) || is(prop_factor, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['prop_factor']] <- prop_factor
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['prop_factor']] <- prop_factor
    else  stop('random parameter not allowed')
  }
  if (hasArg(sigma) && !is.null(subst <- substitute(sigma))) {
    u <- try(is.numeric(sigma) || is.logical(sigma) || is.language(sigma)
	 || is.list(sigma) || is(sigma, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['sigma']] <- sigma
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['sigma']] <- sigma
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPspectral', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPspectral <- new('RMmodelgenerator',
	.Data = RPspectral,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPspecific <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  
  model <- new('RMmodel', call = cl, name = 'RPspecific', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPspecific <- new('RMmodelgenerator',
	.Data = RPspecific,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPtbm <- function(phi, fulldim, reduceddim, layers, lines, linessimufactor, linesimustep, center, points) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(fulldim) && !is.null(subst <- substitute(fulldim))) {
    u <- try(is.numeric(fulldim) || is.logical(fulldim) || is.language(fulldim)
	 || is.list(fulldim) || is(fulldim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['fulldim']] <- fulldim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['fulldim']] <- fulldim
    else  stop('random parameter not allowed')
  }
  if (hasArg(reduceddim) && !is.null(subst <- substitute(reduceddim))) {
    u <- try(is.numeric(reduceddim) || is.logical(reduceddim) || is.language(reduceddim)
	 || is.list(reduceddim) || is(reduceddim, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['reduceddim']] <- reduceddim
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['reduceddim']] <- reduceddim
    else  stop('random parameter not allowed')
  }
  if (hasArg(layers) && !is.null(subst <- substitute(layers))) {
    u <- try(is.numeric(layers) || is.logical(layers) || is.language(layers)
	 || is.list(layers) || is(layers, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['layers']] <- layers
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['layers']] <- layers
    else  stop('random parameter not allowed')
  }
  if (hasArg(lines) && !is.null(subst <- substitute(lines))) {
    u <- try(is.numeric(lines) || is.logical(lines) || is.language(lines)
	 || is.list(lines) || is(lines, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['lines']] <- lines
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['lines']] <- lines
    else  stop('random parameter not allowed')
  }
  if (hasArg(linessimufactor) && !is.null(subst <- substitute(linessimufactor))) {
    u <- try(is.numeric(linessimufactor) || is.logical(linessimufactor) || is.language(linessimufactor)
	 || is.list(linessimufactor) || is(linessimufactor, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['linessimufactor']] <- linessimufactor
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['linessimufactor']] <- linessimufactor
    else  stop('random parameter not allowed')
  }
  if (hasArg(linesimustep) && !is.null(subst <- substitute(linesimustep))) {
    u <- try(is.numeric(linesimustep) || is.logical(linesimustep) || is.language(linesimustep)
	 || is.list(linesimustep) || is(linesimustep, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['linesimustep']] <- linesimustep
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['linesimustep']] <- linesimustep
    else  stop('random parameter not allowed')
  }
  if (hasArg(center) && !is.null(subst <- substitute(center))) {
    u <- try(is.numeric(center) || is.logical(center) || is.language(center)
	 || is.list(center) || is(center, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['center']] <- center
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['center']] <- center
    else  stop('random parameter not allowed')
  }
  if (hasArg(points) && !is.null(subst <- substitute(points))) {
    u <- try(is.numeric(points) || is.logical(points) || is.language(points)
	 || is.list(points) || is(points, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['points']] <- points
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['points']] <- points
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPtbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPtbm <- new('RMmodelgenerator',
	.Data = RPtbm,
	type = 'method for Gauss processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RPbrorig <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg(xi) && !is.null(subst <- substitute(xi))) {
    u <- try(is.numeric(xi) || is.logical(xi) || is.language(xi)
	 || is.list(xi) || is(xi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['xi']] <- xi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['xi']] <- xi
    else  stop('random parameter not allowed')
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPbrorig', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrorig <- new('RMmodelgenerator',
	.Data = RPbrorig,
	type = 'method for Brown-Resnick processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPbrmixed <- function(phi, tcf, xi, mu, s, meshsize, vertnumber, optim_mixed, optim_mixed_tol, optim_mixed_maxpo, lambda, areamat, variobound) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg(xi) && !is.null(subst <- substitute(xi))) {
    u <- try(is.numeric(xi) || is.logical(xi) || is.language(xi)
	 || is.list(xi) || is(xi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['xi']] <- xi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['xi']] <- xi
    else  stop('random parameter not allowed')
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else  stop('random parameter not allowed')
  }
  if (hasArg(meshsize) && !is.null(subst <- substitute(meshsize))) {
    u <- try(is.numeric(meshsize) || is.logical(meshsize) || is.language(meshsize)
	 || is.list(meshsize) || is(meshsize, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['meshsize']] <- meshsize
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['meshsize']] <- meshsize
    else  stop('random parameter not allowed')
  }
  if (hasArg(vertnumber) && !is.null(subst <- substitute(vertnumber))) {
    u <- try(is.numeric(vertnumber) || is.logical(vertnumber) || is.language(vertnumber)
	 || is.list(vertnumber) || is(vertnumber, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['vertnumber']] <- vertnumber
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['vertnumber']] <- vertnumber
    else  stop('random parameter not allowed')
  }
  if (hasArg(optim_mixed) && !is.null(subst <- substitute(optim_mixed))) {
    u <- try(is.numeric(optim_mixed) || is.logical(optim_mixed) || is.language(optim_mixed)
	 || is.list(optim_mixed) || is(optim_mixed, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['optim_mixed']] <- optim_mixed
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['optim_mixed']] <- optim_mixed
    else  stop('random parameter not allowed')
  }
  if (hasArg(optim_mixed_tol) && !is.null(subst <- substitute(optim_mixed_tol))) {
    u <- try(is.numeric(optim_mixed_tol) || is.logical(optim_mixed_tol) || is.language(optim_mixed_tol)
	 || is.list(optim_mixed_tol) || is(optim_mixed_tol, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['optim_mixed_tol']] <- optim_mixed_tol
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['optim_mixed_tol']] <- optim_mixed_tol
    else  stop('random parameter not allowed')
  }
  if (hasArg(optim_mixed_maxpo) && !is.null(subst <- substitute(optim_mixed_maxpo))) {
    u <- try(is.numeric(optim_mixed_maxpo) || is.logical(optim_mixed_maxpo) || is.language(optim_mixed_maxpo)
	 || is.list(optim_mixed_maxpo) || is(optim_mixed_maxpo, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['optim_mixed_maxpo']] <- optim_mixed_maxpo
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['optim_mixed_maxpo']] <- optim_mixed_maxpo
    else  stop('random parameter not allowed')
  }
  if (hasArg(lambda) && !is.null(subst <- substitute(lambda))) {
    u <- try(is.numeric(lambda) || is.logical(lambda) || is.language(lambda)
	 || is.list(lambda) || is(lambda, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['lambda']] <- lambda
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['lambda']] <- lambda
    else  stop('random parameter not allowed')
  }
  if (hasArg(areamat) && !is.null(subst <- substitute(areamat))) {
    u <- try(is.numeric(areamat) || is.logical(areamat) || is.language(areamat)
	 || is.list(areamat) || is(areamat, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['areamat']] <- areamat
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['areamat']] <- areamat
    else  stop('random parameter not allowed')
  }
  if (hasArg(variobound) && !is.null(subst <- substitute(variobound))) {
    u <- try(is.numeric(variobound) || is.logical(variobound) || is.language(variobound)
	 || is.list(variobound) || is(variobound, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['variobound']] <- variobound
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['variobound']] <- variobound
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPbrmixed', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrmixed <- new('RMmodelgenerator',
	.Data = RPbrmixed,
	type = 'method for Brown-Resnick processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPbrshifted <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg(xi) && !is.null(subst <- substitute(xi))) {
    u <- try(is.numeric(xi) || is.logical(xi) || is.language(xi)
	 || is.list(xi) || is(xi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['xi']] <- xi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['xi']] <- xi
    else  stop('random parameter not allowed')
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPbrshifted', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrshifted <- new('RMmodelgenerator',
	.Data = RPbrshifted,
	type = 'method for Brown-Resnick processes',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPbernoulli <- function(phi, stationary_only, threshold) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(stationary_only) && !is.null(subst <- substitute(stationary_only))) {
    u <- try(is.numeric(stationary_only) || is.logical(stationary_only) || is.language(stationary_only)
	 || is.list(stationary_only) || is(stationary_only, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['stationary_only']] <- stationary_only
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['stationary_only']] <- stationary_only
    else  stop('random parameter not allowed')
  }
  if (hasArg(threshold) && !is.null(subst <- substitute(threshold))) {
    u <- try(is.numeric(threshold) || is.logical(threshold) || is.language(threshold)
	 || is.list(threshold) || is(threshold, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['threshold']] <- threshold
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['threshold']] <- threshold
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPbernoulli', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbernoulli <- new('RMmodelgenerator',
	.Data = RPbernoulli,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 11000,
	vdim = -3
	)



RPbrownresnick <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg(xi) && !is.null(subst <- substitute(xi))) {
    u <- try(is.numeric(xi) || is.logical(xi) || is.language(xi)
	 || is.list(xi) || is(xi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['xi']] <- xi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['xi']] <- xi
    else  stop('random parameter not allowed')
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPbrownresnick', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrownresnick <- new('RMmodelgenerator',
	.Data = RPbrownresnick,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPgauss <- function(phi, stationary_only) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(stationary_only) && !is.null(subst <- substitute(stationary_only))) {
    u <- try(is.numeric(stationary_only) || is.logical(stationary_only) || is.language(stationary_only)
	 || is.list(stationary_only) || is(stationary_only, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['stationary_only']] <- stationary_only
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['stationary_only']] <- stationary_only
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPgauss <- new('RMmodelgenerator',
	.Data = RPgauss,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 11000,
	vdim = -3
	)



RPpoisson <- function(phi, intensity) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(intensity) && !is.null(subst <- substitute(intensity))) {
    u <- try(is.numeric(intensity) || is.logical(intensity) || is.language(intensity)
	 || is.list(intensity) || is(intensity, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['intensity']] <- intensity
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['intensity']] <- intensity
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPpoisson', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPpoisson <- new('RMmodelgenerator',
	.Data = RPpoisson,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPschlather <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg(xi) && !is.null(subst <- substitute(xi))) {
    u <- try(is.numeric(xi) || is.logical(xi) || is.language(xi)
	 || is.list(xi) || is(xi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['xi']] <- xi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['xi']] <- xi
    else  stop('random parameter not allowed')
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPschlather', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPschlather <- new('RMmodelgenerator',
	.Data = RPschlather,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 11000,
	vdim = 1
	)



RPopitz <- function(phi, xi, mu, s, alpha) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(xi) && !is.null(subst <- substitute(xi))) {
    u <- try(is.numeric(xi) || is.logical(xi) || is.language(xi)
	 || is.list(xi) || is(xi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['xi']] <- xi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['xi']] <- xi
    else  stop('random parameter not allowed')
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else  stop('random parameter not allowed')
  }
  if (hasArg(alpha) && !is.null(subst <- substitute(alpha))) {
    u <- try(is.numeric(alpha) || is.logical(alpha) || is.language(alpha)
	 || is.list(alpha) || is(alpha, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['alpha']] <- alpha
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['alpha']] <- alpha
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPopitz', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPopitz <- new('RMmodelgenerator',
	.Data = RPopitz,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 11000,
	vdim = 1
	)



RPsmith <- function(shape, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(shape)) submodels[['shape']] <- shape
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg(xi) && !is.null(subst <- substitute(xi))) {
    u <- try(is.numeric(xi) || is.logical(xi) || is.language(xi)
	 || is.list(xi) || is(xi, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['xi']] <- xi
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['xi']] <- xi
    else  stop('random parameter not allowed')
  }
  if (hasArg(mu) && !is.null(subst <- substitute(mu))) {
    u <- try(is.numeric(mu) || is.logical(mu) || is.language(mu)
	 || is.list(mu) || is(mu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['mu']] <- mu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['mu']] <- mu
    else  stop('random parameter not allowed')
  }
  if (hasArg(s) && !is.null(subst <- substitute(s))) {
    u <- try(is.numeric(s) || is.logical(s) || is.language(s)
	 || is.list(s) || is(s, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['s']] <- s
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['s']] <- s
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPsmith', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPsmith <- new('RMmodelgenerator',
	.Data = RPsmith,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPchi2 <- function(phi, f) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(f) && !is.null(subst <- substitute(f))) {
    u <- try(is.numeric(f) || is.logical(f) || is.language(f)
	 || is.list(f) || is(f, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['f']] <- f
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['f']] <- f
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPchi2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPchi2 <- new('RMmodelgenerator',
	.Data = RPchi2,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 11000,
	vdim = -3
	)



RPt <- function(phi, nu) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg(nu) && !is.null(subst <- substitute(nu))) {
    u <- try(is.numeric(nu) || is.logical(nu) || is.language(nu)
	 || is.list(nu) || is(nu, class2='RMmodel'), silent=TRUE)
    if (is.logical(u) && u) par.model[['nu']] <- nu
    else if (substr(deparse(subst), 1, 1)=='R') par.model[['nu']] <- nu
    else  stop('random parameter not allowed')
  }
  
  model <- new('RMmodel', call = cl, name = 'RPt', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPt <- new('RMmodelgenerator',
	.Data = RPt,
	type = 'process',
	domain = 'single variable',
	isotropy = 'non-dimension-reducing',
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 11000,
	vdim = -3
	)



