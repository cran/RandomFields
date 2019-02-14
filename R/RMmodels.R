
## This file has been created automatically by 'rfGenerateModels'.


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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMplus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMplus <- new(CLASS_RM, 
	.Data = RMplus,
	type = c('of manifold type'),
	isotropy = c('submodel dependent'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMmult', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmult <- new(CLASS_RM, 
	.Data = RMmult,
	type = c('of manifold type'),
	isotropy = c('submodel dependent'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMS  <- function(phi, var, scale, Aniso, proj, anisoT) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.model[['var']] <- CheckArg(var, subst, TRUE)
  if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.model[['scale']] <- CheckArg(scale, subst, TRUE)
  if (hasArg('anisoT') && !is.null(subst <- substitute(anisoT))) 
	par.model[['anisoT']] <- CheckArg(anisoT, subst, TRUE)
  if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.model[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
  if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.model[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMS', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMS <- new(CLASS_RM, 
	.Data = RMS,
	type = c('of manifold type', 'of manifold type'),
	isotropy = c('submodel dependent', 'submodel dependent'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMave <- function(phi, A, z, spacetime, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  if (hasArg('z') && !is.null(subst <- substitute(z))) 
	par.model[['z']] <- CheckArg(z, subst, TRUE)
  if (hasArg('spacetime') && !is.null(subst <- substitute(spacetime))) 
	par.model[['spacetime']] <- CheckArg(spacetime, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMave', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMave <- new(CLASS_RM, 
	.Data = RMave,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 1
	)



RMbcw <- function(alpha, beta, c, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckArg(c, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbcw', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbcw <- new(CLASS_RM, 
	.Data = RMbcw,
	type = c('variogram', 'positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMlsfbm <- function(alpha, const, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('const') && !is.null(subst <- substitute(const))) 
	par.model[['const']] <- CheckArg(const, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMlsfbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMlsfbm <- new(CLASS_RM, 
	.Data = RMlsfbm,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMbessel <- function(nu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbessel', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbessel <- new(CLASS_RM, 
	.Data = RMbessel,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('kappa') && !is.null(subst <- substitute(kappa))) 
	par.model[['kappa']] <- CheckArg(kappa, subst, TRUE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, TRUE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
  if (hasArg('sred12') && !is.null(subst <- substitute(sred12))) 
	par.model[['sred12']] <- CheckArg(sred12, subst, TRUE)
  if (hasArg('gamma') && !is.null(subst <- substitute(gamma))) 
	par.model[['gamma']] <- CheckArg(gamma, subst, TRUE)
  if (hasArg('cdiag') && !is.null(subst <- substitute(cdiag))) 
	par.model[['cdiag']] <- CheckArg(cdiag, subst, TRUE)
  if (hasArg('rhored') && !is.null(subst <- substitute(rhored))) 
	par.model[['rhored']] <- CheckArg(rhored, subst, TRUE)
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckArg(c, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbigneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbigneiting <- new(CLASS_RM, 
	.Data = RMbigneiting,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('threshold') && !is.null(subst <- substitute(threshold))) 
	par.model[['threshold']] <- CheckArg(threshold, subst, TRUE)
  if (hasArg('correlation') && !is.null(subst <- substitute(correlation))) 
	par.model[['correlation']] <- CheckArg(correlation, subst, TRUE)
  if (hasArg('centred') && !is.null(subst <- substitute(centred))) 
	par.model[['centred']] <- CheckArg(centred, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbernoulli', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbernoulli <- new(CLASS_RM, 
	.Data = RMbernoulli,
	type = c('tail correlation'),
	isotropy = c('submodel dependent'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbiwm <- function(nudiag, nured12, nu, s, cdiag, rhored, c, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nudiag') && !is.null(subst <- substitute(nudiag))) 
	par.model[['nudiag']] <- CheckArg(nudiag, subst, TRUE)
  if (hasArg('nured12') && !is.null(subst <- substitute(nured12))) 
	par.model[['nured12']] <- CheckArg(nured12, subst, TRUE)
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
  if (hasArg('cdiag') && !is.null(subst <- substitute(cdiag))) 
	par.model[['cdiag']] <- CheckArg(cdiag, subst, TRUE)
  if (hasArg('rhored') && !is.null(subst <- substitute(rhored))) 
	par.model[['rhored']] <- CheckArg(rhored, subst, TRUE)
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckArg(c, subst, TRUE)
  if (hasArg('notinvnu') && !is.null(subst <- substitute(notinvnu))) 
	par.model[['notinvnu']] <- CheckArg(notinvnu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbiwm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbiwm <- new(CLASS_RM, 
	.Data = RMbiwm,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 2
	)



RMbistable <- function(alpha, s, cdiag, rho, rhored, betared, alphadiag, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
  if (hasArg('cdiag') && !is.null(subst <- substitute(cdiag))) 
	par.model[['cdiag']] <- CheckArg(cdiag, subst, TRUE)
  if (hasArg('rho') && !is.null(subst <- substitute(rho))) 
	par.model[['rho']] <- CheckArg(rho, subst, TRUE)
  if (hasArg('rhored') && !is.null(subst <- substitute(rhored))) 
	par.model[['rhored']] <- CheckArg(rhored, subst, TRUE)
  if (hasArg('betared') && !is.null(subst <- substitute(betared))) 
	par.model[['betared']] <- CheckArg(betared, subst, TRUE)
  if (hasArg('alphadiag') && !is.null(subst <- substitute(alphadiag))) 
	par.model[['alphadiag']] <- CheckArg(alphadiag, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbistable', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbistable <- new(CLASS_RM, 
	.Data = RMbistable,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 2
	)



RMblend <- function(multi, blend, thresholds, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(multi)) submodels[['multi']] <- multi
  if (hasArg(blend)) submodels[['blend']] <- blend
  
  if (hasArg('thresholds') && !is.null(subst <- substitute(thresholds))) 
	par.model[['thresholds']] <- CheckArg(thresholds, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMblend', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMblend <- new(CLASS_RM, 
	.Data = RMblend,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbrownresnick <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbrownresnick', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbrownresnick <- new(CLASS_RM, 
	.Data = RMbrownresnick,
	type = c('tail correlation'),
	isotropy = c('submodel dependent'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbr2bg', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbr2bg <- new(CLASS_RM, 
	.Data = RMbr2bg,
	type = c('positive definite'),
	isotropy = c('submodel dependent'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbr2eg', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbr2eg <- new(CLASS_RM, 
	.Data = RMbr2eg,
	type = c('positive definite'),
	isotropy = c('submodel dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMbubble <- function(phi, scaling, z, weight, minscale, barycentre, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(scaling)) submodels[['scaling']] <- scaling
  
  if (hasArg('z') && !is.null(subst <- substitute(z))) 
	par.model[['z']] <- CheckArg(z, subst, TRUE)
  if (hasArg('weight') && !is.null(subst <- substitute(weight))) 
	par.model[['weight']] <- CheckArg(weight, subst, TRUE)
  if (hasArg('minscale') && !is.null(subst <- substitute(minscale))) 
	par.model[['minscale']] <- CheckArg(minscale, subst, TRUE)
  if (hasArg('barycentre') && !is.null(subst <- substitute(barycentre))) 
	par.model[['barycentre']] <- CheckArg(barycentre, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbubble', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbubble <- new(CLASS_RM, 
	.Data = RMbubble,
	type = c('positive definite'),
	isotropy = c('non-dimension-reducing'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMcauchy <- function(gamma, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('gamma') && !is.null(subst <- substitute(gamma))) 
	par.model[['gamma']] <- CheckArg(gamma, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMcauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcauchy <- new(CLASS_RM, 
	.Data = RMcauchy,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMcircular', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcircular <- new(CLASS_RM, 
	.Data = RMcircular,
	type = c('tail correlation'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'Gneiting-Schaback class',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMconstant <- function(M, var) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
  model <- methods::new('RMmodel', call = cl, name = 'RMconstant', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMconstant <- new(CLASS_RM, 
	.Data = RMconstant,
	type = c('positive definite', 'negative definite'),
	isotropy = c('framework dependent', 'isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'submodel dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



iRMfixcov <- function(norm, M, x, raw, var, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(norm)) submodels[['norm']] <- norm
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckArg(x, subst, TRUE)
  if (hasArg('raw') && !is.null(subst <- substitute(raw))) 
	par.model[['raw']] <- CheckArg(raw, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMfixcov', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRMfixcov <- new(CLASS_RM, 
	.Data = iRMfixcov,
	type = c('positive definite', 'positive definite', 'positive definite'),
	isotropy = c('non-dimension-reducing', 'isotropic', 'earth isotropic'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = -1
	)



RMcoxisham <- function(phi, mu, D, beta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, TRUE)
  if (hasArg('D') && !is.null(subst <- substitute(D))) 
	par.model[['D']] <- CheckArg(D, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMcoxisham', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcoxisham <- new(CLASS_RM, 
	.Data = RMcoxisham,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMcubic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcubic <- new(CLASS_RM, 
	.Data = RMcubic,
	type = c('tail correlation'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMcurlfree <- function(phi, which, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('which') && !is.null(subst <- substitute(which))) 
	par.model[['which']] <- CheckArg(which, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMcurlfree', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcurlfree <- new(CLASS_RM, 
	.Data = RMcurlfree,
	type = c('positive definite'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMcutoff <- function(phi, diameter, a, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, TRUE)
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMcutoff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMcutoff <- new(CLASS_RM, 
	.Data = RMcutoff,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = -3
	)



RMdagum <- function(beta, gamma, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
  if (hasArg('gamma') && !is.null(subst <- substitute(gamma))) 
	par.model[['gamma']] <- CheckArg(gamma, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMdagum', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdagum <- new(CLASS_RM, 
	.Data = RMdagum,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'parameter dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMdampedcos <- function(lambda, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMdampedcos', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdampedcos <- new(CLASS_RM, 
	.Data = RMdampedcos,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMderiv <- function(phi, which, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('which') && !is.null(subst <- substitute(which))) 
	par.model[['which']] <- CheckArg(which, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMderiv', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMderiv <- new(CLASS_RM, 
	.Data = RMderiv,
	type = c('positive definite'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMdewijsian <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMdewijsian', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdewijsian <- new(CLASS_RM, 
	.Data = RMdewijsian,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMdivfree <- function(phi, which, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('which') && !is.null(subst <- substitute(which))) 
	par.model[['which']] <- CheckArg(which, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMdivfree', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdivfree <- new(CLASS_RM, 
	.Data = RMdivfree,
	type = c('positive definite'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMepscauchy <- function(alpha, beta, eps, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
  if (hasArg('eps') && !is.null(subst <- substitute(eps))) 
	par.model[['eps']] <- CheckArg(eps, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMepscauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMepscauchy <- new(CLASS_RM, 
	.Data = RMepscauchy,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMexp <- new(CLASS_RM, 
	.Data = RMexp,
	type = c('tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('n') && !is.null(subst <- substitute(n))) 
	par.model[['n']] <- CheckArg(n, subst, TRUE)
  if (hasArg('standardised') && !is.null(subst <- substitute(standardised))) 
	par.model[['standardised']] <- CheckArg(standardised, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMexponential', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMexponential <- new(CLASS_RM, 
	.Data = RMexponential,
	type = c('positive definite'),
	isotropy = c('submodel dependent'),
	domain = c('submodel dependent'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMschlather', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMschlather <- new(CLASS_RM, 
	.Data = RMschlather,
	type = c('tail correlation'),
	isotropy = c('submodel dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMfractdiff <- function(a, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMfractdiff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfractdiff <- new(CLASS_RM, 
	.Data = RMfractdiff,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



RMflatpower <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMflatpower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMflatpower <- new(CLASS_RM, 
	.Data = RMflatpower,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'Bernstein',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMfbm <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMfbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfbm <- new(CLASS_RM, 
	.Data = RMfbm,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMfractgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMfractgauss <- new(CLASS_RM, 
	.Data = RMfractgauss,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgauss <- new(CLASS_RM, 
	.Data = RMgauss,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMgenfbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgenfbm <- new(CLASS_RM, 
	.Data = RMgenfbm,
	type = c('variogram'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMgencauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgencauchy <- new(CLASS_RM, 
	.Data = RMgencauchy,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'parameter dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMbicauchy <- function(alpha, beta, s, rho, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
  if (hasArg('rho') && !is.null(subst <- substitute(rho))) 
	par.model[['rho']] <- CheckArg(rho, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMbicauchy', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMbicauchy <- new(CLASS_RM, 
	.Data = RMbicauchy,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 2
	)



RMgengneiting <- function(kappa, mu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('kappa') && !is.null(subst <- substitute(kappa))) 
	par.model[['kappa']] <- CheckArg(kappa, subst, TRUE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMgengneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgengneiting <- new(CLASS_RM, 
	.Data = RMgengneiting,
	type = c('positive definite', 'positive definite', 'positive definite', 'positive definite', 'positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic', 'isotropic', 'spherical isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('orig') && !is.null(subst <- substitute(orig))) 
	par.model[['orig']] <- CheckArg(orig, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMgneiting', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgneiting <- new(CLASS_RM, 
	.Data = RMgneiting,
	type = c('positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = 1
	)



RMgennsst <- function(phi, psi, dim_u, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(psi)) submodels[['psi']] <- psi
  
  if (hasArg('dim_u') && !is.null(subst <- substitute(dim_u))) 
	par.model[['dim_u']] <- CheckArg(dim_u, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMgennsst', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMgennsst <- new(CLASS_RM, 
	.Data = RMgennsst,
	type = c('positive definite', 'positive definite'),
	isotropy = c('symmetric', 'symmetric'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMhyperbolic <- function(nu, lambda, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMhyperbolic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMhyperbolic <- new(CLASS_RM, 
	.Data = RMhyperbolic,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMiaco', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMiaco <- new(CLASS_RM, 
	.Data = RMiaco,
	type = c('positive definite'),
	isotropy = c('space-isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMidmodel <- function(phi, vdim, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('vdim') && !is.null(subst <- substitute(vdim))) 
	par.model[['vdim']] <- CheckArg(vdim, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMidmodel', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMidmodel <- new(CLASS_RM, 
	.Data = RMidmodel,
	type = c('of manifold type'),
	isotropy = c('framework dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMkolmogorov <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMkolmogorov', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMkolmogorov <- new(CLASS_RM, 
	.Data = RMkolmogorov,
	type = c('variogram'),
	isotropy = c('vector-isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('beta') && !is.null(subst <- substitute(beta))) 
	par.model[['beta']] <- CheckArg(beta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMlgd', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMlgd <- new(CLASS_RM, 
	.Data = RMlgd,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMmastein', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmastein <- new(CLASS_RM, 
	.Data = RMmastein,
	type = c('positive definite'),
	isotropy = c('space-isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  if (hasArg('theta') && !is.null(subst <- substitute(theta))) 
	par.model[['theta']] <- CheckArg(theta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMma', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMma <- new(CLASS_RM, 
	.Data = RMma,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMintexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMintexp <- new(CLASS_RM, 
	.Data = RMintexp,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMmatrix <- function(C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, M, vdim, var, scale, Aniso, proj) {
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
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
  if (hasArg('vdim') && !is.null(subst <- substitute(vdim))) 
	par.model[['vdim']] <- CheckArg(vdim, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMmatrix', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmatrix <- new(CLASS_RM, 
	.Data = RMmatrix,
	type = c('of manifold type'),
	isotropy = c('submodel dependent'),
	domain = c('submodel dependent'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMmatern <- function(nu, notinvnu, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('notinvnu') && !is.null(subst <- substitute(notinvnu))) 
	par.model[['notinvnu']] <- CheckArg(notinvnu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMmatern', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmatern <- new(CLASS_RM, 
	.Data = RMmatern,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('parameter dependent'),
	operator = FALSE,
	monotone = 'submodel dependent monotonicity',
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
  
  if (hasArg('theta') && !is.null(subst <- substitute(theta))) 
	par.model[['theta']] <- CheckArg(theta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMmqam', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmqam <- new(CLASS_RM, 
	.Data = RMmqam,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMmultiquad <- function(delta, tau, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
  if (hasArg('tau') && !is.null(subst <- substitute(tau))) 
	par.model[['tau']] <- CheckArg(tau, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMmultiquad', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmultiquad <- new(CLASS_RM, 
	.Data = RMmultiquad,
	type = c('positive definite'),
	isotropy = c('spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMnatsc <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMnatsc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnatsc <- new(CLASS_RM, 
	.Data = RMnatsc,
	type = c('positive definite', 'tail correlation'),
	isotropy = c('isotropic', 'isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMnsst <- function(phi, psi, delta, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(psi)) submodels[['psi']] <- psi
  
  if (hasArg('delta') && !is.null(subst <- substitute(delta))) 
	par.model[['delta']] <- CheckArg(delta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMnsst', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnsst <- new(CLASS_RM, 
	.Data = RMnsst,
	type = c('positive definite'),
	isotropy = c('space-isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMparswm <- function(nudiag, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('nudiag') && !is.null(subst <- substitute(nudiag))) 
	par.model[['nudiag']] <- CheckArg(nudiag, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMparswm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMparswm <- new(CLASS_RM, 
	.Data = RMparswm,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMpenta', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpenta <- new(CLASS_RM, 
	.Data = RMpenta,
	type = c('positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMaskey', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMaskey <- new(CLASS_RM, 
	.Data = RMaskey,
	type = c('positive definite', 'positive definite', 'tail correlation'),
	isotropy = c('isotropic', 'spherical isotropic', 'isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMpower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpower <- new(CLASS_RM, 
	.Data = RMpower,
	type = c('shape function', 'negative definite', 'positive definite', 'tail correlation'),
	isotropy = c('submodel dependent', 'submodel dependent', 'submodel dependent', 'isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMprod <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMprod', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMprod <- new(CLASS_RM, 
	.Data = RMprod,
	type = c('positive definite'),
	isotropy = c('non-dimension-reducing'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
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
  
  if (hasArg('theta') && !is.null(subst <- substitute(theta))) 
	par.model[['theta']] <- CheckArg(theta, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMqam', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMqam <- new(CLASS_RM, 
	.Data = RMqam,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMqexp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMqexp <- new(CLASS_RM, 
	.Data = RMqexp,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMscale <- function(phi, scaling, penalty, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(scaling)) submodels[['scaling']] <- scaling
  if (hasArg(penalty)) submodels[['penalty']] <- penalty
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMscale', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMscale <- new(CLASS_RM, 
	.Data = RMscale,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('kernel'),
	operator = TRUE,
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
  
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
  if (hasArg('diag') && !is.null(subst <- substitute(diag))) 
	par.model[['diag']] <- CheckArg(diag, subst, TRUE)
  if (hasArg('rhored') && !is.null(subst <- substitute(rhored))) 
	par.model[['rhored']] <- CheckArg(rhored, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMschur', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMschur <- new(CLASS_RM, 
	.Data = RMschur,
	type = c('positive definite'),
	isotropy = c('framework dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RMdelay <- function(phi, s, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMdelay', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMdelay <- new(CLASS_RM, 
	.Data = RMdelay,
	type = c('positive definite'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMsinepower <- function(alpha, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMsinepower', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMsinepower <- new(CLASS_RM, 
	.Data = RMsinepower,
	type = c('positive definite'),
	isotropy = c('spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RMspheric <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMspheric', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMspheric <- new(CLASS_RM, 
	.Data = RMspheric,
	type = c('tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMstable', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstable <- new(CLASS_RM, 
	.Data = RMstable,
	type = c('positive definite', 'tail correlation', 'positive definite'),
	isotropy = c('isotropic', 'isotropic', 'spherical isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'parameter dependent monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMintrinsic <- function(phi, diameter, rawR, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, TRUE)
  if (hasArg('rawR') && !is.null(subst <- substitute(rawR))) 
	par.model[['rawR']] <- CheckArg(rawR, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMintrinsic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMintrinsic <- new(CLASS_RM, 
	.Data = RMintrinsic,
	type = c('positive definite', 'positive definite'),
	isotropy = c('isotropic', 'spherical isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('z') && !is.null(subst <- substitute(z))) 
	par.model[['z']] <- CheckArg(z, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMstein', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstein <- new(CLASS_RM, 
	.Data = RMstein,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('single variable'),
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
  
  if (hasArg('S') && !is.null(subst <- substitute(S))) 
	par.model[['S']] <- CheckArg(S, subst, TRUE)
  if (hasArg('z') && !is.null(subst <- substitute(z))) 
	par.model[['z']] <- CheckArg(z, subst, TRUE)
  if (hasArg('M') && !is.null(subst <- substitute(M))) 
	par.model[['M']] <- CheckArg(M, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMstp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMstp <- new(CLASS_RM, 
	.Data = RMstp,
	type = c('positive definite'),
	isotropy = c('symmetric'),
	domain = c('kernel'),
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
  
  if (hasArg('fulldim') && !is.null(subst <- substitute(fulldim))) 
	par.model[['fulldim']] <- CheckArg(fulldim, subst, TRUE)
  if (hasArg('reduceddim') && !is.null(subst <- substitute(reduceddim))) 
	par.model[['reduceddim']] <- CheckArg(reduceddim, subst, TRUE)
  if (hasArg('layers') && !is.null(subst <- substitute(layers))) 
	par.model[['layers']] <- CheckArg(layers, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMtbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtbm <- new(CLASS_RM, 
	.Data = RMtbm,
	type = c('of manifold type'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -1,
	vdim = -3
	)



RMsum <- function(phi, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMsum', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMsum <- new(CLASS_RM, 
	.Data = RMsum,
	type = c('negative definite'),
	isotropy = c('non-dimension-reducing'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



iRMcov <- function(gamma, x, a, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(gamma)) submodels[['gamma']] <- gamma
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckArg(x, subst, TRUE)
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMcov', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRMcov <- new(CLASS_RM, 
	.Data = iRMcov,
	type = c('positive definite'),
	isotropy = c('cartesian system'),
	domain = c('kernel'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = 1
	)



RMvector <- function(phi, a, Dspace, var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  if (hasArg('Dspace') && !is.null(subst <- substitute(Dspace))) 
	par.model[['Dspace']] <- CheckArg(Dspace, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMvector', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMvector <- new(CLASS_RM, 
	.Data = RMvector,
	type = c('positive definite', 'positive definite'),
	isotropy = c('cartesian system', 'cartesian system'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RMwave <- function(var, scale, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMwave', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMwave <- new(CLASS_RM, 
	.Data = RMwave,
	type = c('positive definite'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, TRUE)
  if (hasArg('notinvnu') && !is.null(subst <- substitute(notinvnu))) 
	par.model[['notinvnu']] <- CheckArg(notinvnu, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMwhittle', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMwhittle <- new(CLASS_RM, 
	.Data = RMwhittle,
	type = c('positive definite'),
	isotropy = c('parameter dependent'),
	domain = c('parameter dependent'),
	operator = FALSE,
	monotone = 'normal mixture',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMnugget <- function(tol, vdim, var, Aniso, proj) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('tol') && !is.null(subst <- substitute(tol))) 
	par.model[['tol']] <- CheckArg(tol, subst, TRUE)
  if (hasArg('vdim') && !is.null(subst <- substitute(vdim))) 
	par.model[['vdim']] <- CheckArg(vdim, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMnugget', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMnugget <- new(CLASS_RM, 
	.Data = RMnugget,
	type = c('tail correlation'),
	isotropy = c('parameter dependent'),
	domain = c('parameter dependent'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -2
	)



RMtrend <- function(mean) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('mean') && !is.null(subst <- substitute(mean))) 
	par.model[['mean']] <- CheckArg(mean, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMtrend', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtrend <- new(CLASS_RM, 
	.Data = RMtrend,
	type = c('trend'),
	isotropy = c('parameter dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMangle <- function(angle, lat.angle, ratio, diag) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('angle') && !is.null(subst <- substitute(angle))) 
	par.model[['angle']] <- CheckArg(angle, subst, TRUE)
  if (hasArg('lat.angle') && !is.null(subst <- substitute(lat.angle))) 
	par.model[['lat.angle']] <- CheckArg(lat.angle, subst, TRUE)
  if (hasArg('ratio') && !is.null(subst <- substitute(ratio))) 
	par.model[['ratio']] <- CheckArg(ratio, subst, TRUE)
  if (hasArg('diag') && !is.null(subst <- substitute(diag))) 
	par.model[['diag']] <- CheckArg(diag, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMangle', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMangle <- new(CLASS_RM, 
	.Data = RMangle,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
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
  
  if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
    if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.general[['scale']] <- CheckArg(scale, subst, TRUE)
    if (hasArg('Aniso') && !is.null(subst <- substitute(Aniso))) 
	par.general[['Aniso']] <- CheckArg(Aniso, subst, TRUE)
    if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.general[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  model <- methods::new('RMmodel', call = cl, name = 'RMball', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMball <- new(CLASS_RM, 
	.Data = RMball,
	type = c('shape function'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'monotone',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



iRMcovariate <- function(norm, data, x, raw, addNA, factor, var) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(norm)) submodels[['norm']] <- norm
  
  if (hasArg('data') && !is.null(subst <- substitute(data))) 
	par.model[['data']] <- CheckArg(data, subst, TRUE)
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckArg(x, subst, TRUE)
  if (hasArg('raw') && !is.null(subst <- substitute(raw))) 
	par.model[['raw']] <- CheckArg(raw, subst, TRUE)
  if (hasArg('addNA') && !is.null(subst <- substitute(addNA))) 
	par.model[['addNA']] <- CheckArg(addNA, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckArg(factor, subst, TRUE)
    if (hasArg('var') && !is.null(subst <- substitute(var))) 
	par.general[['var']] <- CheckArg(var, subst, TRUE)
  model <- methods::new('RMmodel', call = cl, name = 'RMcovariate', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRMcovariate <- new(CLASS_RM, 
	.Data = iRMcovariate,
	type = c('shape function', 'shape function', 'shape function', 'trend'),
	isotropy = c('non-dimension-reducing', 'isotropic', 'earth isotropic', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = FALSE,
	maxdim = Inf,
	vdim = -1
	)



RMeaxxa <- function(E, A) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('E') && !is.null(subst <- substitute(E))) 
	par.model[['E']] <- CheckArg(E, subst, TRUE)
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMeaxxa', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMeaxxa <- new(CLASS_RM, 
	.Data = RMeaxxa,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
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
  
  if (hasArg('E') && !is.null(subst <- substitute(E))) 
	par.model[['E']] <- CheckArg(E, subst, TRUE)
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMetaxxa', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMetaxxa <- new(CLASS_RM, 
	.Data = RMetaxxa,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 10,
	vdim = 3
	)



RMid <- function() {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  

  model <- methods::new('RMmodel', call = cl, name = 'RMid', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMid <- new(CLASS_RM, 
	.Data = RMid,
	type = c('shape function'),
	isotropy = c('framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMtrafo <- function(phi, new) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('new') && !is.null(subst <- substitute(new))) 
	par.model[['new']] <- CheckChar(new, subst, ISO_NAMES, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMtrafo', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtrafo <- new(CLASS_RM, 
	.Data = RMtrafo,
	type = c('of manifold type'),
	isotropy = c('parameter dependent'),
	domain = c('parameter dependent'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RMpolygon <- function(lambda) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMpolygon', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMpolygon <- new(CLASS_RM, 
	.Data = RMpolygon,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
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
  
  if (hasArg('A') && !is.null(subst <- substitute(A))) 
	par.model[['A']] <- CheckArg(A, subst, TRUE)
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMrational', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrational <- new(CLASS_RM, 
	.Data = RMrational,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
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
  
  if (hasArg('speed') && !is.null(subst <- substitute(speed))) 
	par.model[['speed']] <- CheckArg(speed, subst, TRUE)
  if (hasArg('phi') && !is.null(subst <- substitute(phi))) 
	par.model[['phi']] <- CheckArg(phi, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMrotat', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrotat <- new(CLASS_RM, 
	.Data = RMrotat,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
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
  
  if (hasArg('phi') && !is.null(subst <- substitute(phi))) 
	par.model[['phi']] <- CheckArg(phi, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMrotation', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMrotation <- new(CLASS_RM, 
	.Data = RMrotation,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
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
  
  if (hasArg('p') && !is.null(subst <- substitute(p))) 
	par.model[['p']] <- CheckArg(p, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMsign', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMsign <- new(CLASS_RM, 
	.Data = RMsign,
	type = c('shape function'),
	isotropy = c('framework dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RMm2r <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- methods::new('RMmodel', call = cl, name = 'RMm2r', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMm2r <- new(CLASS_RM, 
	.Data = RMm2r,
	type = c('shape function'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = 3,
	vdim = 1
	)



RMm3b <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- methods::new('RMmodel', call = cl, name = 'RMm3b', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMm3b <- new(CLASS_RM, 
	.Data = RMm3b,
	type = c('shape function'),
	isotropy = c('isotropic'),
	domain = c('single variable'),
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
  

  model <- methods::new('RMmodel', call = cl, name = 'RMmps', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmps <- new(CLASS_RM, 
	.Data = RMmps,
	type = c('shape function'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RMtruncsupport <- function(phi, radius) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('radius') && !is.null(subst <- substitute(radius))) 
	par.model[['radius']] <- CheckArg(radius, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMtruncsupport', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMtruncsupport <- new(CLASS_RM, 
	.Data = RMtruncsupport,
	type = c('shape function'),
	isotropy = c('framework dependent'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'submodel dependent monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = 1
	)



RRdeterm <- function(mean) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('mean') && !is.null(subst <- substitute(mean))) 
	par.model[['mean']] <- CheckArg(mean, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RRdeterm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRdeterm <- new(CLASS_RM, 
	.Data = RRdeterm,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = FALSE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RRgauss <- function(mu, sd, log) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('sd') && !is.null(subst <- substitute(sd))) 
	par.model[['sd']] <- CheckArg(sd, subst, FALSE)
  if (hasArg('log') && !is.null(subst <- substitute(log))) 
	par.model[['log']] <- CheckArg(log, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RRgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRgauss <- new(CLASS_RM, 
	.Data = RRgauss,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
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
  
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('scale') && !is.null(subst <- substitute(scale))) 
	par.model[['scale']] <- CheckArg(scale, subst, FALSE)
  if (hasArg('pow') && !is.null(subst <- substitute(pow))) 
	par.model[['pow']] <- CheckArg(pow, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RRloc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRloc <- new(CLASS_RM, 
	.Data = RRloc,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



RRmcmc <- function(phi, mcmc_n, sigma, normed, maxdensity, rand.loc, gibbs) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('mcmc_n') && !is.null(subst <- substitute(mcmc_n))) 
	par.model[['mcmc_n']] <- CheckArg(mcmc_n, subst, FALSE)
  if (hasArg('sigma') && !is.null(subst <- substitute(sigma))) 
	par.model[['sigma']] <- CheckArg(sigma, subst, FALSE)
  if (hasArg('normed') && !is.null(subst <- substitute(normed))) 
	par.model[['normed']] <- CheckArg(normed, subst, FALSE)
  if (hasArg('maxdensity') && !is.null(subst <- substitute(maxdensity))) 
	par.model[['maxdensity']] <- CheckArg(maxdensity, subst, FALSE)
  if (hasArg('rand.loc') && !is.null(subst <- substitute(rand.loc))) 
	par.model[['rand.loc']] <- CheckArg(rand.loc, subst, FALSE)
  if (hasArg('gibbs') && !is.null(subst <- substitute(gibbs))) 
	par.model[['gibbs']] <- CheckArg(gibbs, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RRmcmc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRmcmc <- new(CLASS_RM, 
	.Data = RRmcmc,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -1
	)



RRrectangular <- function(phi, safety, minsteplen, maxsteps, parts, maxit, innermin, outermax, mcmc_n, normed, approx, onesided) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('safety') && !is.null(subst <- substitute(safety))) 
	par.model[['safety']] <- CheckArg(safety, subst, FALSE)
  if (hasArg('minsteplen') && !is.null(subst <- substitute(minsteplen))) 
	par.model[['minsteplen']] <- CheckArg(minsteplen, subst, FALSE)
  if (hasArg('maxsteps') && !is.null(subst <- substitute(maxsteps))) 
	par.model[['maxsteps']] <- CheckArg(maxsteps, subst, FALSE)
  if (hasArg('parts') && !is.null(subst <- substitute(parts))) 
	par.model[['parts']] <- CheckArg(parts, subst, FALSE)
  if (hasArg('maxit') && !is.null(subst <- substitute(maxit))) 
	par.model[['maxit']] <- CheckArg(maxit, subst, FALSE)
  if (hasArg('innermin') && !is.null(subst <- substitute(innermin))) 
	par.model[['innermin']] <- CheckArg(innermin, subst, FALSE)
  if (hasArg('outermax') && !is.null(subst <- substitute(outermax))) 
	par.model[['outermax']] <- CheckArg(outermax, subst, FALSE)
  if (hasArg('mcmc_n') && !is.null(subst <- substitute(mcmc_n))) 
	par.model[['mcmc_n']] <- CheckArg(mcmc_n, subst, FALSE)
  if (hasArg('normed') && !is.null(subst <- substitute(normed))) 
	par.model[['normed']] <- CheckArg(normed, subst, FALSE)
  if (hasArg('approx') && !is.null(subst <- substitute(approx))) 
	par.model[['approx']] <- CheckArg(approx, subst, FALSE)
  if (hasArg('onesided') && !is.null(subst <- substitute(onesided))) 
	par.model[['onesided']] <- CheckArg(onesided, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RRrectangular', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRrectangular <- new(CLASS_RM, 
	.Data = RRrectangular,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
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
  
  if (hasArg('spacedim') && !is.null(subst <- substitute(spacedim))) 
	par.model[['spacedim']] <- CheckArg(spacedim, subst, FALSE)
  if (hasArg('balldim') && !is.null(subst <- substitute(balldim))) 
	par.model[['balldim']] <- CheckArg(balldim, subst, FALSE)
  if (hasArg('R') && !is.null(subst <- substitute(R))) 
	par.model[['R']] <- CheckArg(R, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RRspheric', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRspheric <- new(CLASS_RM, 
	.Data = RRspheric,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable'),
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
  
  if (hasArg('min') && !is.null(subst <- substitute(min))) 
	par.model[['min']] <- CheckArg(min, subst, FALSE)
  if (hasArg('max') && !is.null(subst <- substitute(max))) 
	par.model[['max']] <- CheckArg(max, subst, FALSE)
  if (hasArg('normed') && !is.null(subst <- substitute(normed))) 
	par.model[['normed']] <- CheckArg(normed, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RRunif', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RRunif <- new(CLASS_RM, 
	.Data = RRunif,
	type = c('distribution family'),
	isotropy = c('cartesian system'),
	domain = c('single variable', 'kernel'),
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
  
  if (hasArg('p') && !is.null(subst <- substitute(p))) 
	par.model[['p']] <- CheckArg(p, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RMmppplus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RMmppplus <- new(CLASS_RM, 
	.Data = RMmppplus,
	type = c('point-shape function'),
	isotropy = c('framework dependent'),
	domain = c('single variable', 'kernel'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



iRFcov <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- methods::new('RMmodel', call = cl, name = 'RFcov', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRFcov <- new(CLASS_RM, 
	.Data = iRFcov,
	type = c('interface'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



iRFcovmatrix <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- methods::new('RMmodel', call = cl, name = 'RFcovmatrix', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRFcovmatrix <- new(CLASS_RM, 
	.Data = iRFcovmatrix,
	type = c('interface'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



iRFloglikelihood <- function(phi, data, estimate_variance, betas_separate, ignore_trend) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('data') && !is.null(subst <- substitute(data))) 
	par.model[['data']] <- CheckArg(data, subst, FALSE)
  if (hasArg('estimate_variance') && !is.null(subst <- substitute(estimate_variance))) 
	par.model[['estimate_variance']] <- CheckArg(estimate_variance, subst, FALSE)
  if (hasArg('betas_separate') && !is.null(subst <- substitute(betas_separate))) 
	par.model[['betas_separate']] <- CheckArg(betas_separate, subst, FALSE)
  if (hasArg('ignore_trend') && !is.null(subst <- substitute(ignore_trend))) 
	par.model[['ignore_trend']] <- CheckArg(ignore_trend, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RFloglikelihood', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRFloglikelihood <- new(CLASS_RM, 
	.Data = iRFloglikelihood,
	type = c('interface'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = FALSE,
	maxdim = -3,
	vdim = -3
	)



iRFpseudovariogra <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- methods::new('RMmodel', call = cl, name = 'RFpseudovariogra', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRFpseudovariogra <- new(CLASS_RM, 
	.Data = iRFpseudovariogra,
	type = c('interface'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



iRFvariogram <- function(phi) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  

  model <- methods::new('RMmodel', call = cl, name = 'RFvariogram', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRFvariogram <- new(CLASS_RM, 
	.Data = iRFvariogram,
	type = c('interface'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -3
	)



iRFsimulate <- function(phi, checkonly, setseed, env) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('checkonly') && !is.null(subst <- substitute(checkonly))) 
	par.model[['checkonly']] <- CheckArg(checkonly, subst, FALSE)
  if (hasArg('setseed') && !is.null(subst <- substitute(setseed))) 
	par.model[['setseed']] <- CheckArg(setseed, subst, FALSE)
  if (hasArg('env') && !is.null(subst <- substitute(env))) 
	par.model[['env']] <- CheckArg(env, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RFsimulate', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

iRFsimulate <- new(CLASS_RM, 
	.Data = iRFsimulate,
	type = c('interface'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = NA,
	simpleArguments = FALSE,
	maxdim = -3,
	vdim = -3
	)



RPtrend <- function(phi, mean) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('mean') && !is.null(subst <- substitute(mean))) 
	par.model[['mean']] <- CheckArg(mean, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPtrend', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPtrend <- new(CLASS_RM, 
	.Data = RPtrend,
	type = c('process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPaverage <- function(phi, shape, boxcox, intensity, method) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(shape)) submodels[['shape']] <- shape
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('intensity') && !is.null(subst <- substitute(intensity))) 
	par.model[['intensity']] <- CheckArg(intensity, subst, FALSE)
  if (hasArg('method') && !is.null(subst <- substitute(method))) 
	par.model[['method']] <- CheckArg(method, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPaverage', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPaverage <- new(CLASS_RM, 
	.Data = RPaverage,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPcoins <- function(phi, shape, boxcox, intensity, method) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(shape)) submodels[['shape']] <- shape
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('intensity') && !is.null(subst <- substitute(intensity))) 
	par.model[['intensity']] <- CheckArg(intensity, subst, FALSE)
  if (hasArg('method') && !is.null(subst <- substitute(method))) 
	par.model[['method']] <- CheckArg(method, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPcoins', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcoins <- new(CLASS_RM, 
	.Data = RPcoins,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPcirculant <- function(phi, boxcox, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('force') && !is.null(subst <- substitute(force))) 
	par.model[['force']] <- CheckArg(force, subst, FALSE)
  if (hasArg('mmin') && !is.null(subst <- substitute(mmin))) 
	par.model[['mmin']] <- CheckArg(mmin, subst, FALSE)
  if (hasArg('strategy') && !is.null(subst <- substitute(strategy))) 
	par.model[['strategy']] <- CheckArg(strategy, subst, FALSE)
  if (hasArg('maxGB') && !is.null(subst <- substitute(maxGB))) 
	par.model[['maxGB']] <- CheckArg(maxGB, subst, FALSE)
  if (hasArg('maxmem') && !is.null(subst <- substitute(maxmem))) 
	par.model[['maxmem']] <- CheckArg(maxmem, subst, FALSE)
  if (hasArg('tolIm') && !is.null(subst <- substitute(tolIm))) 
	par.model[['tolIm']] <- CheckArg(tolIm, subst, FALSE)
  if (hasArg('tolRe') && !is.null(subst <- substitute(tolRe))) 
	par.model[['tolRe']] <- CheckArg(tolRe, subst, FALSE)
  if (hasArg('trials') && !is.null(subst <- substitute(trials))) 
	par.model[['trials']] <- CheckArg(trials, subst, FALSE)
  if (hasArg('useprimes') && !is.null(subst <- substitute(useprimes))) 
	par.model[['useprimes']] <- CheckArg(useprimes, subst, FALSE)
  if (hasArg('dependent') && !is.null(subst <- substitute(dependent))) 
	par.model[['dependent']] <- CheckArg(dependent, subst, FALSE)
  if (hasArg('approx_step') && !is.null(subst <- substitute(approx_step))) 
	par.model[['approx_step']] <- CheckArg(approx_step, subst, FALSE)
  if (hasArg('approx_maxgrid') && !is.null(subst <- substitute(approx_maxgrid))) 
	par.model[['approx_maxgrid']] <- CheckArg(approx_maxgrid, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPcirculant', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcirculant <- new(CLASS_RM, 
	.Data = RPcirculant,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = -3
	)



RPcutoff <- function(phi, boxcox, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid, diameter, a) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('force') && !is.null(subst <- substitute(force))) 
	par.model[['force']] <- CheckArg(force, subst, FALSE)
  if (hasArg('mmin') && !is.null(subst <- substitute(mmin))) 
	par.model[['mmin']] <- CheckArg(mmin, subst, FALSE)
  if (hasArg('strategy') && !is.null(subst <- substitute(strategy))) 
	par.model[['strategy']] <- CheckArg(strategy, subst, FALSE)
  if (hasArg('maxGB') && !is.null(subst <- substitute(maxGB))) 
	par.model[['maxGB']] <- CheckArg(maxGB, subst, FALSE)
  if (hasArg('maxmem') && !is.null(subst <- substitute(maxmem))) 
	par.model[['maxmem']] <- CheckArg(maxmem, subst, FALSE)
  if (hasArg('tolIm') && !is.null(subst <- substitute(tolIm))) 
	par.model[['tolIm']] <- CheckArg(tolIm, subst, FALSE)
  if (hasArg('tolRe') && !is.null(subst <- substitute(tolRe))) 
	par.model[['tolRe']] <- CheckArg(tolRe, subst, FALSE)
  if (hasArg('trials') && !is.null(subst <- substitute(trials))) 
	par.model[['trials']] <- CheckArg(trials, subst, FALSE)
  if (hasArg('useprimes') && !is.null(subst <- substitute(useprimes))) 
	par.model[['useprimes']] <- CheckArg(useprimes, subst, FALSE)
  if (hasArg('dependent') && !is.null(subst <- substitute(dependent))) 
	par.model[['dependent']] <- CheckArg(dependent, subst, FALSE)
  if (hasArg('approx_step') && !is.null(subst <- substitute(approx_step))) 
	par.model[['approx_step']] <- CheckArg(approx_step, subst, FALSE)
  if (hasArg('approx_maxgrid') && !is.null(subst <- substitute(approx_maxgrid))) 
	par.model[['approx_maxgrid']] <- CheckArg(approx_maxgrid, subst, FALSE)
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, FALSE)
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckArg(a, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPcutoff', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPcutoff <- new(CLASS_RM, 
	.Data = RPcutoff,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RPintrinsic <- function(phi, boxcox, force, mmin, strategy, maxGB, maxmem, tolIm, tolRe, trials, useprimes, dependent, approx_step, approx_maxgrid, diameter, rawR) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('force') && !is.null(subst <- substitute(force))) 
	par.model[['force']] <- CheckArg(force, subst, FALSE)
  if (hasArg('mmin') && !is.null(subst <- substitute(mmin))) 
	par.model[['mmin']] <- CheckArg(mmin, subst, FALSE)
  if (hasArg('strategy') && !is.null(subst <- substitute(strategy))) 
	par.model[['strategy']] <- CheckArg(strategy, subst, FALSE)
  if (hasArg('maxGB') && !is.null(subst <- substitute(maxGB))) 
	par.model[['maxGB']] <- CheckArg(maxGB, subst, FALSE)
  if (hasArg('maxmem') && !is.null(subst <- substitute(maxmem))) 
	par.model[['maxmem']] <- CheckArg(maxmem, subst, FALSE)
  if (hasArg('tolIm') && !is.null(subst <- substitute(tolIm))) 
	par.model[['tolIm']] <- CheckArg(tolIm, subst, FALSE)
  if (hasArg('tolRe') && !is.null(subst <- substitute(tolRe))) 
	par.model[['tolRe']] <- CheckArg(tolRe, subst, FALSE)
  if (hasArg('trials') && !is.null(subst <- substitute(trials))) 
	par.model[['trials']] <- CheckArg(trials, subst, FALSE)
  if (hasArg('useprimes') && !is.null(subst <- substitute(useprimes))) 
	par.model[['useprimes']] <- CheckArg(useprimes, subst, FALSE)
  if (hasArg('dependent') && !is.null(subst <- substitute(dependent))) 
	par.model[['dependent']] <- CheckArg(dependent, subst, FALSE)
  if (hasArg('approx_step') && !is.null(subst <- substitute(approx_step))) 
	par.model[['approx_step']] <- CheckArg(approx_step, subst, FALSE)
  if (hasArg('approx_maxgrid') && !is.null(subst <- substitute(approx_maxgrid))) 
	par.model[['approx_maxgrid']] <- CheckArg(approx_maxgrid, subst, FALSE)
  if (hasArg('diameter') && !is.null(subst <- substitute(diameter))) 
	par.model[['diameter']] <- CheckArg(diameter, subst, FALSE)
  if (hasArg('rawR') && !is.null(subst <- substitute(rawR))) 
	par.model[['rawR']] <- CheckArg(rawR, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPintrinsic', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPintrinsic <- new(CLASS_RM, 
	.Data = RPintrinsic,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 13,
	vdim = 1
	)



RPdirect <- function(phi, boxcox) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPdirect', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPdirect <- new(CLASS_RM, 
	.Data = RPdirect,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPhyperplane <- function(phi, boxcox, superpos, maxlines, mar_distr, mar_param, additive) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('superpos') && !is.null(subst <- substitute(superpos))) 
	par.model[['superpos']] <- CheckArg(superpos, subst, FALSE)
  if (hasArg('maxlines') && !is.null(subst <- substitute(maxlines))) 
	par.model[['maxlines']] <- CheckArg(maxlines, subst, FALSE)
  if (hasArg('mar_distr') && !is.null(subst <- substitute(mar_distr))) 
	par.model[['mar_distr']] <- CheckArg(mar_distr, subst, FALSE)
  if (hasArg('mar_param') && !is.null(subst <- substitute(mar_param))) 
	par.model[['mar_param']] <- CheckArg(mar_param, subst, FALSE)
  if (hasArg('additive') && !is.null(subst <- substitute(additive))) 
	par.model[['additive']] <- CheckArg(additive, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPhyperplane', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPhyperplane <- new(CLASS_RM, 
	.Data = RPhyperplane,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 2,
	vdim = 1
	)



RPnugget <- function(phi, boxcox, tol, vdim) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('tol') && !is.null(subst <- substitute(tol))) 
	par.model[['tol']] <- CheckArg(tol, subst, FALSE)
  if (hasArg('vdim') && !is.null(subst <- substitute(vdim))) 
	par.model[['vdim']] <- CheckArg(vdim, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPnugget', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPnugget <- new(CLASS_RM, 
	.Data = RPnugget,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = TRUE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -2
	)



RPsequential <- function(phi, boxcox, back_steps, initial) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('back_steps') && !is.null(subst <- substitute(back_steps))) 
	par.model[['back_steps']] <- CheckArg(back_steps, subst, FALSE)
  if (hasArg('initial') && !is.null(subst <- substitute(initial))) 
	par.model[['initial']] <- CheckArg(initial, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPsequential', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPsequential <- new(CLASS_RM, 
	.Data = RPsequential,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RPspectral <- function(phi, boxcox, sp_lines, sp_grid, prop_factor, sigma) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('sp_lines') && !is.null(subst <- substitute(sp_lines))) 
	par.model[['sp_lines']] <- CheckArg(sp_lines, subst, FALSE)
  if (hasArg('sp_grid') && !is.null(subst <- substitute(sp_grid))) 
	par.model[['sp_grid']] <- CheckArg(sp_grid, subst, FALSE)
  if (hasArg('prop_factor') && !is.null(subst <- substitute(prop_factor))) 
	par.model[['prop_factor']] <- CheckArg(prop_factor, subst, FALSE)
  if (hasArg('sigma') && !is.null(subst <- substitute(sigma))) 
	par.model[['sigma']] <- CheckArg(sigma, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPspectral', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPspectral <- new(CLASS_RM, 
	.Data = RPspectral,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPspecific <- function(phi, boxcox) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPspecific', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPspecific <- new(CLASS_RM, 
	.Data = RPspecific,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPtbm <- function(phi, boxcox, fulldim, reduceddim, layers, lines, linessimufactor, linesimustep, center, points) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('fulldim') && !is.null(subst <- substitute(fulldim))) 
	par.model[['fulldim']] <- CheckArg(fulldim, subst, FALSE)
  if (hasArg('reduceddim') && !is.null(subst <- substitute(reduceddim))) 
	par.model[['reduceddim']] <- CheckArg(reduceddim, subst, FALSE)
  if (hasArg('layers') && !is.null(subst <- substitute(layers))) 
	par.model[['layers']] <- CheckArg(layers, subst, FALSE)
  if (hasArg('lines') && !is.null(subst <- substitute(lines))) 
	par.model[['lines']] <- CheckArg(lines, subst, FALSE)
  if (hasArg('linessimufactor') && !is.null(subst <- substitute(linessimufactor))) 
	par.model[['linessimufactor']] <- CheckArg(linessimufactor, subst, FALSE)
  if (hasArg('linesimustep') && !is.null(subst <- substitute(linesimustep))) 
	par.model[['linesimustep']] <- CheckArg(linesimustep, subst, FALSE)
  if (hasArg('center') && !is.null(subst <- substitute(center))) 
	par.model[['center']] <- CheckArg(center, subst, FALSE)
  if (hasArg('points') && !is.null(subst <- substitute(points))) 
	par.model[['points']] <- CheckArg(points, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPtbm', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPtbm <- new(CLASS_RM, 
	.Data = RPtbm,
	type = c('method for Gauss process', 'method for Gauss process'),
	isotropy = c('non-dimension-reducing', 'non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -3,
	vdim = -1
	)



RPloggaussnormed <- function(variogram, prob, optimize_p, nth, burn.in, rejection) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(variogram)) submodels[['variogram']] <- variogram
  
  if (hasArg('prob') && !is.null(subst <- substitute(prob))) 
	par.model[['prob']] <- CheckArg(prob, subst, FALSE)
  if (hasArg('optimize_p') && !is.null(subst <- substitute(optimize_p))) 
	par.model[['optimize_p']] <- CheckArg(optimize_p, subst, FALSE)
  if (hasArg('nth') && !is.null(subst <- substitute(nth))) 
	par.model[['nth']] <- CheckArg(nth, subst, FALSE)
  if (hasArg('burn.in') && !is.null(subst <- substitute(burn.in))) 
	par.model[['burn.in']] <- CheckArg(burn.in, subst, FALSE)
  if (hasArg('rejection') && !is.null(subst <- substitute(rejection))) 
	par.model[['rejection']] <- CheckArg(rejection, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPloggaussnormed', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPloggaussnormed <- new(CLASS_RM, 
	.Data = RPloggaussnormed,
	type = c('normed process (non-negative values with maximum value being 0 or 1)'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = -3
	)



RPbrorig <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPbrorig', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrorig <- new(CLASS_RM, 
	.Data = RPbrorig,
	type = c('method for Brown-Resnick process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPbrmixed <- function(phi, tcf, xi, mu, s, meshsize, vertnumber, optim_mixed, optim_mixed_tol, lambda, areamat, variobound) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  if (hasArg('meshsize') && !is.null(subst <- substitute(meshsize))) 
	par.model[['meshsize']] <- CheckArg(meshsize, subst, FALSE)
  if (hasArg('vertnumber') && !is.null(subst <- substitute(vertnumber))) 
	par.model[['vertnumber']] <- CheckArg(vertnumber, subst, FALSE)
  if (hasArg('optim_mixed') && !is.null(subst <- substitute(optim_mixed))) 
	par.model[['optim_mixed']] <- CheckArg(optim_mixed, subst, FALSE)
  if (hasArg('optim_mixed_tol') && !is.null(subst <- substitute(optim_mixed_tol))) 
	par.model[['optim_mixed_tol']] <- CheckArg(optim_mixed_tol, subst, FALSE)
  if (hasArg('lambda') && !is.null(subst <- substitute(lambda))) 
	par.model[['lambda']] <- CheckArg(lambda, subst, FALSE)
  if (hasArg('areamat') && !is.null(subst <- substitute(areamat))) 
	par.model[['areamat']] <- CheckArg(areamat, subst, FALSE)
  if (hasArg('variobound') && !is.null(subst <- substitute(variobound))) 
	par.model[['variobound']] <- CheckArg(variobound, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPbrmixed', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrmixed <- new(CLASS_RM, 
	.Data = RPbrmixed,
	type = c('method for Brown-Resnick process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
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
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPbrshifted', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrshifted <- new(CLASS_RM, 
	.Data = RPbrshifted,
	type = c('method for Brown-Resnick process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPbrownresnick <- function(phi, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPbrownresnick', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbrownresnick <- new(CLASS_RM, 
	.Data = RPbrownresnick,
	type = c('method for Brown-Resnick process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPbernoulli <- function(phi, stationary_only, threshold) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('stationary_only') && !is.null(subst <- substitute(stationary_only))) 
	par.model[['stationary_only']] <- CheckArg(stationary_only, subst, FALSE)
  if (hasArg('threshold') && !is.null(subst <- substitute(threshold))) 
	par.model[['threshold']] <- CheckArg(threshold, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPbernoulli', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPbernoulli <- new(CLASS_RM, 
	.Data = RPbernoulli,
	type = c('normed process (non-negative values with maximum value being 0 or 1)'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPgauss <- function(phi, boxcox, stationary_only) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('stationary_only') && !is.null(subst <- substitute(stationary_only))) 
	par.model[['stationary_only']] <- CheckArg(stationary_only, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPgauss', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPgauss <- new(CLASS_RM, 
	.Data = RPgauss,
	type = c('method for Gauss process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPpoisson <- function(phi, intensity) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('intensity') && !is.null(subst <- substitute(intensity))) 
	par.model[['intensity']] <- CheckArg(intensity, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPpoisson', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPpoisson <- new(CLASS_RM, 
	.Data = RPpoisson,
	type = c('Poisson'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
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
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPschlather', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPschlather <- new(CLASS_RM, 
	.Data = RPschlather,
	type = c('Schlather'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RPopitz <- function(phi, xi, mu, s, alpha) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  if (hasArg('alpha') && !is.null(subst <- substitute(alpha))) 
	par.model[['alpha']] <- CheckArg(alpha, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPopitz', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPopitz <- new(CLASS_RM, 
	.Data = RPopitz,
	type = c('Schlather'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



RPsmith <- function(shape, tcf, xi, mu, s) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(shape)) submodels[['shape']] <- shape
  if (hasArg(tcf)) submodels[['tcf']] <- tcf
  
  if (hasArg('xi') && !is.null(subst <- substitute(xi))) 
	par.model[['xi']] <- CheckArg(xi, subst, FALSE)
  if (hasArg('mu') && !is.null(subst <- substitute(mu))) 
	par.model[['mu']] <- CheckArg(mu, subst, FALSE)
  if (hasArg('s') && !is.null(subst <- substitute(s))) 
	par.model[['s']] <- CheckArg(s, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPsmith', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPsmith <- new(CLASS_RM, 
	.Data = RPsmith,
	type = c('Smith'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 4,
	vdim = 1
	)



RPchi2 <- function(phi, boxcox, f) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('f') && !is.null(subst <- substitute(f))) 
	par.model[['f']] <- CheckArg(f, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPchi2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPchi2 <- new(CLASS_RM, 
	.Data = RPchi2,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



RPt <- function(phi, boxcox, nu) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  if (hasArg(phi)) submodels[['phi']] <- phi
  
  if (hasArg('boxcox') && !is.null(subst <- substitute(boxcox))) 
	par.model[['boxcox']] <- CheckArg(boxcox, subst, FALSE)
  if (hasArg('nu') && !is.null(subst <- substitute(nu))) 
	par.model[['nu']] <- CheckArg(nu, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'RPt', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

RPt <- new(CLASS_RM, 
	.Data = RPt,
	type = c('process'),
	isotropy = c('non-dimension-reducing'),
	domain = c('single variable'),
	operator = TRUE,
	monotone = 'mismatch in monotonicity',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = -3
	)



R.minus <- function(x, y, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckMaths(factor, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.minus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.minus <- new(CLASS_RM, 
	.Data = R.minus,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.plus <- function(x, y, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, FALSE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckMaths(factor, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.plus', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.plus <- new(CLASS_RM, 
	.Data = R.plus,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.div <- function(x, y, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, FALSE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckMaths(factor, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.div', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.div <- new(CLASS_RM, 
	.Data = R.div,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.mult <- function(x, y, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, FALSE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckMaths(factor, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.mult', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.mult <- new(CLASS_RM, 
	.Data = R.mult,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.const <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.const', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.const <- new(CLASS_RM, 
	.Data = R.const,
	type = c('mathematical operator', 'trend', 'negative definite', 'tail correlation'),
	isotropy = c('framework dependent', 'framework dependent', 'framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.p <- function(proj, new, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('proj') && !is.null(subst <- substitute(proj))) 
	par.model[['proj']] <- CheckMixed(proj, subst, PROJECTION_NAMES)
  if (!(hasArg('new') && !is.null(subst <- substitute(new)))) new <- UNREDUCED
	par.model[['new']] <- CheckChar(new, subst, ISO_NAMES, FALSE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckMaths(factor, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.p', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.p <- new(CLASS_RM, 
	.Data = R.p,
	type = c('mathematical operator', 'trend'),
	isotropy = c('parameter dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = Inf,
	vdim = 1
	)



R.c <- function(a, b, c, d, e, f, g, h, i, j, l, m, n, o, p, q, ncol, factor) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('a') && !is.null(subst <- substitute(a))) 
	par.model[['a']] <- CheckMaths(a, subst, TRUE)
  if (hasArg('b') && !is.null(subst <- substitute(b))) 
	par.model[['b']] <- CheckMaths(b, subst, TRUE)
  if (hasArg('c') && !is.null(subst <- substitute(c))) 
	par.model[['c']] <- CheckMaths(c, subst, TRUE)
  if (hasArg('d') && !is.null(subst <- substitute(d))) 
	par.model[['d']] <- CheckMaths(d, subst, TRUE)
  if (hasArg('e') && !is.null(subst <- substitute(e))) 
	par.model[['e']] <- CheckMaths(e, subst, TRUE)
  if (hasArg('f') && !is.null(subst <- substitute(f))) 
	par.model[['f']] <- CheckMaths(f, subst, TRUE)
  if (hasArg('g') && !is.null(subst <- substitute(g))) 
	par.model[['g']] <- CheckMaths(g, subst, TRUE)
  if (hasArg('h') && !is.null(subst <- substitute(h))) 
	par.model[['h']] <- CheckMaths(h, subst, TRUE)
  if (hasArg('i') && !is.null(subst <- substitute(i))) 
	par.model[['i']] <- CheckMaths(i, subst, TRUE)
  if (hasArg('j') && !is.null(subst <- substitute(j))) 
	par.model[['j']] <- CheckMaths(j, subst, TRUE)
  if (hasArg('l') && !is.null(subst <- substitute(l))) 
	par.model[['l']] <- CheckMaths(l, subst, TRUE)
  if (hasArg('m') && !is.null(subst <- substitute(m))) 
	par.model[['m']] <- CheckMaths(m, subst, TRUE)
  if (hasArg('n') && !is.null(subst <- substitute(n))) 
	par.model[['n']] <- CheckMaths(n, subst, TRUE)
  if (hasArg('o') && !is.null(subst <- substitute(o))) 
	par.model[['o']] <- CheckMaths(o, subst, TRUE)
  if (hasArg('p') && !is.null(subst <- substitute(p))) 
	par.model[['p']] <- CheckMaths(p, subst, TRUE)
  if (hasArg('q') && !is.null(subst <- substitute(q))) 
	par.model[['q']] <- CheckMaths(q, subst, TRUE)
  if (hasArg('ncol') && !is.null(subst <- substitute(ncol))) 
	par.model[['ncol']] <- CheckMaths(ncol, subst, TRUE)
  if (hasArg('factor') && !is.null(subst <- substitute(factor))) 
	par.model[['factor']] <- CheckMaths(factor, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.c', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.c <- new(CLASS_RM, 
	.Data = R.c,
	type = c('shape function', 'trend'),
	isotropy = c('submodel dependent', 'submodel dependent'),
	domain = c('submodel dependent'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = -1
	)



R.is <- function(x, is, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('is') && !is.null(subst <- substitute(is))) 
	par.model[['is']] <- CheckChar(is, subst, EQ_NAMES, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.is', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.is <- new(CLASS_RM, 
	.Data = R.is,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = 1,
	vdim = 1
	)



R.asin <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.asin', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.asin <- new(CLASS_RM, 
	.Data = R.asin,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.atan <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.atan', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.atan <- new(CLASS_RM, 
	.Data = R.atan,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.atan2 <- function(y, x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.atan2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.atan2 <- new(CLASS_RM, 
	.Data = R.atan2,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.cos <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.cos', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.cos <- new(CLASS_RM, 
	.Data = R.cos,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.sin <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.sin', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.sin <- new(CLASS_RM, 
	.Data = R.sin,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.tan <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.tan', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.tan <- new(CLASS_RM, 
	.Data = R.tan,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.asinh <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.asinh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.asinh <- new(CLASS_RM, 
	.Data = R.asinh,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.atanh <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.atanh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.atanh <- new(CLASS_RM, 
	.Data = R.atanh,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.cosh <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.cosh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.cosh <- new(CLASS_RM, 
	.Data = R.cosh,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.sinh <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.sinh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.sinh <- new(CLASS_RM, 
	.Data = R.sinh,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.tanh <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.tanh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.tanh <- new(CLASS_RM, 
	.Data = R.tanh,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.log <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.log', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.log <- new(CLASS_RM, 
	.Data = R.log,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.expm1 <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.expm1', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.expm1 <- new(CLASS_RM, 
	.Data = R.expm1,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.log1p <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.log1p', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.log1p <- new(CLASS_RM, 
	.Data = R.log1p,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.exp2 <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.exp2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.exp2 <- new(CLASS_RM, 
	.Data = R.exp2,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.log2 <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.log2', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.log2 <- new(CLASS_RM, 
	.Data = R.log2,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.hypot <- function(x, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.hypot', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.hypot <- new(CLASS_RM, 
	.Data = R.hypot,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.cbrt <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.cbrt', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.cbrt <- new(CLASS_RM, 
	.Data = R.cbrt,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.ceil <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.ceil', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.ceil <- new(CLASS_RM, 
	.Data = R.ceil,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.floor <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.floor', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.floor <- new(CLASS_RM, 
	.Data = R.floor,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fmod <- function(x, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.fmod', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fmod <- new(CLASS_RM, 
	.Data = R.fmod,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.round <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.round', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.round <- new(CLASS_RM, 
	.Data = R.round,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.trunc <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.trunc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.trunc <- new(CLASS_RM, 
	.Data = R.trunc,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.erfc <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.erfc', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.erfc <- new(CLASS_RM, 
	.Data = R.erfc,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.lgamma <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.lgamma', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.lgamma <- new(CLASS_RM, 
	.Data = R.lgamma,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.remainder <- function(x, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.remainder', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.remainder <- new(CLASS_RM, 
	.Data = R.remainder,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fdim <- function(x, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.fdim', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fdim <- new(CLASS_RM, 
	.Data = R.fdim,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fmax <- function(x, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.fmax', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fmax <- new(CLASS_RM, 
	.Data = R.fmax,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fmin <- function(x, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, TRUE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, TRUE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.fmin', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fmin <- new(CLASS_RM, 
	.Data = R.fmin,
	type = c('shape function', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.gamma <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.gamma', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.gamma <- new(CLASS_RM, 
	.Data = R.gamma,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.exp <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.exp', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.exp <- new(CLASS_RM, 
	.Data = R.exp,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.erf <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.erf', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.erf <- new(CLASS_RM, 
	.Data = R.erf,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.fabs <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.fabs', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.fabs <- new(CLASS_RM, 
	.Data = R.fabs,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.acos <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.acos', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.acos <- new(CLASS_RM, 
	.Data = R.acos,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.acosh <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.acosh', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.acosh <- new(CLASS_RM, 
	.Data = R.acosh,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.pow <- function(x, y) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  if (hasArg('y') && !is.null(subst <- substitute(y))) 
	par.model[['y']] <- CheckMaths(y, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.pow', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.pow <- new(CLASS_RM, 
	.Data = R.pow,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



R.sqrt <- function(x) {
  cl <- match.call()
  submodels <- par.general <- par.model <- list() 
  
  if (hasArg('x') && !is.null(subst <- substitute(x))) 
	par.model[['x']] <- CheckMaths(x, subst, FALSE)
  
  model <- methods::new('RMmodel', call = cl, name = 'R.sqrt', 
  		submodels = submodels, 
  		par.model = par.model, par.general = par.general)
  return(model)
}

R.sqrt <- new(CLASS_RM, 
	.Data = R.sqrt,
	type = c('mathematical operator', 'trend'),
	isotropy = c('framework dependent', 'framework dependent'),
	domain = c('single variable'),
	operator = FALSE,
	monotone = 'not monotone',
	finiterange = FALSE,
	simpleArguments = TRUE,
	maxdim = -2,
	vdim = 1
	)



