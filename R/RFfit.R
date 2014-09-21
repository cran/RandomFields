### !!!!!!!!!!!! ACHTUNG !!!!!!!!!!!! TREND als cov-fct muss
### noch programmiert werden !!!

# RFsimulate:  Not implemented yet: If \code{model} is a formula or of class
#    \command{\dQuote{\link{RFformula}}},
#    the corresponding linear mixed model of the type
 #   \deqn{response = W*b + Z*u + e} is simulated

##   source("~/R/RF/RandomFields/R/MLES.R")

## PrintLevels
## 0 : no message
## 1 : important error messages
## 2 : warnings
## 3 : minium debugging information
## 5 : extended debugging information

## jetzt nur noch global naturalscaling (ja / nein)
## spaeter eine Funktion schreibbar, die den naturscaling umwandelt;
##   im prinzipt CMbuild, aber ruechwaers mit 1/newscale und eingefuegt
##   in eventuell schon vorhandene $ operatoren


#Beim paper lesen im Zug nach Muenchen heute morgen ist mir eine Referenz zu einem R Paket "mlegp: Maximum likelihood estimates of Gaussian processes" aufgefallen. Ist Dir aber sicher schon bekannt! 

#  stop("")
  # problem: natscale; im moment 2x implementiert, 1x mal ueber
  # scale/aniso (user) und einmal gedoppelt -- irgendwas muss raus

## LSQ variogram fuer trend = const.
## kann verbessert werden, insb. fuer fixed effects, aber auch eingeschraenkt
## fuer random effects -> BA/MA


## REML fehlt

## users.guess muss in eine List von meheren Vorschlaegen umgewandelt werden !!! Und dann muss RFfit recursiver call mit allen bisherigen Werden laufen !!


## NAs in data mit mixed model grundsaetzlich nicht kombinierbar !
## NAs in data mit trend (derzeit) nicht kombinierbar

## zentrale C -Schnittstellen
##    .C("PutValuesAtNA", RegNr, param, PACKAGE="RandomFields")

## bins bei Distances automatisch


## bei repet sind die Trends/fixed effects gleich, es muessen aber die
## random effects unterschiedlich sein.
## bei list(data) werden auch trend/fixed effects unterschiedlich geschaetzt.


## Erweiterungen: Emilio's Bi-MLE, Covarianz-Matrix-INversion per fft oder
## per INLA, grosse Datensaetze spalten in kleinere "unabhaengige".


###################################
## !!! Mixed Model Equations !!! ##
###################################


RFfit <-
  function(model, x, y=NULL, z=NULL, T=NULL,  grid, data, 
           lower=NULL, upper=NULL, 
           bc_lambda, ## if missing then no BoxCox-Trafo
           methods, # "reml", "rml1"),
           sub.methods,
           ## "internal" : name should not be changed; should always be last
           ##              method!
           optim.control=NULL,
           users.guess=NULL,  
           distances=NULL, dim,
           transform=NULL,
           ##type = c("Gauss", "BrownResnick", "Smith", "Schlather",
           ##             "Poisson"),
           ... )
{

  #Print(RFoptions()$fit); xxxxxx###
  .C("NoCurrentRegister")

  RFoptOld <- internal.rfoptions(xyz=!is.null(y),...,
                                 RELAX=isFormulaModel(model))
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]

  if (RFopt$general$vdim_close_together)
    stop("'vdim_close_together' must be FALSE")

 # Print(RFoptions()$fit, RFoptOld, RFopt$fit); xxxxxx###
  
  model <- PrepareModel2(model, ...) 
  if (!is.null(lower)) lower <- PrepareModel2(lower, ...)
  if (!is.null(upper)) upper <- PrepareModel2(upper, ...)
  if (!is.null(users.guess)) users.guess <- PrepareModel2(users.guess, ...)
  if (!is.null(optim.control$parscale) && !is.numeric(optim.control$parscale))
    optim.control$parscale <- PrepareModel2(optim.control$parscale, ...)

    
  if (model[[1]] %in% c("RPpoisson", "poisson")) {
    fit.poisson()
  } else if (model[[1]] %in% c("BRmixed", "BRshifted", "BRmixedIntern",
                               "RFbrownresnick")) {
    fit.br()
  } else if (model[[1]] %in% c("RPschlather", "extremalgauss")) {
    fit.extremal.gauss()
  } else if (model[[1]] %in% c("RPsmith", "smith")) {
    fit.smith()
  } else if (model[[1]] %in% c("RPbernoulli", "binaryprocess")) {
    fit.bernoulli()    
  } else {
    do.call("rffit.gauss",
            c(list(model=model, y=y, z=z, T=T, data=data,
                   lower=lower, upper=upper, 
                   users.guess=users.guess,  
                   distances=distances,
                   optim.control=optim.control,
                   transform=transform,
                   recall = FALSE),
              if (!missing(x)) list(x=x),
              if (!missing(grid)) list(grid = grid),
              if (!missing(bc_lambda)) list(bc_lambda=bc_lambda),
              if (!missing(methods))  list(mle.methods = methods),
              if (!missing(sub.methods)) list(lsq.methods=sub.methods),
              ## "internal" : name should not be changed; should always
                ## be last method!
              if (!missing(dim)) list(dimensions=dim)
              ))
  }
}
