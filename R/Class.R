##########################################################################
## classes for 2- and higher-dimensional data objects, based on ##########
## Spatial classes from 'sp'-package                            ##########


setClass("RFspatialGridDataFrame", contains ="SpatialGridDataFrame",
         representation(.RFparams="list") )
setValidity("RFspatialGridDataFrame", 
            function(object) {
              return(check.validity.n.vdim(object))
            })


setClass("RFspatialPointsDataFrame", contains ="SpatialPointsDataFrame",
         representation(.RFparams="list") )
setValidity("RFspatialPointsDataFrame", 
            function(object) {
              return(check.validity.n.vdim(object))
            })


## classes for 1-dimensional data objects ################################

setClass("RFgridDataFrame", 
         representation(data="data.frame", grid="GridTopology",
                        .RFparams="list"))
setValidity("RFgridDataFrame", 
            function(object) {
              if (length(object@grid@cells.dim) > 1)
                return("this class is only for 1-dimensional coordinates")
              if (nrow(object@data) != object@grid@cells.dim)
                return("data must have the same length as the grid length")
              return(check.validity.n.vdim(object))
            })

setClass("RFpointsDataFrame", 
         representation(data="data.frame", coords="matrix",
                        .RFparams="list"),
         prototype(data=data.frame(NULL), coords=NULL, .RFparams=list()))
setValidity("RFpointsDataFrame", 
            function(object) {
              if (nrow(object@data) != length(object@coords))
                return("data and coords must have the same length")
              if (ncol(object@coords)!=1)
                return("coords must have exactly 1 column, otherwise use class 'RFspatialPointsDataFrame'")
              return(check.validity.n.vdim(object))
            })



setClassUnion("RFsp", c("RFspatialGridDataFrame", "RFspatialPointsDataFrame",
                        "RFgridDataFrame", "RFpointsDataFrame"))




check.validity.n.vdim <- function(object) {
  if (!all(c("n", "vdim") %in% names(object@.RFparams)))
    return("slot '.RFparams' must contain 'n' and 'vdim'")
  stopifnot((object@.RFparams$n +
             (!is.null(object@.RFparams$has.variance) && object@.RFparams$has.variance)
             ) *
            object@.RFparams$vdim == ncol(object@data))
  return(TRUE)
}




## definition of class ZF_MODEL
setClass('RMmodel', 
         representation(
                        # call='RMexp(var=1, sclae=1, Aniso=id, proj=id)
                        call = "language",
                        # name='RMexp'
                        name = "character",
                        # submodels=NULL, submodels=list(RMmodel1, RMmodel2)
                        submodels = "list",
                        # model specific parameter 
                        par.model = "list",
                        # var=1, scale=1, Aniso=id, proj=id 
                        par.general = "list"
                        )
         )

## rules for validity checking of ZF_MODEL objects
setValidity('RMmodel', 
            function(object){
              #return(TRUE)
              isRMmodel <- function(x) is(x, class='RMmodel')
              
              isNUMorDEFAULT <- function(x) {
                # Print("class.R", x, ZF_DEFAULT_STRING)
                (is.numeric(x) ||
                 class(x) == "function" || class(x) == "call" ||
                 is.environment(x) ||
                 x==ZF_DEFAULT_STRING || is.logical(x) ||
                 is.matrix(x) || is.list(x))
              }
              isRMmodelorNUMorDEFAULT <- function(x)
                isRMmodel(x) || isNUMorDEFAULT(x)
              
              if (length(object@submodels) > 0)
                if (!all(unlist(lapply(object@submodels, FUN = isRMmodel))))
                  return("submodels must be of class 'RMmodel'")

              passed.params <- c(object@par.model, object@par.general)
              
              if (length(passed.params) > 0) {
                #if (any(unlist(lapply(passed.params, FUN = isRMmodel))))
                #  return("parameters must NOT be of class 'RMmodel'; probably the model has less submodels than you have passed or you might have mixed up e.g., RPgauss with RMgauss")
                if (!all(unlist(lapply(passed.params,
                                       FUN = isRMmodelorNUMorDEFAULT)) ))
                  return("all parameters must be of class numeric or logical or RMmodel") 
              }
                                    
              if(!is.null(object@par.general$var) &&
                 !isRMmodel(object@par.general$var) &&
                 !is.na(object@par.general$var) &&
                 !(object@par.general$var==ZF_DEFAULT_STRING))
                if (object@par.general$var < 0)
                  return("negative variance")
              
              if(!is.null(object@par.general$scale) &&
                 !isRMmodel(object@par.general$scale) &&
                 !is.na(object@par.general$scale) &&
                 !(object@par.general$scale==ZF_DEFAULT_STRING))
                if(object@par.general$scale < 0)
                  return("negative scale")
              
              if(!is.null(object@par.general$Aniso) &&
                 !isRMmodel(object@par.general$Aniso) &&
                 !is.na(object@par.general$Aniso) &&
                 !all(object@par.general$Aniso==ZF_DEFAULT_STRING))
                if(!is.matrix(object@par.general$Aniso))
                  return("'Aniso' must be a matrix")
              
              if(!is.null(object@par.general$proj) &&
                 !isRMmodel(object@par.general$proj) &&
                 !is.na(object@par.general$proj) &&
                 !all(object@par.general$proj==ZF_DEFAULT_STRING)){
                if(!is.vector(object@par.general$proj) ||
                   any(object@par.general$proj <= 0) ||
                   any(object@par.general$proj !=
                       as.integer(object@par.general$proj))
                   )
                  return("proj must be a vector of non-negative integers")
              }
              
              return(TRUE)
            })


## definition of class 'RMmodelgenerator' ################################
setClass('RMmodelgenerator', contains ="function",
         representation(
                        type = "character",
                        domain = "character",
                        isotropy = "character",
                        operator = "logical",
                        normalmix = "logical", # is normal mixture?
                        finiterange = "logical",
                        maxdim = "numeric",      # max possible dimension
                        vdim = "numeric"         # ??
                        )
         )


## definition of class 'RMmodelExt'
setClass('RMmodelExt',  contains='RMmodel',
         representation(likelihood = "numeric",
                        trend = "numeric",
                        residuals = "ANY"
                        )
         )


## definition of class 'RFempVariog'
setClass("RFempVariog", 
         representation(centers = "ANY",
                        emp.vario = "ANY",
                        var = "ANY",
                        sd = "ANY",
                        n.bin = "ANY",
                        phi.centers = "ANY",
                        theta.centers = "ANY",
                        T = "ANY",
                        call = "ANY"
                        )
         )

## rules for validity checking of 'RFempVariog' objects
setValidity("RFempVariog", 
            function(object){
              if(!(is.null(object@call)) && !(is(object@call, "language")))
                return("slot 'call' must be NULL or of class 'language'")
              return(TRUE)
            })




## definition of class 'RFfit'
setClass("RFfit", 
         representation(ev="list",
                        table = "matrix",
                        lowerbounds ='RMmodel',
                        upperbounds ='RMmodel',
                        transform = "list",
                        #vario = "character",
                        autostart = 'RMmodelExt',
                        users.guess = 'RMmodelExt', # Martin: 2.4.: eingefuegt
                        self = 'RMmodelExt',
                        plain = 'RMmodelExt',
                        sqrt.nr = 'RMmodelExt',
                        sd.inv = 'RMmodelExt',
                        internal1 = 'RMmodelExt',
                        internal2 = 'RMmodelExt',
                        internal3 = 'RMmodelExt',
                        ml = 'RMmodelExt'#,
                        #ml.residuals = "ANY" # matrix or RFsp
                        )
         )






## generic S4 method for 'plot'

setGeneric("plot")
