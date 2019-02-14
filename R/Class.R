
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2012 -- 2015 Alexander Malinowski & Martin Schlather
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




setClassUnion("RFspatialDataFrame",
              c("RFspatialGridDataFrame", "RFspatialPointsDataFrame"))
setClassUnion("RFdataFrame", c("RFgridDataFrame", "RFpointsDataFrame"))
setClassUnion("RFsp", c("RFspatialGridDataFrame", "RFspatialPointsDataFrame",
                       "RFgridDataFrame", "RFpointsDataFrame"))

check.validity.n.vdim <- function(object) {
  if (!all(c("n", "vdim") %in% names(object@.RFparams)))
    return("slot '.RFparams' must contain 'n' and 'vdim'")
  var.given <- !is.null(object@.RFparams$has.variance) &&
    object@.RFparams$has.variance
  nc <- (object@.RFparams$n + var.given) * object@.RFparams$vdim

  if (nc != ncol(object@data)) {
    stop("number of data at each location (=", ncol(object@data),
         ") does not match the expected ones (=", nc,
         ").\n  The latter is based on the information or assumption that there are\n  ",
         object@.RFparams$n, " repetition(s) of ",
         object@.RFparams$vdim, " variable(s)",
         if (var.given) " and the variances",
         " given.",
         if (object@.RFparams$n * object@.RFparams$vdim == 1) "\nEither wrong parameters have been given for 'RFspatialGridDataFrame' or 'RFspatialPointsDataFrame' or an explicite transformation of the data from 'sp' objects to 'RFsp' objects has not been performed with 'sp2RF'.")
  }
  return(TRUE)
}




## definition of class CLASS_CLIST
setClass(CLASS_CLIST, 
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

## rules for validity checking of CLASS_CLIST objects
isModel <- function(model) return(is.list(model) && is.character(model[[1]]) )
isRMmodel <- function(x) is(x, class2=CLASS_CLIST)
setValidity(CLASS_CLIST, 
            function(object){
              
              isNUMorDEFAULT <- function(x) {
                # Print("class.R", x, RM_DEFAULT)
                (is.numeric(x) || is(x, "function") || is(x, "call") ||
		  is.environment(x) || is.character(x) ||
                 (length(x) ==1 && x==RM_DEFAULT) || is.logical(x) ||
                 is.matrix(x) || is.list(x))
              }
              isRMmodelorNUMorDEFAULT <- function(x)
                isRMmodel(x) || isNUMorDEFAULT(x)
              
              if (length(object@submodels) > 0)
                if (!all(unlist(lapply(object@submodels, FUN = isRMmodel))))
                  return(paste("submodels must be of class '", CLASS_CLIST, "'",
			       sep=""))

              passed.params <- c(object@par.model, object@par.general)
              if (length(passed.params) > 0) {
                if (!all(unlist(lapply(passed.params,
                                       FUN = isRMmodelorNUMorDEFAULT)) ))
                  return("all parameters must be of class numeric or logical or RMmodel") 
              }
                                    
              if(!is.null(object@par.general$var) &&
                 !isRMmodel(object@par.general$var)) {
                if (length(object@par.general$var) > 1) return("not scalar")
                if (!is.na(object@par.general$var) &&
                    (object@par.general$var!=RM_DEFAULT)){
                  if (object@par.general$var < 0)
                    return("negative variance")
                  }
              }
              
              if(!is.null(object@par.general$scale) &&
                 !isRMmodel(object@par.general$scale)) {
                if (length(object@par.general$scale) > 1) return("not scalar")
                if (!is.na(object@par.general$scale) &&
                    object@par.general$scale!=RM_DEFAULT)
                if(object@par.general$scale <= 0)
                  return("non-positive scale")
              }
              
              if(!is.null(object@par.general$Aniso) &&
                 !isRMmodel(object@par.general$Aniso) &&
                 !all(object@par.general$Aniso==RM_DEFAULT) &&
                 length(object@par.general$Aniso) > 1)
                if(!is.matrix(object@par.general$Aniso))
                  return("'Aniso' must be a matrix")
              
              if(!is.null(object@par.general$proj) &&
              #   !isRMmodel(object@par.general$proj) &&
              #   !is.na(object@par.general$proj) &&
                 !all(object@par.general$proj==RM_DEFAULT)){
                if(!is.vector(object@par.general$proj) ||
                   any(is.na(object@par.general$proj)) ||
                   any(object@par.general$proj == 0) ||
                   any(object@par.general$proj < -2) ||
                   any(object@par.general$proj !=
                       as.integer(object@par.general$proj))
                   )
                  return("proj must be a vector of non-negative integers")
              }
             return(TRUE)
            })


## definition of class CLASS_RM ################################
setClass(CLASS_RM, contains ="function",
         representation(
                        type = "character",
                        domain = "character",
                        isotropy = "character",
                        operator = "logical",
                        monotone = "character",
                        finiterange = "logical",
                        maxdim = "numeric",      # max possible dimension
                        simpleArguments = "logical",
                        vdim = "numeric"         # ??
                        )
         )


## definition of class 'RMmodelFit'
setClass(CLASS_FIT,  contains=CLASS_CLIST,
         representation(
             formel = "ANY",
             likelihood = "numeric",
             variab = "matrix",
             param = "matrix",
             globalvariance  = "ANY",
             covariat = "ANY",
             hessian = "ANY",
             AIC = "numeric",
             AICc = "numeric",
             BIC = "numeric",
             residuals = "ANY"
             )
         )


## definition of class 'RFempVariog'
setClass("RFempVariog", 
         representation(centers = "ANY",
                        empirical = "ANY",
                        var = "ANY",
                        sd = "ANY",
                        n.bin = "ANY",
                        phi.centers = "ANY",
                        theta.centers = "ANY",
                        T = "ANY",
                        vdim = "ANY",
                        coordunits = "character",
                        varunits = "character",
                        call = "ANY",
                        method = "ANY"
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
         representation(Z = "list",
                        ev="RFempVariog",
                        table = "data.frame",
                        n.variab = "integer",
                        n.param = "integer",
                        n.covariates = "integer",
                        lowerbounds = CLASS_CLIST,
                        upperbounds = CLASS_CLIST,
                        transform = "list",
                        #vario = "character",
                        coordunits = "character",
                        varunits = "character",
                        number.of.data = "integer",
                        modelinfo = "data.frame",
                        number.of.parameters = "integer",
                        p.proj = "integer",
                        v.proj = "integer",
                        x.proj = "ANY",  ## integer or logical
                        fixed = "ANY",
                        true.tsdim = "integer",
                        true.vdim = "integer",
                        report = "character",
                        submodels = "ANY",
                        autostart = CLASS_FIT,
                        users.guess = CLASS_FIT, # Martin: 2.4.: eingefuegt
                        self = CLASS_FIT,
                        plain = CLASS_FIT,
                        sqrt.nr = CLASS_FIT,
                        sd.inv = CLASS_FIT,
                        internal1 = CLASS_FIT,
                        internal2 = CLASS_FIT,
                        internal3 = CLASS_FIT,
                        ml = CLASS_FIT
                        #ml.residuals = "ANY" # matrix or RFsp
                        )
         )

## generic S4 method for 'plot'
setGeneric("plot")
