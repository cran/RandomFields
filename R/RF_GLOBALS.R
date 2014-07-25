


## @FILE-STARP******************************************************************
# @NAME		ZF_GLOBALS
# @DESCRIPTION	Any value that is used throughout the randomfield package
#               has to appear here
# @AUTHOR	Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
#               Martin Schlather
# @DATE		26.08.2011 (Gross), 2012 -- 2013 (Schlather)
#
# @FILE-END*********************************************************************




###############################################################################
##                        DEFINITIONS OF SYMBOLS                             ##
###############################################################################

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_PLUS/MAL/SYMBOLS
# @DESCRIPTION	The + operator in any valid model formula
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_PLUS <- '+'
ZF_PLUS <- c("RMplus", ZF_SYMBOLS_PLUS)

ZF_SELECT <- c("RMselect", "select")
ZF_PLUSSELECT <- c(ZF_PLUS, ZF_SELECT)

ZF_SYMBOLS_MULT <- '*'
ZF_MULT <- c("RMmult", ZF_SYMBOLS_MULT)

## Special Models
DOLLAR <- c("$", "RMS")
ZF_DOLLAR <- rev(DOLLAR)



# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_AT
# @DESCRIPTION	The former @ operator in any valid model formula used to create fixed effects must ffs not be "*"
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_AT <- "@"

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_L_PAR
# @DESCRIPTION	Left parenthesis
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_L_PAR <- "("

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_SYMBOLS_R_PAR
# @DESCRIPTION	Right parenthesis
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		26.08.2011
# @GLOBAL-END*******************************************************************
ZF_SYMBOLS_R_PAR <- ")"



###############################################################################
##                        DEFINITIONS OF MODELNAMES                          ##
###############################################################################

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_FIXED and other names of special models
# @DESCRIPTION	The function name of fixed effects
# @AUTHOR       Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
#               Martin Schlather
# @DATE		29.08.2011 (Gross) and 2013--2013 (Schlather)
# @GLOBAL-END*******************************************************************
ZF_FIXED <- "RMfixed"
ZF_INTERNALMIXED <- "internalRMmixed"
ZF_TREND <- c("RMtrend", "trend")
ZF_TRENDFCT <- paste(ZF_TREND[1], "(", sep="")
ZF_DISTR <- c('RRdistr', 'Distr')
ZF_USER <- c('RMuser', 'U')
ZF_COORD <- "RMcoord"
ZF_MODEL <- "RMmodel"

ZF_MIXED <- c( "RMmixed", "mixed") 
ZF_NUGGET <- c("RMnugget", "nugget")
ZF_MODELEXT <- "RMmodelFit"


# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_MODEL_FACTORY
# @DESCRIPTION	Each covariance model is an object of this class
# @AUTHOR		Sebastian Gross <sebastian.gross@stud.uni-goettingen.de>
# @DATE		29.08.2011
# @GLOBAL-END*******************************************************************
ZF_MODEL_FACTORY <- "RMmodelgenerator"

# @GLOBAL-STARP*****************************************************************
# @NAME		ZF_DEFAULT_STRING
# @DESCRIPTION	var, scale, Aniso, proj  get this value assigned by functions
#               of class 'RMmodelgenerator', if no such arguments are passed
#               by the user. NULL is passed to C-level
# @AUTHOR       A Malinowski <malinows@math.uni-goettingen.de>
# @DATE		29.08.2011
# @GLOBAL-END*******************************************************************
ZF_DEFAULT_STRING <- "RFdefault"
ZF_MODEL_PREFIX <- "RM"



###############################################################################
##                        STATIONARITY AND ISOTROPY                          ##
###############################################################################



#ZF_NULL <- '<null>'
#RC_TYPE_DOMAINS <-
#  cbind(c('domain', 'variogram', 'process', 'gauss method',
#          'point-shape function', 'distribution', 'shape function',
#          'trend', 'interface', 'undefined', 'of other type'),
#        c('covariance', 'gen. variogram', ZF_NULL, ZF_NULL,
#          ZF_NULL, ZF_NULL, ZF_NULL,
#          ZF_NULL, 'kernel interface', 'undefined kernel', 'of other kernel'),
#        c('param dep., pos. definite', 'param dep., neg. definite',
#          ZF_NULL, ZF_NULL,
#          ZF_NULL, ZF_NULL, ZF_NULL,
#          ZF_NULL, 'param. dep. interface', 'undefined, parametric function',
#          'of other param. dep. function')       
#        )        
#RC_TYPE_DOM <- RC_TYPE_DOMAINS[RC_TYPE_DOMAINS != ZF_NULL]


TRANS_INV <- as.integer(0)
KERNEL <- as.integer(1)
PREVMODELD <- as.integer(2)
DOMAIN_MISMATCH <- as.integer(3)
RC_DOMAIN <- c('single variable', 'kernel', 'framework dependent', 'mismatch')
     

RC_ISOTROPIC <- as.integer(0)
RC_SPACEISOTROPIC <- as.integer(1)
ZERORC_SPACEISOTROPIC <- as.integer(2)
VECTORRC_ISOTROPIC <- as.integer(3)
SYMMETRIC <- as.integer(4)
RC_CARTESIAN_COORD <- as.integer(5)
EARTH_COORD <- as.integer(6)
SPHERICAL_COORD <- as.integer(7)
CYLINDER_COORD <- as.integer(8)
UNREDUCED <- as.integer(9)
PREVMODELI <- as.integer(10)
RC_ISOTROPY <- c("isotropic", "space-isotropic", "zero-space-isotropic",
                 "vector-isotropic", "symmetric", "cartesian system",
                 "earth system", "spherical system", "cylinder system",
                 "non-dimension-reducing", "parameter dependent", "<mismatch>")
MON_MISMATCH <- 5
MON_PARAMETER <- MON_MISMATCH - 1
NOTMONOTONE <- MON_MISMATCH + 0
BERNSTEIN <-  MON_MISMATCH + 5
RC_MONOTONE <- c("mismatch in monotonicity",
                 "submodel dependent monotonicity",
                "previous model dependent monotonicity",
                "parameter dependent monotonicity",
                "not monotone", "monotone", "Gneiting-Schaback class", 
                "normal mixture", "completely monotone", "Bernstein")

## Coding of stationarity and isotropy in RF.h, see also RC_DOMAIN
## and RC_ISOTROPY


# mostly unused
TcfType <- as.integer(0)
PosDefType <- as.integer(1)
NegDefType <- as.integer(2) ## the only one currently needed
ProcessType <- as.integer(3)
MethodType <- as.integer(4)  #/* Gauss Only */
BrType <- as.integer(5)
PointShapeType <- as.integer(6)
.RandomType <- as.integer(7) ## NICHT GROSS SCHREIBEN, da von NAMESPACE ERFASSR
ShapeType <- as.integer(8)
TrendType <- as.integer(9)
InterfaceType <- as.integer(10) 
UndefinedType <- as.integer(11) 
OtherType <- as.integer(12) 
RC_TYPE <- c("tail correlation function", "positive definite",
             "negative definite", "process", 
    "method for Gauss processes", "method for Brown-Resnick processes",
    "shifted shape function",
    "distribution family", "shape function", "trend", "interface",
             "undefined", "other type")
isPosDef <- function(type) {
  (is.numeric(type) && (type==TcfType || type == PosDefType ||
                        type==UndefinedType)) ||
  (is.character(type) &&
   (type == RC_TYPE[TcfType+1] || type == RC_TYPE[PosDefType + 1]
    || type == RC_TYPE[UndefinedType + 1]))
}
isNegDef <- function(type) {
  isPosDef(type) ||
  (is.numeric(type) && type == NegDefType) ||
  (is.character(type) && type == RC_TYPE[NegDefType+1])
}

COORD_SYSTEMS <- c("auto", "cartesian", "earth")
ZF_CARTESIAN_COORD_NAMES <- c("auto", "cartesian")


#RC_TYPE_PREFIX <- .Call("GetCathegoryNames")
                            #c("RM", "RP", "", "RL", "RM", "RM", "RM", "RM")


###############################################################################
##                  TYPES AND OTHER CHARACTERISING FLAGS                     ##
###############################################################################

## SetAndGetModelInfos:
## Flags that characterise parameters, see RF.h
## this is only important in MLE
VARPARAM <- 0  ## to be consistent with the C definitions in RF.h
SIGNEDVARPARAM <- 1
SDPARAM <- 2
SIGNEDSDPARAM <- 3
SCALEPARAM <- 4
DIAGPARAM <- 5
ANISOPARAM <- 6
INTEGERPARAM <- 7 ## NEU
ANYPARAM <- 8
TRENDPARAM <- 9
NUGGETVAR <- 10
MIXEDVAR <- 11
.REGRESSION <- 12
CRITICALPARAM <- 13
ZF_TYPEOF_PARAM <- c("var", "signed var", "sd", "signed sd", "scale", "diag",
                     "aniso", "integer", "unspecfd", "trend",  "nugget",
                     "mixed var", "regress", "any")


## Flags that characterise a component within a mixed model definition
DetTrendEffect <- 0  ## trend, nichts wird geschaetzt
DeterministicEffect <- 1 ## nichts wird geschaetzt
FixedTrendEffect <- 2
FixedEffect <-  3 ## trend is also converted to FixedEffect ?
.RandomEffect <- 4 ## b is random, no variance is estimated; cov. matrix
.RVarEffect <-   5 ## b is random, variance is estimated; covariance matrix
LargeEffect <- 6 ## wie RandomEffect, aber gross, somit keine Optimierung
LVarEffect <- 7 ## wie RVarEffekt, aber gross, somit keine Optimierung
SpaceEffect <-  8 ## spatial covariance model for random effect
SpVarEffect <-  9 ## spatial covariance model for random effect
Primitive <-    10 ## (but not simple) primitive and remaining might be worth to distinguish
Simple <-      11 ## if C->primitive and domain and isotropic  and vdim=1
.RemainingError <- 12




###############################################################################
##                              OTHER DEFINITIONS                            ##
###############################################################################

## Levels of printing
# cross check with PL_C_* in RF.h !!
PL.IMPORPANT <- 1
PL.SUBIMPORPANT <- 2
PL.RECURSIVE.CALL <- 3
PL.RECURSIVE.DETAILS <- 4
PL.FCTN.STRUCTURE <- 5
PL.FCTN.ERRORS <- 6 ## only those that are caught internally
PL.FCTN.DETAILS <- 7
PL.FCTN.SUBDETAILS <- 8


## predefined registers
## 0 -- (MAXFIELDS-1) : register that be used freely
MAXFIELDS <- as.integer(10)
MODEL.USER <- as.integer(MAXFIELDS + 0)
MODEL.UNUSED <- as.integer(MAXFIELDS + 1)
MODEL.INTERN <- as.integer(MAXFIELDS + 2)
MODEL.SPLIT <- as.integer(MAXFIELDS + 3)
MODEL.GUI <- as.integer(MAXFIELDS + 4)
MODEL.MLE <- as.integer(MAXFIELDS + 5)
MODEL.MLESPLIT <- as.integer(MAXFIELDS + 6)
MODEL.MLETREND <- as.integer(MAXFIELDS + 7)
MODEL.BOUNDS <- as.integer(MAXFIELDS + 8)
MODEL.KRIGE <- as.integer(MAXFIELDS + 9)
MODEL.COND <- as.integer(MAXFIELDS + 10)
MODEL.ERR <- as.integer(MAXFIELDS + 11)
MODEL.MAX <- as.integer(MODEL.BOUNDS + 1)

MaxNameCharacter <- as.integer(200)

GETMODEL_AS_SAVED <- as.integer(0)
GETMODEL_DEL_NATSC <- as.integer(1)
GETMODEL_SOLVE_NATSC <- as.integer(2)
GETMODEL_DEL_MLE <- as.integer(3)
GETMODEL_SOLVE_MLE <- as.integer(4)


LSQMETHODS <- c("self", "plain", "sqrt.nr", "sd.inv", "internal") 
MLMETHODS <- c("ml") # "reml", "rml1"),

DUPFALSE <- FALSE
