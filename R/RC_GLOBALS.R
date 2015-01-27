


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

isPosDef <- function(type) {
  if (is.character(type)) type <- pmatch(type, TYPENAMES, duplicates.ok=TRUE)-1
  ##  .C("isPosDef", as.integer(type))$type
  type==TcfType | type == PosDefType | type==UndefinedType
}
isVariogram <- function(type) { 
  if (is.character(type)) type <- pmatch(type, TYPENAMES, duplicates.ok=TRUE)-1
  ##  .C("isNefDef", as.integer(type))$type
  isPosDef(type) | type == VariogramType
}




###############################################################################
##                              OTHER DEFINITIONS                            ##
###############################################################################

LSQMETHODS <- c("self", "plain", "sqrt.nr", "sd.inv", "internal") 
MLMETHODS <- c("ml") # "reml", "rml1"),
