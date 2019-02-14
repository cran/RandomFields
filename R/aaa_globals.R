## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
## Sebastian Gross
##
##
## Copyright (C) 2017 - 2018 Martin Schlather
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



COORD_NAMES_CART <- c("x", "y", "z", "T")
COORD_NAMES_GENERAL <- c("coords.x", "coords.T") ## check general_coordinates if changed
COORD_NAMES_EARTH <- c("longitude", "latitude", "height", "time")

SYMBOL_L_PAR <- "("
SYMBOL_R_PAR <- ")"
SYMBOL_PLUS <- '+'
SYMBOL_MULT <- '*'
DOLLAR <- c("$", "RMS")

RM_PLUS <- c("RMplus", SYMBOL_PLUS)
RM_MULT <- c("RMmult", SYMBOL_MULT)
#RM_MIXED <- c( "RMmixed", "mixed") 
##RM_INTERNALMIXED <- "internalRMmixed"
RM_TREND <- c("RMtrend", "trend")
RM_DISTR <- c('RRdistr', 'Distr')
RM_USER <- c('RMuser', 'U')
RM_NUGGET <- c("RMnugget", "nugget")
RM_COVARIATE <- "RMcovariate"
RM_DEFAULT <- "RFdefault"
RM_DECLARE <- "RMdeclare"
R_C <- "R.c"
R_CONST <- "R.const"
# R_CC <- c(R_C, R_CONST)

CLASS_CLIST <- 'RMmodel'
CLASS_RM <- 'RMmodelgenerator'
CLASS_CLISTFIT <- "RMmodelFit"
CLASS_FIT <- 'RMmodelFit'


ERRMIXED <- "Please use 'RMmixed' and not '@' from version 3.1.54 on"
syntaxError <- "Malformed model expression -- maybe you have used a wrong or obsolete definition, or just used an incorrect option name or incorrect ordering of the arguments. See ?RMmodel for the model definition. Check manual for further information (RMmodel, RFsimulate)"

 
isPosDef <- function(type) {
  if (is.character(type)) type <- pmatch(type, TYPE_NAMES, duplicates.ok=TRUE)-1
  ##  .C(C_isPosDef, as.integer(type))$type
  type==TcfType | type == PosDefType | type == ManifoldType
}
isVariogram <- function(type) { 
  if (is.character(type)) type <- pmatch(type, TYPE_NAMES, duplicates.ok=TRUE)-1
  ##  .C(C_isNefDef, as.integer(type))$type
  isPosDef(type) | type == VariogramType
}

LSQMETHODS <- c("self", "plain", "sqrt.nr", "sd.inv", "internal") 
MLMETHODS <- c("ml") # "reml", "rml1"),

par.storage <- ".RandomFields.par"
.RandomFields.env <- new.env()


## ACHTUNG! Die in C definierten PL_* haben andere Bedeutung
PL_IMPORTANT 	<- as.integer(1)
PL_SUBIMPORTANT 	<- as.integer(2)
PL_DETAILSUSER <- as.integer(3)## currently unused
PL_RECURSIVE 	<- as.integer(4)
PL_STRUCTURE 	<- as.integer(5)
PL_ERRORS 	<- as.integer(6)

PL_FCTN_DETAILS 	<- as.integer(7)
PL_FCTN_SUBDETAILS 	<- as.integer(8)

PL_COV_STRUCTURE 	<- as.integer(7)
PL_DIRECT_SEQU 	<- as.integer(8)
PL_DETAILS 	<- as.integer(9)
PL_SUBDETAILS 	<- as.integer(10)


