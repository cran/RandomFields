
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
### Diese Datei wandelt RMmodel in eine Liste um
##
## Copyright (C) 2015 -- 2016 Alexander Malinowski, Martin Schlather
##               2017 -- 2017 Martin Schlather
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



reps <- function(n, sign=",") paste(rep(sign, n), collapse="")


## FUNCTION-STARP***********************************************************************************
## NAME		parseModel, isGenuineCovModel, extractSummands
## PARAM	$model - list, formula		
## AUTHOR		Sebastian Gross 
### DATE		29.08.2011
## changed by Martin Schlather 2015 -- 2017
## FUNCTION-END*************************************************************************************

parseModel <- function(model, ..., x=NULL) {
  
  ## check whether $model is already in list syntax
  if (is.list(model)) return(model)

  ## check whether $model has RMmodel syntax
  if (isRMmodel(model)) return(buildCovList(model, x=x))
  
  ## check whether $model has correct formula syntax
  if (!is(model, "formula")) stop(syntaxError)
  
  
  ## extract tokens/summands
  tmpList <- list()
  ## ignore rest of the formula
  rightSide <- tail(as.character(model), 1)
  chars <- strsplit(rightSide, "")[[1]]
 	
  ## toggles parenthesis, eg. whether we confront a toplevel plus or not
  parToggle <- 0
  
  token <- ""
  for (char in chars) {
    if (char == SYMBOL_PLUS && parToggle == 0) {
      tmpList <- c(tmpList, token)
      token <- ""			
    } else {
      if (char != " ") token <- paste(token, char, sep="")		
      if (char == SYMBOL_L_PAR) parToggle <- parToggle+ 1
      if (char == SYMBOL_R_PAR) parToggle <- parToggle- 1
    }
  }
  tmpList <- c(tmpList, token)

  summands <- vector("list", length(tmpList))
  for (i in 1:length(tmpList)) summands[[i]] <- removeParenthesis(tmpList[[i]])

  ## if (length(summands)==1) return(summands)  ## sonst steht unten paste(NULL), was "" gibt
  
  paramEnv <- new.env(parent=.GlobalEnv)
  dots <- list(...)
  assign("dots", dots, envir=paramEnv)
  if (length(dots)>0) 
    eval(parse(text = paste(names(dots), "<- dots[[", 1:length(dots), "]]",
                            collapse=";")),
         envir=paramEnv)
  trendfct <- paste(RM_TREND[1], "(", sep="")
  isGenuineCovModel <- lapply(summands, FUN=function(name)
    ## Martin: habe hier RFtrend hinzugefuegt!
    substr(name, 1, length(trendfct)) != trendfct &&
    (is(try(eval(parse(text=name),envir=paramEnv), silent=TRUE), CLASS_CLIST) ||
     (regexpr("^[[:alnum:]_]+\\([[:print:]]*\\)$", name) == 1 && 
      (exists(fun <- strsplit(name, "\\(")[[1]][1]) && is(get(fun), CLASS_RM))
     )))
  isGenuineCovModel <- unlist(isGenuineCovModel)
                                    
  summands <- c(summands[!isGenuineCovModel],
                list(paste(unlist(summands[isGenuineCovModel]),
                           collapse=SYMBOL_PLUS)))

  if (length(summands) == 1) {
    listModel <- buildFactorList(summands[[1]], Env=paramEnv, ..x..=x)
  } else {
    listModel <- list(SYMBOL_PLUS)		
    ##last <- getLastCovIndex(summands)
    for (i in 1:length(summands))
      listModel[[i+1]] <- buildFactorList(summands[[i]], Env=paramEnv, ..x..=x)
    firstgenuine <- sum(!isGenuineCovModel) + 2
    if (firstgenuine <= length(listModel) && listModel[[firstgenuine]][[1]] == SYMBOL_PLUS) {
      listModel <- c(listModel[1:(firstgenuine - 1)], listModel[[firstgenuine]][-1])
      names(listModel) <- NULL
    }
  }
  
  return(listModel)
}



## FUNCTION-STARP***********************************************************************************
## NAME		buildCovList
## PARAM		$model - RMmodel
## RETURN		list
## REQUIRE	none
## ENSURE		isListModel(output) == TRUE
## SEE		RMmodel, devel-doc
## AUTHOR		Sebastian Gross 
#                        modified 2015--2017 by Martin Schlather
### DATE		26.08.2011
## FUNCTION-END*************************************************************************************
buildCovList <- function(model, x=NULL) {
  if (is.atomic(model) || is.list(model) ||
      is.language(model) || is.environment(model))
    return(model) ## for recursive calling
  
  if (!isRMmodel(model)) stop('model must be of class ', CLASS_CLIST) 
  
  li <- c(list(model@name),
          lapply(model@par.model[model@par.model != RM_DEFAULT],
                 FUN=buildCovList, x=x),
          lapply(model@submodels, FUN=buildCovList, x=x) )
  
  if (li[[1]] == RM_PLUS[1]) li[[1]] <- SYMBOL_PLUS
  if (li[[1]] == RM_MULT[1]) li[[1]] <- SYMBOL_MULT


 
  ##  par.general.is.default <-
  ##    unlist(lapply(model@par.general, FUN=function(x) x==RM_DEFAULT))
  if (length(model@par.general)>0 &
      !all(model@par.general==RM_DEFAULT)) {
    li <- c(DOLLAR[1],
            lapply(model@par.general[model@par.general != RM_DEFAULT],
                   FUN=buildCovList, x=x),
            list(li))
 #   if (length(pos <- which(names(li)=="Aniso")) > 0)
  #    ## in c-level, parameter is called 'A'
  #    names(li)[pos] <- "A"
  }
  
  return(li)           
}




## FUNCTION-STARP***********************************************************************************
## NAME		removeParenthesis
## PARAM		$string - string
## RETURN		string
## REQUIRE	none
## ENSURE		The returned string is one not enclosed by parenthesis
## AUTHOR		Sebastian Gross 
## DATE		26.08.2011
## FUNCTION-END*************************************************************************************
removeParenthesis <- function(string)
{
	# split the string
	chars <- strsplit(string, "")[[1]]

	while (head(chars, 1) == SYMBOL_L_PAR &&
               tail(chars, 1) == SYMBOL_R_PAR)
	{
		chars <- chars[-1]
		chars <- chars[-length(chars)]
	}

	# rejoin the string
	string <- ""
	for (char in chars)
	{
		string <- paste(string, char, sep="")
	}
	return (string)
}

## FUNCTION-STARP***********************************************************************************
## NAME		buildFactorList, isFormalCovModel, catch
## PARAM		$summand - string
## AUTHOR		Sebastian Gross 
### DATE		29.08.2011
## changed by Martin Schlather 2015 -- 2017
## FUNCTION-END*************************************************************************************


buildFactorList <- function(summand, Env, ..x..) { #, last
  ## remove parenthesis
  factorA <- removeParenthesis(summand)

  ## formerly partially 'catch'
  X <- eval(parse(text=factorA), envir=Env)

  
  ##  Print(X, factorA, !isFormalCovModel(factorA))
  isFormalCovModel <- is(X, CLASS_CLIST) ||
    (regexpr("^[[:alnum:]_]+\\([[:print:]]*\\)$", factorA) == 1 &&
     exists(fun <- strsplit(factorA, "\\(")[[1]][1]) && is(get(fun), CLASS_RM))
  
  if (!isFormalCovModel) { ## (( && last))    
    if (is.factor(X)) {
      lev <- levels(X)
      if (length(lev) > MAXSUB^2 + 1)
        stop("max number of factors limited to ", MAXSUB^2)
      L <- list(RM_PLUS[1])
      i <- 2
      while (i <= length(lev)) {
        last <- min(length(lev), i + MAXSUB -1)
        plusList <- list(RM_PLUS[1])
        while (i <= last) {
          tmpList <- list(RM_COVARIATE)         
          tmpList[[COVARIATE_C_NAME]] <- as.numeric(X == lev[i])
          tmpList[[COVARIATE_X_NAME]] <- ..x..
          tmpList[[COVARIATE_ADDNA_NAME]] <- TRUE
          plusList[[length(plusList) + 1]] <- tmpList
          i <- i + 1;
        }
        L[[length(L) + 1]] <- plusList
      }
      return(if (length(lev) == 2) tmpList
             else if (length(lev) <= MAXSUB + 1) plusList
                                     else  L)
#            else list("RMtrend", mean=if (length(lev) <= MAXSUB + 1) plusList
#                                     else list("RMtrend", L)))
    } else if (is.vector(X)) {
      if (length(X) > 1) {
       tmpList <- list(RM_COVARIATE)         
       tmpList[[COVARIATE_C_NAME]] <- X
       tmpList[[COVARIATE_X_NAME]] <- ..x..
       tmpList[[COVARIATE_ADDNA_NAME]] <- TRUE       
      } else {
        tmpList <- list(R_CONST)
        tmpList[[CONST_A_NAME]] <-
          if (is.finite(X) && X == 1) NA else as.numeric(X)
      }
      return(tmpList) 
    } else stop(ERRMIXED)
    
  } else {
    tmpList <- buildCovList(X, x=..x..)
    return(tmpList)
  }
}

