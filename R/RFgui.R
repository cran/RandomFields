
## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
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


RFgui <- function(data, x, y,
                  same.algorithm = TRUE,
                  ev, bin = NULL,
                  xcov, ycov,
                  sim_only1dim=FALSE,
                  wait = 0, 
                  ...) {
  if (!interactive()) {
    warning("'RFgui' can be used only in an interactive mode")
    return(NULL)
  }
  wait <- as.integer(wait)
  Env <- if (wait >= 0) environment() else .GlobalEnv
  assign("RFgui.model", NULL, envir=Env)
  if (exists(".RFgui.exit", .GlobalEnv)) rm(".RFgui.exit", envir=.GlobalEnv)
 
  rfgui.intern(x=x, y=y, same.alg=same.algorithm,
                          ev=ev, bin=bin, xcov=xcov, ycov=ycov,
                          data=data, sim_only1dim=sim_only1dim, 
                          parent.ev = Env,
                          ...)
 
  if (wait >= 0) {
    while (!exists(".RFgui.exit", envir=Env))
      RandomFieldsUtils::sleep.micro(wait)
    res <- get("RFgui.model", envir=Env)
    if (is.null(res)) return(res) 
    if (RFoptions()$general$spConform) {
      RFvariogram(model=res, 0, get(".RFgui.y", envir=Env))
      res <- list2RMmodel(GetModel(RFvariogram))
    } else {
      class(res) <- "RM_model"
      invisible(res)
    }
  } else invisible(NULL)
}

rfgui.intern <- function(data, x, y,
                         same.alg = FALSE,
                         fixed.rs=NULL,                         
                         ev,
                         bin,
                         xcov, ycov,
                         sim_only1dim=FALSE,                         
                         parent.ev=NULL,
                         printlevel=0,...) {
  circ.trials <- 1
  circ.force <- TRUE
  circ.min <- -2

  if (missing(y)) y <- NULL
  if (missing(ev)) ev <- NULL
  if (missing(xcov)) xcov <- NULL
  if (missing(ycov)) ycov <- NULL
    
#  tcltk::tclRequire("BWidget", warn=!FALSE)

  ENVIR <- environment()
  assign("model", NULL, envir = ENVIR) # orignal: model als parameter uebergeben

  RFoptOld <- if (same.alg)
    internal.rfoptions(storing=FALSE, printlevel=printlevel - 10,
                       circulant.trials=circ.trials,
                       circulant.force=circ.force,
                       circulant.mmin=circ.min, ...,
                       graphics.height=-1
                       )
  else  internal.rfoptions(storing=FALSE, printlevel=printlevel - 10, ...)
  assign("RFopt.old", RFoptOld[[1]], envir=ENVIR)
  RFopt <- RFoptOld[[2]]
  rm("RFoptOld")

  guiReg <- MODEL_GUI
  guiOpt <- RFopt$gui

  tkDestroy <- tcltk::tkdestroy
  tkValue <- tcltk::tclvalue
  "tkValue<-" <- do.call("::", list("tcltk", "tclvalue<-"))
  tkLabel <- tcltk::tklabel
  tkEntry <- tcltk::tkentry
  tkScale <- tcltk::tkscale
  tkBind <- tcltk::tkbind
  tkGridConf <- tcltk::tkgrid.configure
  tkVar <- tcltk::tclVar
  Round <- function(x) base::round(x, digits=2)
  tkGrid <- tcltk::tkgrid
  tkPlot <- tkrplot::tkrplot
  tkCheckbutton <- tcltk::tkcheckbutton
  tkRadiobutton <- tcltk::tkradiobutton
  tkConfigure <- tcltk::tkconfigure
  Tcl <- tcltk::tcl
  tkRreplot <- tkrplot::tkrreplot
  tkButton <- tcltk::tkbutton
  tkCombobox <- tcltk::ttkcombobox
  tkRemove <- tcltk::tkgrid.remove
  

  
  OnModelSelected <- function(...)
  { 
    # delete die alten Parameterw*hler
    if(exists("baseModel", envir=ENVIR)) {
      baseParam <- get("baseModel", envir=ENVIR)$k
      if(length(baseParam) > 0) {
        for (i in 1:length(baseParam)) {
          n <- paste("slParam", i, "Value",   sep="") 
          baseParam[i] <- as.numeric(tkValue(get(n, envir=ENVIR)))
        }
        assign(paste("remember",selModelNum,sep=""), baseParam, envir=ENVIR)
        for (i in 1:length(baseParam)) {
          tkDestroy(get(paste("slParam", i, sep=""), envir=ENVIR))
          tkDestroy(get(paste("slParam", i, "Name", sep=""),envir=ENVIR))
          tkDestroy(get(paste("entryParam", i, sep=""), envir=ENVIR))
        }
      }
    }

    modelChoiceNum <- as.numeric(tkValue(Tcl(comboBox,"current")))
    if(modelChoiceNum == -1) return(0)

    # nun zum neuen Model
    modelChoice <- models[modelChoiceNum+1]
    selModelNum <- .C(C_GetModelNr, as.character(modelChoice), nr=integer(1),
		      PACKAGE="RandomFields")$nr

    selModelCountPar <- .C(C_GetNrParameters, selModelNum, k=integer(1),
                           PACKAGE="RandomFields")$k
    dim <- as.integer(2 - sim_only1dim)  
    newmodel <- list(modelChoice, k=rep(NA, times=selModelCountPar))
    modelParam <- try(.Call(C_SetAndGetModelInfo, guiReg,
                        list("Dummy", newmodel), dim,
                        FALSE, FALSE, FALSE, dim,
                        as.integer(10), ## ehemals RFoptions(short=10)
                        TRUE, TRUE, PACKAGE="RandomFields")$minmax)
    if (class(modelParam) == "try-error") return(0)
    
    assign("selModelNum",selModelNum, envir=ENVIR)
    if (exists("baseModel", where=ENVIR)) remove("baseModel", envir=ENVIR)
    if (selModelCountPar == 0) {
      assign("baseModel", list(modelChoice), ENVIR)
      Plot()
      return(0)
    }
    
   
    baseParam <- rep(NA, times=selModelCountPar)
    if(exists(paste("remember", selModelNum, sep=""), envir=ENVIR)) 
      baseParam <- get(paste("remember", selModelNum, sep=""), envir=ENVIR)

    ## selModelCountPar > 0 hier !!
##    openeps <- 1e-10
    for (i in 1:selModelCountPar) {
      baseParam[i] <- 
        if (!is.na(baseParam[i])) baseParam[i] else
        if (modelParam[i,3] == INTEGERPARAM) modelParam[i,2] else
        if (modelParam[i,1] >=0) 0.25 * sum(sqrt(modelParam[i,1:2]))^2 + 0.1
        else 0.5 * (modelParam[i,1] + modelParam[i,2])
 
      #Slider fuer den neuen Parameter 
      slParamValue <- tkVar(baseParam[i])
      entryParamValue <- tkVar(tkValue(slParamValue))
      # name <- unlist(strsplit(attr(modelParam, "dimnames")[[1]][i],"\\."))[2]
      # slParamName <- tkLabel(tt,text=paste(toupper(substring(name, 1,1)), substring(name, 2), sep=""))

      txt <- unlist(strsplit(attr(modelParam, "dimnames")[[1]][i],"\\."))[2]
      slParamName <- tkLabel(tt, text=txt)
      slParam <- tkScale(tt, command = Plot,
                         from= modelParam[i,1], 
                         to = modelParam[i,2],
                         showvalue=FALSE, variable=slParamValue,
                                ## neg value needed to get precise bounds:
                         resolution=if (modelParam[i, 3]==INTEGERPARAM) -1 else
                                    diff(modelParam[i,2:1])/numberSteps, 
                         orient="horizontal", length=length.slider, width=18)
      entryParam <- tkEntry(tt,width=size.entry,textvariable=entryParamValue)
      tkBind(entryParam, "<Return>", OnAddParamEntryChanged)
      
      assign("slParam", slParam, envir=ENVIR)
      assign(paste("slParam", i, sep=""), slParam, envir=ENVIR)
      assign(paste("slParam", i, "Name", sep=""), slParamName, envir=ENVIR)
      assign(paste("slParam", i, "Value", sep=""), slParamValue, envir=ENVIR)
      assign("entryParam", entryParam, envir=ENVIR)
      assign(paste("entryParam", i, sep=""), entryParam, envir=ENVIR)
      assign(paste("entryParam", i, "Value", sep=""), entryParamValue,
             envir=ENVIR)
    }
    baseModel <- list(modelChoice, k=baseParam)

    assign("baseModel", baseModel, ENVIR)
    position()
  }

  OnPlotVarCovChanged <- function(...)
  {
    if((as.character(tkValue(plotVarCov)) == "Variogram") && !is.null(ev)) {
      #Print("here", cbPlotEV);
      
      ## tkConfigure(cbPlotEV, disabled=FALSE)
      #Print("here end")
    } else {
      tkValue(plotEV) <- "0"
    }
  #    tkConfigure(cbPlotEV, disabled=TRUE)
    
    tkRreplot(imgVar)
  }

  OnplotEVChanged <- function(...)
  {
    if((as.character(tkValue(plotVarCov)) == "Covariance") || is.null(ev))
      tkValue(plotEV) <- "0"
    tkRreplot(imgVar)
  }

  
  EntryChanges <- function(var, value, strictpos=TRUE, factor=2) {
    value <- as.numeric(tkValue(value))
    if (is.na(strictpos)) {
      to <- if (value > 0) value * factor else value / factor
      from <- if (value < 0) value * factor else value / factor
    } else {
      if (value < 0) stop("negative values not allowed")
      if (strictpos) {
        value <- log(value)
        to <- value + log(factor)
        from <- value - log(factor)
      } else {
        to <- if (value == 0) 1 else value * factor
        from <- if (value > 1) value / factor else 0
      }
    }
    resolution <-  0.01 * (to - from)
    tkConfigure(var, to=to, from = from, resolution = -resolution)    
    return(value)
  }
 
  
  OnScaleEntryChanged <- function(...) { 
    tkValue(slScaleValue) <- EntryChanges(slScale, entryScaleValue)
    Plot()
  }


  OnVarEntryChanged <- function(...) {
    tkValue(slVarianceValue) <- EntryChanges(slVariance,
                                                     entryVarianceValue)
    Plot()
  } 

   OnNuggetEntryChanged <- function(...) {
     tkValue(slNuggetValue) <-
       EntryChanges(slNugget, entryNuggetValue, strict=FALSE)
     Plot()
   }  

  OnRotationEntryChanged <- function(...)
  { 
    tkValue(slRotationValue) <- as.numeric(tkValue(entryRotationValue))
    Plot()
  }  

  OnRadiusEntryChanged <- function(...)
  {
    tkValue(slScaleAValue) <- EntryChanges(slScaleA, entryScaleAValue)
    tkValue(slScaleBValue) <- EntryChanges(slScaleB, entryScaleBValue)
    Plot()
  }  

  OnAddParamEntryChanged <- function(...)
  { 
    baseModel <- get("baseModel", envir=ENVIR)
    if(length(baseModel$k) > 0)
      for (i in 1:length(baseModel$k)) {
        slParamValue <- get(paste("slParam", i, "Value", sep=""), envir=ENVIR)
        value <- get(paste("entryParam", i, "Value", sep=""), envir=ENVIR)
        tkValue(slParamValue) <-
          Round(as.numeric(tkValue(value)))
      }
    Plot()
  } 
    
  GetGuiModel <- function() {
    variance <- exp(as.numeric(tkValue(slVarianceValue)))
    nugget <- as.numeric(tkValue(slNuggetValue))
    baseParam <- baseModel$k
    if(length(baseModel$k) > 0)
      for (i in 1:length(baseModel$k)) { 
        baseParam[i] <- as.numeric(tkValue(get(paste("slParam", i, "Value",
                                                      sep=""), envir=ENVIR)))
        entryParamValue <-
          get(paste("entryParam", i, "Value", sep=""), envir=ENVIR)
        tkValue(entryParamValue) <- Round(baseParam[i])
      }

    baseModel$k <- baseParam
   
    if(!as.numeric(tkValue(showAniso))) {
      scale <- exp(as.numeric(tkValue(slScaleValue)))
      newmodel <- list(ZF_SYMBOLS_PLUS,
                    list(DOLLAR[1], var=variance, scale=scale, baseModel),
                    list(DOLLAR[1], var=nugget, list(ZF_NUGGET[1])))
    } else {
      a <-  as.numeric(tkValue(slRotationValue))
      r <- c(exp(as.numeric(tkValue(slScaleAValue))),
             exp(as.numeric(tkValue(slScaleBValue))))
      u <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), ncol=2 )
      aniso <- u %*% (1/r * t(u))
      newmodel <- list(ZF_SYMBOLS_PLUS,
                    list(DOLLAR[1], var=variance, aniso=aniso, baseModel),
                    list(DOLLAR[1], var=nugget, list(ZF_NUGGET[1])))
    }
    return(newmodel)
  }

  plotFunction <- function(...)
  {

    #Print(tkValue(plotEV), cbPlotEV)

    plotev = as.numeric(tkValue(plotEV))
    par(cex=0.6, bg="lightgrey", mar=c(3,3,1,1))
    if(!exists("baseModel",envir=ENVIR)) {
      if(!is.null(ev) && plotev) {
        notNA <- !is.nan(ev@emp.vario)
        xm <- c(min(ev@centers[notNA]), max(ev@centers[notNA]))
        ym <- c(min(ev@emp.vario[notNA]), max(ev@emp.vario[notNA])*1.1)
        lab <- xylabs("", NULL)
        plot(ev@centers[!is.nan(ev@emp.vario)],
             ev@emp.vario[!is.nan(ev@emp.vario)], pch=19, xlab=lab$x) 
        return(0)
      } 
      plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
      return(0)
    }

    #baseModel <- get("baseModel",envir=ENVIR)
    tkValue(entryScaleValue) <-
      Round(exp(as.numeric(tkValue(slScaleValue))))
    tkValue(entryVarianceValue) <-
      Round(exp(as.numeric(tkValue(slVarianceValue))))
    tkValue(entryNuggetValue) <-
      Round(as.numeric(tkValue(slNuggetValue)))
    tkValue(entryScaleAValue) <-
      Round(exp(as.numeric(tkValue(slScaleAValue))))
    tkValue(entryScaleBValue) <-
      Round(exp(as.numeric(tkValue(slScaleBValue))))
    tkValue(entryRotationValue) <-
      Round(as.numeric(tkValue(slRotationValue)))


    newmodel <- GetGuiModel()
    assign("RFgui.model", newmodel, envir=parent.ev)
    
     if(as.numeric(tkValue(showAniso))) {
      x1 <- rep(xcov, each=length(ycov))
      x2 <- rep(ycov, times=length(xcov))

      cv <- RFvariogram(x=as.matrix(expand.grid(xcov, ycov)),
                        model=newmodel, 
                        practicalrange = tkValue(cbPracRangeVal) != "0")
      dim(cv) <- c(length(ycov),length(xcov))

      cv00 <- cv[1,1]
      if (xcov[1] == 0 && ycov[1] == 0) {
        zlim <- c(0, 1.1*max(cv))
        cv[1,1] <- NA
      }
      tranMatrix <- persp(x=xcov, y=ycov, z=cv,
                          theta = as.numeric(tkValue(slTurnPlotValue)),
                          zlim = zlim, phi = 0, xlab = "x", ylab = "y",
                          zlab = as.character(tkValue(plotVarCov)),
                          col = "lightblue", ltheta = 120, shade = 0.75,
                          ticktype = "detailed")
      if (xcov[1] == 0 && ycov[1] == 0)
        points(trans3d(xcov[1], ycov[1], cv00, pmat = tranMatrix), pch =16)

      assign("model", newmodel, envir = ENVIR)
      return(0)
    }
    
    cv <- xcov
    if(as.character(tkValue(plotVarCov)) == "Covariance") {
     
      cv <- RFcov(x=xcov, model=newmodel,
                  practicalrange = tkValue(cbPracRangeVal) != "0")
    }
    if(as.character(tkValue(plotVarCov)) == "Variogram") {      
      pr.dummy <- tkValue(cbPracRangeVal) != "0"
 
      cv <- RFvariogram(x=xcov, model=newmodel,
                        practicalrange = pr.dummy)
     }

    if(!is.null(ev) && plotev) {
      xm <- range(ev@centers, na.rm=TRUE)
      ym <- range(ev@emp.vario, na.rm=TRUE) * c(1, 1.1)
    } else {
      xm <- range(xcov, na.rm=TRUE)
      ym <- range(cv, na.rm=TRUE) * c(1, 1.1)
    }

    lab <- xylabs("", NULL)
    plot(xcov[2:length(xcov)], cv[2:length(xcov)], type="l",
         xlab=lab$x, ylab="", xlim=xm, ylim=ym)
    points(xcov[1], cv[1])

    # plot emp.vario
    if(!is.null(ev) && plotev)    
      points(ev@centers[!is.nan(ev@emp.vario)],
             ev@emp.vario[!is.nan(ev@emp.vario)], pch=19)

     assign("model", newmodel, envir = ENVIR)
  } # function

  plotSimulation <- function(...) {
    par(cex=0.6, bg="lightgrey", mar=c(3,3,1,1))
    if(!exists("model", envir=ENVIR)) {
      plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
      return(0)
    }

    simu.model <- get("model", envir=ENVIR)
    if (!is.null(simu.model)) {
      par(cex=0.6, bg="lightgrey")
      if (guiOpt$simu_method != "any method")
        simu.model <- list(guiOpt$simu_method, simu.model)
      
      yy <- (if (get("simDim", envir = ENVIR) =="sim1Dim") NULL else
             if (length(y)==0) x else y)
      pr <-  tkValue(cbPracRangeVal) != "0"
      z <- try(RFsimulate(simu.model,x=x, grid=TRUE, 
                          y=if (get("simDim", envir = ENVIR)=="sim1Dim") NULL
                          else if (length(y)==0) x else y,
                          seed = fixed.rs,
                          register=guiReg, spConform=TRUE,
                          practicalrange =
                          tkValue(cbPracRangeVal) != "0"),
               silent=!TRUE)
 
     if (class(z) == "try-error") {
         plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=!FALSE, xlab="",
               ylab="",
             cex.main=2, col.main="brown",
              main=paste("\n\n\n\n\n\n\n\n",
               if (guiOpt$simu_method == "any method")
               "No suitable simulation algorithm found"
               else paste("Simulation method '", guiOpt$simu_method,
                          "'\ndoes not work for this specification", sep=""),
               ".\n\nSet 'simu_method = \"any method\"'",
               "\nor set 'same.algorithm=FALSE'",
              "\nor see RFoptions() for controlling parameters",
               if (guiOpt$simu_method != "any method")
               paste("\nof '", guiOpt$simu_method, "'", sep=""),
               sep=""))
     } else {
        plot(z, cex=.5, legend=FALSE, xlab=NULL)
     }
    }
  }

  
  Plot <- function(...) {    
#
    #tkConfigure(labelOccupancy,textvariable=tkVar("Busy"))    
    tkRreplot(imgVar)
    if (as.numeric(tkValue(simAlways))) tkRreplot(imgSim)
   #tkConfigure(labelOccupancy,textvariable=tkVar("Free"))
  }

  OnChangeIsotropie <- function(...)
  {
    if(as.numeric(tkValue(showAniso))) {
      tkRemove(slScale)
      tkRemove(entryScale)
      tkRemove(labelScale)
    }else {
      tkRemove(slScaleA)
      tkRemove(entryScaleA)
      tkRemove(labelScaleA)
      tkRemove(slScaleB)
      tkRemove(entryScaleB)
      tkRemove(labelScaleB)
      tkRemove(slRotation)
      tkRemove(entryRotation)
      tkRemove(labelRotation)
      tkRemove(slTurnplot)
    }
    position()
    Plot()
  }

  OnTurnPlot <- function(...)
  {
    tkRreplot(imgVar)
  }

  OnNewSimu <- function(...) 
  {
    assign("fixed.rs", Round(runif(1,1,100000)), envir=ENVIR)  
    tkRreplot(imgSim)
  }

  OnSimDimChanged <- function(...)
  {
  
    if(!sim_only1dim) {
      assign("simDim", tkValue(rb2DimValue), envir = ENVIR) 
      tkRreplot(imgSim)
      return (0)
    }

    if(as.numeric(tkValue(showAniso))) 
    {
      tkValue(rb2DimValue) <-"sim2Dim"
      return(0)
    } 
  }

  OnReturn <- function(...)
  {
    # hier muss eine rueckgabe stehen, emp vario und model mit parametern
   #   Print(GetGuiModel())
    RFoptions(LIST=get("RFopt.old", envir=ENVIR))
    ##remove("RFopt.old", envir=ENVIR)
    assign(".RFgui.exit", TRUE, envir=parent.ev)
    tkDestroy(tt)    
  }

  position <- function(...) {  
    #--- PLOT  ---------------------------------------------------------
    tkGridConf(imgVar, rowspan=image.rowspan, columnspan=image.colspan,
               column=col.var, row=1,sticky="w") 
    tkGridConf(imgSim, rowspan=image.rowspan, columnspan=image.colspan,
               column=col.sim, row=1,sticky="w")

    #--- DropDown-ComboBox for model selection -------------------------
    tkGridConf(labModelSelect, column=col.sl, row=row.sl)
    row.sl <- row.sl+1
    tkGridConf(comboBox, column=col.sl, row=row.sl, sticky = "e")
    row.sl <- row.sl+1
 
    #--- Radiobutton zur Frage Variogram oder Covarianzfunktion --------
    tkGridConf(rbCovariance, column=col.var+image.colspan-1,
                     row=image.rowspan+1, sticky="w")
    tkGridConf(labelCovariance, column=col.var+image.colspan-1,
                     row=image.rowspan+1, sticky="e")
    tkGridConf(rbVariogram, column=col.var+image.colspan-1,
                     row=image.rowspan+2, sticky="w") 
    tkGridConf(labelVariogram, column=col.var+image.colspan-1,
                     row=image.rowspan+2, sticky="e")

    #--- Checkbox show the empirical variogram --------------------------
    tkGridConf(cbPlotEV, column=col.var, row=image.rowspan+1, sticky="w")
    tkGridConf(labelPlotEV, column=col.var, row=image.rowspan+1,
                     sticky="e")
    
    #--- Radiobutton: select dimension for simulation -------------------
    tkGridConf(rbSim1Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+1, sticky="w")
    tkGridConf(labelSim1Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+1, sticky="e")
    tkGridConf(rbSim2Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+2, sticky="w") 
    tkGridConf(labelSim2Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+2, sticky="e")

    #--- Checkbox simulate on slider movement --------------------------
    tkGridConf(cbSimAlways, column=col.sim, row=image.rowspan+1,
                     sticky="w")
    tkGridConf(labelSimAlways, column=col.sim, row=image.rowspan+1,
                     sticky="e")
  
    #--- Checkboxes Practical Range  and anisotropy option -------------
    if(length(y)!=0) {
      tkGridConf(cbAnisotropy, column=col.sl, row=row.sl, sticky="w")
      tkGridConf(labelAniso, column=col.sl, row=row.sl)
      row.sl=row.sl+1
    }
  
    tkGridConf(cbPracRange, column=col.sl, row=row.sl, sticky="w")
    tkGridConf(labelPracRange, column=col.sl, row=row.sl)
    row.sl=row.sl+1

    #--- Parameterwaehler ----------------------------------------------
   
    if(!as.numeric(tkValue(showAniso))) {
      #Slider Scale
      tkGridConf(labelScale,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkGridConf(slScale,  column=col.sl, row=row.sl, sticky="w")
      tkGridConf(entryScale, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1
    }else {
      #Slider Rotation
      tkGridConf(labelRotation,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkGridConf(slRotation,  column=col.sl, row=row.sl, sticky="w")
      tkGridConf(entryRotation, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1

      #Slider Radius
      tkGridConf(labelScaleA,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkGridConf(slScaleA,  column=col.sl, row=row.sl, sticky="w")
      tkGridConf(entryScaleA, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1

      tkGridConf(labelScaleB,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkGridConf(slScaleB,  column=col.sl, row=row.sl, sticky="w")
      tkGridConf(entryScaleB, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1

      #Slider turn the now 2dim covarianz plot
      tkGridConf(slTurnplot, column=col.var,columnspan=image.colspan,
                       row=image.rowspan+3)
    }

    #Slider variance
    tkGridConf(labelVariance,  column=col.sl, row=row.sl)
    row.sl <- row.sl+1
    tkGridConf(slVariance,  column=col.sl, row=row.sl, sticky="w")
    tkGridConf(entryVariance, column=col.sl, row=row.sl, sticky="e")
    row.sl <- row.sl+1

    #Slider nugget
    tkGridConf(labelNugget,  column=col.sl, row=row.sl)
    row.sl <- row.sl+1
    tkGridConf(slNugget,  column=col.sl, row=row.sl, sticky="w")
    tkGridConf(entryNugget, column=col.sl, row=row.sl, sticky="e")
    row.sl <- row.sl+1

    if(exists("baseModel", envir=ENVIR)) {
      baseModel <- get("baseModel",envir=ENVIR)
      baseParam <- baseModel$k
      if(length(baseModel$k) > 0)
        for (i in 1:length(baseModel$k)) { 
          tkGridConf(get(paste("slParam",i,"Name",sep=""),
                               get("slParam", envir=ENVIR),
                               envir=ENVIR), column=col.sl, row=row.sl)
          row.sl <- row.sl+1 
          tkGridConf(get(paste("slParam", i, sep=""),
                               get("slParam", envir=ENVIR),
                               envir=ENVIR),
                           column=col.sl, row=row.sl, sticky="w")  
          tkGridConf(get(paste("entryParam", i, sep=""),
                               get("entryParam", envir=ENVIR),
                               envir=ENVIR),
                           column=col.sl, row=row.sl, sticky="e")
          row.sl <- row.sl+1
        }
    }

    #--- Buttons - new simulation (new seed), return ---------------
    tkGridConf(buttonNewSimu, column=col.sl,
                     row=max(row.sl,image.rowspan+1), sticky="e")
    row.sl=row.sl+1
    tkGridConf(buttonReturn, column=col.sl,
                     row=max(row.sl,image.rowspan+2), sticky="e")
    row.sl=row.sl+1
    #--- Beschaeftigungsindikator ---------------------------------------
#    tkGridConf(labelOccupancy, row=row.last, column=col.sl, sticky="e")    
  } ## end fct position


 
  PRINT <- FALSE

  if (missing(data)) {
    if (missing(x)) x <- seq(1, 5, len=guiOpt$size[if (sim_only1dim) 1 else 2] )
  } else {
    S <- # if (exists("UnifyData")) UnifyData(data=data, RFopt=RFopt) else
         StandardizeData(data=data, RFopt=RFopt)
    if (S$matrix.indep.of.x.assumed)
      stop("data must contain the information about the locations of the data")
    else {
      if (!missing(x)) message("coordinates detected in the data")
    }
    if (missing(x)) {
      r <- apply(S$coord[[1]]$x, 2, range)    
      len <- guiOpt$size[if (length(r) == 2) 1 else 2] 
      x <- seq(r[1], r[2], len=len)
      if (length(r) > 2) {
         y <- seq(r[3], r[4], len=len)
      } else sim_only1dim <- TRUE
    }
  }
  assign(".RFgui.y", if (length(y)==0) NULL else 0, envir=parent.ev)

  if(!missing(data) && !is.null(data)) {
    if (!is.null(ev)) stop("if 'data' is given, 'ev' may not be given.")   
    ev <- RFempiricalvariogram(data=data, phi=1, bin=bin, vdim=1)
  }
  
  if (any(diff(x) <= 0)) 
    stop("x should be a sequence of increasing numbers")

  if (length(y)!=0 && any(diff(y) <= 0))
    stop("y should be a sequence of increasing numbers")
  
  if(is.null(xcov)) {
    if(is.null(ev))
      xcov <- seq(0,15,0.1)
    else
      xcov <-seq(min(0,0.9*min(ev@centers)), max(ev@centers),
                 by=diff(range(ev@centers))/100)
  }
  
  if(is.null(ycov) && length(y)!=0) {
    if(is.null(ev))
      ycov <- seq(0,15,by=0.1)
    else
      ycov <-seq(min(0, 0.9*min(ev@centers)), max(ev@centers),
                 by=diff(range(ev@centers))/100)
  }


  if (exists("baseModel", where=ENVIR)) remove("baseModel", envir=ENVIR)
             
  if(is.null(fixed.rs)) {
    if (!exists(".Random.seed")) runif(1)
    fixed.rs <- .Random.seed 
  }

  # get all model names
  models <- if (sim_only1dim) rfgui_Names1 else rfgui_Names2
  models <- models[models != "RMnugget"]
  
  #-------------------------------------------------------------------
  # Start Values and ranges
  #-------------------------------------------------------------------
  cbPracRangeVal <- tkVar(RFopt$general$practicalrange)
  simAlways <- tkVar(as.integer(guiOpt$alwaysSimulate))
  plotVarCov <- tkVar("Variogram")

  plotEV <- tkVar(ifelse(is.null(ev) && tkValue(plotVarCov)=="Variogram",
                          "0", "1"))
  showAniso <- tkVar("0")
  slTurnPlotValue <- tkVar("0")
  numberSteps <- 50

  if (is.null(ev)) {
    ## Nugget  
    nugget <- 0
    nuggetMin <- 0
    nuggetMax <- 10
    ## Variance
    variance <- 1
    ## Scale
    scale <- 1
  } else {
    ## Nugget  
    idx1 <- !is.nan(ev@emp.vario)
    idx2 <- 2:min(6,length(is.nan(ev@emp.vario)))
    nugget <- max(0, lm(ev@emp.vario[idx1][idx2]
                        ~ ev@centers[idx1][idx2])$coefficients[1])
    nuggetMin <- 0 ## nugget/10
    nuggetMax <- max(ev@emp.vario[idx1])
    ## Variance
    variance <- quantile(ev@emp.vario[idx1], probs=0.7)
     ## Scale
    scale <-  0.3*max(ev@centers)
  }
  ## Nugget  
  slNuggetValue <- tkVar(nugget)
   ## Variance
  varianceMin <- Round(log(0.01))
  varianceMax <- log(max(1e-10, nuggetMax))
  slVarianceValue <- tkVar(log(variance))
  ## Scale
  scaleMin <- Round(log(0.1*scale))  
  scaleMax <- Round(log(10*scale))   
  slScaleValue <- tkVar(log(scale))
  
  ## die direkte eingabe muss als variable getrennt von den schiebern laufen
  entryScaleValue <- tkVar(scale)
  entryVarianceValue <- tkVar(variance)
  entryNuggetValue <- tkVar(nugget) 
  ## die direkte eingabe muss als variable getrennt von den schiebern laufen
  entryScaleValue <- tkVar(scale)
  entryVarianceValue <- tkVar(variance)
  entryNuggetValue <- tkVar(nugget) 

  # die anisotropie tierchen
  slRotationValue <- tkVar("0")
  entryRotationValue <- tkVar("0")
  anisoScale = "1"
  slScaleAValue <- tkVar(anisoScale)
  slScaleBValue <- tkVar(anisoScale)
  entryScaleAValue <- tkVar(anisoScale)
  entryScaleBValue <- tkVar(anisoScale)
  radiusMax <- 2

  #------------------------------------------------------------------
  # GUI
  #------------------------------------------------------------------
  tt <- tcltk::tktoplevel()
  tcltk::tktitle(tt) <- "U Diffusion Gui"
  tcltk::tkwm.protocol(tt, "WM_DELETE_WINDOW", OnReturn)
  
  # some position variables
  # assuming all sliders are in the same column and consecutive rows
  image.rowspan <- 15
  image.colspan <- 6
  col.sim <- 1
  col.var <- col.sim+image.colspan
  col.sl <- col.var+image.colspan
  row.last <- image.rowspan+4
  image.colspan <- 5
  row.sl <- 1 # giving the position of the first label
  size.entry <- 4
  length.slider <- 130
  # size of the variogram/cf and simulation plot
  plothscale <- 0.8    # Horizontal scaling
  plotvscale <- 0.8    # Vertical scaling

  tkGrid(tkLabel(tt, text="", width=1, height=0), column=0, row=0)
  tkGrid(tkLabel(tt, text="", width=1), column=col.sim+image.colspan,row=1)
  tkGrid(tkLabel(tt, text="", width=1), column=col.var+image.colspan,row=1)
  tkGrid(tkLabel(tt, text="", width=1), column=col.sl+image.colspan, row=1)

  #--- DropDown-ComboBox for model selection -------------------------
  labModelSelect <- tkLabel(tt,text="Model Selection")
  textModell <- tkVar("Please select a model...")
  comboBox <- tkCombobox(tt,textvariable=textModell, state="readonly",
                          values=models)
  tkBind(comboBox, "<<ComboboxSelected>>", OnModelSelected)

  #--- PLOT  ---------------------------------------------------------
  imgVar <- tkPlot(tt,fun=plotFunction,hscale=plothscale,vscale=plotvscale)
  imgSim <- tkPlot(tt,fun=plotSimulation,hscale=plothscale,vscale=plotvscale)
    
  #--- Beschaeftigungsindikator -------------------------------------
 # labelOccText <- tkVar("Free")
 # labelOccupancy <- tkLabel(tt,text=tkValue(labelOccText))

  #--- anistropy -----------------------------------------------------
  cbAnisotropy <- tkCheckbutton(tt, variable=showAniso,
                                command=OnChangeIsotropie)
  labelAniso <-  tkLabel(tt,text="Anisotropy")

  #--- Radiobutton zur Frage Variogram oder Covarianzfunktion --------
  rbVariogram <- tkRadiobutton(tt, command=OnPlotVarCovChanged)
  rbCovariance <- tkRadiobutton(tt, command=OnPlotVarCovChanged)
  tkConfigure(rbVariogram,variable=plotVarCov,value="Variogram")
  tkConfigure(rbCovariance,variable=plotVarCov,value="Covariance")
  labelCovariance <- tkLabel(tt,text="Covariance")
  labelVariogram <- tkLabel(tt,text="Variogram")

  #--- Checkbox plot empirical variogram    --------------------------
  ## checkbuttion setzt die variable und fuehrt dann noch zusaetzlich
  ## command aus.
  cbPlotEV <- tkCheckbutton(tt, variable=plotEV, command=OnplotEVChanged)
  labelPlotEV <- tkLabel(tt,text="Plot empirical variogram")

  #--- Radiobutton: select dimension for simulation ------------------
  rbSim1Dim <- tkRadiobutton(tt, command=OnSimDimChanged)
  rbSim2Dim <- tkRadiobutton(tt, command=OnSimDimChanged)
  rb2DimValue <- tkVar(if (sim_only1dim) "sim1Dim" else "sim2Dim")
  tkConfigure(rbSim1Dim,variable=rb2DimValue, value="sim1Dim")
  tkConfigure(rbSim2Dim,variable=rb2DimValue, value="sim2Dim")
  labelSim1Dim <- tkLabel(tt,text="1 dim")
  labelSim2Dim <- tkLabel(tt,text=if (sim_only1dim) "1 dim" else "2 dim")
  assign("simDim", as.character(tkValue(rb2DimValue)), envir=ENVIR)

  #--- Button - new simulation (new seed) ----------------------------
  buttonNewSimu <- tkButton(tt,text="New Simulation",command=OnNewSimu)
  
  #--- Button - Return -----------------------------------------------
  buttonReturn <- tkButton(tt,text="      Return       ",command=OnReturn)

  #--- Checkbox simulate on slider movement --------------------------
  cbSimAlways <- tkCheckbutton(tt, variable=simAlways)
  labelSimAlways <- tkLabel(tt,text="Simulate always")

  #--- Checkbox Practical Range --------------------------------------
  cbPracRange <- tkCheckbutton(tt, variable=cbPracRangeVal, command=Plot)
  tkConfigure(cbPracRange,variable=cbPracRangeVal)
  labelPracRange <- tkLabel(tt,text="Practical Range")

  #--- Slider turn covariance plot in anisotropic case ----------------
  slTurnplot <- tkScale(tt, command = OnTurnPlot, from=0, to=360,
                        showvalue=TRUE, variable=slTurnPlotValue,
                        resolution=-360/numberSteps, orient="horizontal",
                        length=2*length.slider, width=18)
 # labelRotation <- tkLabel(tt,text="Rotation")

  #--- Parameterwaehler ----------------------------------------------
  # Aniso Rotation
  slRotation <- tkScale(tt, command = Plot, from=0, to=0.5*pi,
                        showvalue=FALSE, variable=slRotationValue,
                        resolution=-0.5*pi/numberSteps, orient="horizontal",
                        length=length.slider, width=18)
  labelRotation <- tkLabel(tt,text="Rotation")
  entryRotation <- tkEntry(tt,width=size.entry, textvariable=entryRotationValue)
  tkBind(entryRotation, "<Return>", OnRotationEntryChanged)

  # Aniso Scales in first and second axis
  slScaleA <- tkScale(tt, command = Plot, from=scaleMin, to=scaleMax,
                      showvalue=FALSE, variable=slScaleAValue,
                      resolution=-radiusMax/numberSteps, orient="horizontal",
                      length=length.slider, width=18)
  labelScaleA <- tkLabel(tt,text="first axis scale")
  entryScaleA <- tkEntry(tt,width=size.entry, textvariable=entryScaleAValue)
  tkBind(entryScaleA, "<Return>", OnRadiusEntryChanged)

  slScaleB <- tkScale(tt, command = Plot, from=scaleMin, to=scaleMax,
                      showvalue=FALSE, variable=slScaleBValue,
                      resolution=-radiusMax/numberSteps, orient="horizontal",
                      length=length.slider, width=18)
  labelScaleB <- tkLabel(tt,text="second axis scale")
  entryScaleB <- tkEntry(tt,width=size.entry, textvariable=entryScaleBValue)
  tkBind(entryScaleB, "<Return>", OnRadiusEntryChanged)

  # Scale
  slScale <- tkScale(tt, command = Plot, from=scaleMin, to=scaleMax,
                            showvalue=FALSE, variable=slScaleValue,
                            resolution=-(scaleMax-scaleMin)/numberSteps,
                            orient="horizontal", length=length.slider, width=18)
  labelScale <- tkLabel(tt,text="Scale")
  entryScale <- tkEntry(tt,width=size.entry, textvariable=entryScaleValue)
  tkBind(entryScale, "<Return>", OnScaleEntryChanged)

  #Slider variance
  slVariance <- tkScale(tt, command = Plot, from=varianceMin,
                               to=varianceMax,
                               showvalue=FALSE, variable=slVarianceValue,
                               resolution=(varianceMin-varianceMax)/numberSteps,
                               orient="horizontal", length=length.slider, width=18)
  labelVariance <- tkLabel(tt,text="Variance")
  entryVariance <- tkEntry(tt,width=size.entry, textvariable=entryVarianceValue)
  tkBind(entryVariance, "<Return>", OnVarEntryChanged)

  #Slider nugget
  slNugget <- tkScale(tt, command = Plot, from=nuggetMin, to=nuggetMax,
                             showvalue=FALSE, variable=slNuggetValue,
                             resolution=-(nuggetMax-nuggetMin)/numberSteps,
                             orient="horizontal", length=length.slider, width=18)
  labelNugget <- tkLabel(tt,text="Nugget")
  entryNugget <- tkEntry(tt,width=size.entry, textvariable=entryNuggetValue)
  tkBind(entryNugget, "<Return>", OnNuggetEntryChanged)

  position()
}
