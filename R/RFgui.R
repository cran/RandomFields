

RFgui <- function(data, x, y,
                  same.algorithm = FALSE,
                  ev, 
                  xcov, ycov,
                  sim_only1dim=FALSE,
                  wait = 0,
                  ...) {
  if (!interactive())
    return("'RFgui' might be used only in an interactive mode")
  wait <- as.integer(wait)
  Env <- if (wait >= 0) environment() else .GlobalEnv
  assign("RFgui.model", NULL, envir=Env)
  if (exists(".RFgui.exit", .GlobalEnv)) rm(".RFgui.exit", envir=.GlobalEnv)

  #Print(same.algorithm)
  
  rfgui.intern(x=x, y=y, same.alg=same.algorithm,
               ev=ev, xcov=xcov, ycov=ycov,
               data=data, sim_only1dim=sim_only1dim, 
               parent.ev = Env,
              ...)
  
  if (wait >= 0) {
    while (!exists(".RFgui.exit", envir=Env)) .C("sleepMicro", as.integer(wait))
    res <- get("RFgui.model", envir=Env)
    if (RFoptions()$general$spConform) { 
      RFvariogram(model=res, x=0)
      res <- list2RMmodel(GetModel(RFvariogram))
    }
    invisible(res)
  } else invisible(NULL)
}

rfgui.intern <- function(data, x, y,
                         same.alg = FALSE,
                         fixed.rs=NULL,                         
                         ev, 
                         xcov, ycov,
                         sim_only1dim=FALSE,                         
                         parent.ev=NULL, ...) {
  simuMethod <- "circulant"
  printlevel <- 0
  circ.trials <- 1
  circ.force <- TRUE
  circ.min <- -2

  if (missing(y)) y <- NULL
  if (missing(ev)) ev <- NULL
  if (missing(xcov)) xcov <- NULL
  if (missing(ycov)) ycov <- NULL
  
  if (!do.call("require", list("tcltk", quietly=TRUE)) ||
      !do.call("require", list("tkrplot", quietly=TRUE)))
    stop("To use 'RFgui' the tcl/tk packages 'tcltk' and 'tkrplot' must be installed.")
  
  tclRequire("BWidget", warn=!FALSE)

  ENVIR <- environment()
  assign("model", NULL, envir = ENVIR) # orignal: model als parameter uebergeben

  RFoptOld <- if (same.alg)
    internal.rfoptions(storing=FALSE, printlevel=printlevel,
                       circulant.trials=circ.trials,
                       gui.simu_method = simuMethod,
                       circulant.force=circ.force,
                       circulant.mmin=circ.min, ...)
  else  internal.rfoptions(storing=FALSE, printlevel=printlevel, ...)
  #Print(same.alg, RFoptOld)
  assign("RFopt.old", RFoptOld[[1]], envir=ENVIR)
  RFopt <- RFoptOld[[2]]
  rm("RFoptOld")
  guiReg <- GetModelRegister("gui")
  guiOpt <- RFopt$gui
  stopifnot(guiReg == RFopt$general$guiregister)
  
  OnModelSelected <- function(...)
  { 
    # delete die alten Parameterwähler
    if(exists("baseModel", envir=ENVIR)) {
      baseParam <- get("baseModel", envir=ENVIR)$k
      if(length(baseParam) > 0) {
      for (i in 1:length(baseParam)) 
        baseParam[i] <- as.numeric(tclvalue(get(paste("slParam", i, "Value",
                                                      sep=""), envir=ENVIR)))
        assign(paste("remember",selModelNum,sep=""), baseParam, envir=ENVIR)
        for (i in 1:length(baseParam)) {
          tkdestroy(get(paste("slParam", i, sep=""), envir=ENVIR))
          tkdestroy(get(paste("slParam", i, "Name", sep=""), envir=ENVIR))
          tkdestroy(get(paste("entryParam", i, sep=""), envir=ENVIR))
        }
      }
    }

    if (exists("baseModel", where=ENVIR)) remove("baseModel", envir=ENVIR)
    modelChoiceNum <- as.numeric(tclvalue(tcl(comboBox,"current")))
    if(modelChoiceNum == -1) 
      return(0)

    # nun zum neuen Model
    modelChoice <- models[modelChoiceNum+1]
      selModelNum <- .C("GetModelNr", as.character(modelChoice), nr=integer(1),
		      PACKAGE="RandomFields")$nr
    assign("selModelNum",selModelNum, envir=ENVIR)
  
    selModelCountPar <- .C("GetNrParameters", selModelNum, k=integer(1),
                           PACKAGE="RandomFields", DUP = FALSE)$k

    ##Print(selModelCountPar, newmodel)
    
    if (selModelCountPar == 0) {
      assign("baseModel", list(modelChoice), ENVIR)
      plotDensity()
      return(0)
    }
    
    dim <- as.integer(2 - sim_only1dim)     
    newmodel <- list(modelChoice, k=rep(NA, times=selModelCountPar))
    modelParam <- .Call("SetAndGetModelInfo", guiReg,
                        list("Dummy", newmodel), dim,
                        FALSE, FALSE, FALSE, dim,
                        as.integer(10), ## ehemals RFoptions(short=10)
                        TRUE, TRUE, PACKAGE="RandomFields")$minmax
   
    baseParam <- rep(NA, times=selModelCountPar)
    if(exists(paste("remember",selModelNum,sep=""), envir=ENVIR)) 
      baseParam <- get(paste("remember",selModelNum,sep=""), envir=ENVIR)

    ## selModelCountPar > 0 hier !!
    for (i in 1:selModelCountPar) {
      baseParam[i] <- 
        if (!is.na(baseParam[i])) baseParam[i] else
        if (modelParam[i,3] == INTEGERLARAM) modelParam[i,2] else
        if (modelParam[i,1] >=0) 0.25 * sum(sqrt(modelParam[i,1:2]))^2 + 0.1
        else 0.5 * (modelParam[i,1] + modelParam[i,2])
 
      #Slider fuer den neuen Parameter 
      slParamValue <- tclVar(baseParam[i])
      entryParamValue <- tclVar(tclvalue(slParamValue))
      # name <- unlist(strsplit(attr(modelParam, "dimnames")[[1]][i],"\\."))[2]
      # slParamName <- tklabel(tt,text=paste(toupper(substring(name, 1,1)), substring(name, 2), sep=""))
      txt <- unlist(strsplit(attr(modelParam, "dimnames")[[1]][i],"\\."))[2]
      slParamName <- tklabel(tt, text=txt)
      slParam <- tkscale(tt, command = plotDensity, from=modelParam[i,1],
                         to=modelParam[i,2],
                         showvalue=FALSE, variable=slParamValue,
                         resolution=if (modelParam[i,3] == INTEGERLARAM) 1 else
                                (modelParam[i,2]-modelParam[i,1])/numberSteps, 
                         orient="horizontal", length=length.slider, width=18)
      entryParam <- tkentry(tt,width=size.entry,textvariable=entryParamValue)
      tkbind(entryParam, "<Return>", OnAddParamEntryChanged)
      
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

    #Print(modelChoice, baseParam, baseModel)

    assign("baseModel", baseModel, ENVIR)
    position()
  }

  OnPlotVarCovChanged <- function(...)
  {
    if((as.character(tclvalue(plotVarCov)) == "Variogram") && !is.null(ev))
      tkconfigure(cbPlotEV, disabled=FALSE)
    else
      tclvalue(plotEV) = "0" 
  #    tkconfigure(cbPlotEV, disabled=TRUE)
    tkrreplot(imgVar)
  }

  OnplotEVChanged <- function(...)
  {
    if((as.character(tclvalue(plotVarCov)) == "Covariance") || is.null(ev))
      tclvalue(plotEV) = "0"
    else {
       tkrreplot(imgVar)
    }
  }

  OnScaleEntryChanged <- function(...)
  { 
    newscale <- log(as.numeric(tclvalue(entryScaleValue)))
    tkconfigure(slScale, to=newscale+abs(newscale),
                from=min(scaleMin,newscale-2),
                resolution=0.01 *
                ((2 * newscale - scaleMin) - min(scaleMin, newscale - 2)))
    tclvalue(slScaleValue) <- log(as.numeric(tclvalue(entryScaleValue)))
  }  

  OnVarEntryChanged <- function(...)
  { 
    newvariance <- log(as.numeric(tclvalue(entryVarianceValue)))
    tkconfigure(slVariance,
                to=newvariance+abs(newvariance),
                from=min(varianceMin,newvariance-2),
                resolution=0.01 *
                (2 * newvariance - varianceMin) -
                min(varianceMin, newvariance - 2))
    tclvalue(slVarianceValue) <- newvariance
  } 

   OnNuggetEntryChanged <- function(...)
  { 
    newnugget <- max(0,as.numeric(tclvalue(entryNuggetValue)))
    tkconfigure(slNugget, to=if(newnugget==0) 1 else 2*newnugget,
                from=0, resolution=0.02*newnugget)
    tclvalue(slNuggetValue) <- newnugget
  } 
#  XX_OnNuggetEntryChanged <- function(...)
#  { 
##    newnugget <- max(0,as.numeric(tclvalue(entryNuggetValue)))
  ##   tkconfigure(slNugget, to=2*newnugget, from=0, resolution=0.02*newnugget)
#    tclvalue(slNuggetValue) <- newnugget
#  } 

  OnRotationEntryChanged <- function(...)
  { 
    tclvalue(slRotationValue) <- as.numeric(tclvalue(entryRotationValue))
  }  

  OnRadiusEntryChanged <- function(...)
  { 
    newscale <- log(as.numeric(tclvalue(entryScaleAValue)))
    tkconfigure(slScaleA, to=newscale+abs(newscale),
                from=min(scaleMin,newscale-2),
                resolution=0.01*
                ((2 * newscale - scaleMin) - min(scaleMin, newscale - 2))) 
    tclvalue(slScaleAValue) <-newscale
    newscale <- log(as.numeric(tclvalue(entryScaleBValue)))
    tkconfigure(slScaleB, to=newscale+abs(newscale),
                from=min(scaleMin,newscale-2),
                resolution=0.01*
                ((2 * newscale - scaleMin) - min(scaleMin, newscale - 2)))
    tclvalue(slScaleBValue) <- newscale
  }  

  OnAddParamEntryChanged <- function(...)
  { 
    baseModel <- get("baseModel", envir=ENVIR)
    if(length(baseModel$k) > 0)
      for (i in 1:length(baseModel$k)) {
        slParamValue <- get(paste("slParam", i, "Value", sep=""), envir=ENVIR)
        value <- get(paste("entryParam", i, "Value", sep=""), envir=ENVIR)
        tclvalue(slParamValue) <- round(as.numeric(tclvalue(value)), digits=2)
      }
    plotDensity()
  } 
    
  GetGuiModel <- function() {
    variance <- exp(as.numeric(tclvalue(slVarianceValue)))
    nugget <- as.numeric(tclvalue(slNuggetValue))

    #Print(baseModel)
    
    baseParam <- baseModel$k
    if(length(baseModel$k) > 0)
      for (i in 1:length(baseModel$k)) { 
        baseParam[i] <- as.numeric(tclvalue(get(paste("slParam", i, "Value",
                                                      sep=""), envir=ENVIR)))
        entryParamValue <-
          get(paste("entryParam", i, "Value", sep=""), envir=ENVIR)
        tclvalue(entryParamValue) <- round(baseParam[i], digits=2)
      }

    baseModel$k <- baseParam
   
    if(!as.numeric(tclvalue(showAniso))) {
      scale <- exp(as.numeric(tclvalue(slScaleValue)))
      newmodel <- list(ZF_SYMBOLS_PLUS,
                    list(DOLLAR[1], var=variance, scale=scale, baseModel),
                    list(DOLLAR[1], var=nugget, list(ZF_NUGGET[1])))
    } else {
      a <-  as.numeric(tclvalue(slRotationValue))
      r <- c(exp(as.numeric(tclvalue(slScaleAValue))),
             exp(as.numeric(tclvalue(slScaleBValue))))
      u <- matrix(c(cos(a), sin(a), -sin(a), cos(a)), ncol=2 )
      aniso <- u %*% (1/r * t(u))
      newmodel <- list(ZF_SYMBOLS_PLUS,
                    list(DOLLAR[1], var=variance, aniso=aniso, baseModel),
                    list(DOLLAR[1], var=nugget, list(ZF_NUGGET[1])))
    }

    #Print(newmodel)
    
    return(newmodel)
  }

  plotFunction <- function(...)
  {

    plotev = as.numeric(tclvalue(plotEV))
    par(cex=0.6, bg="lightgrey", mar=c(3,3,1,1))
    if(!exists("baseModel",envir=ENVIR)) {
      if(!is.null(ev) && plotev) {
        notNA <- !is.nan(ev@emp.vario)
        xm <- c(min(ev@centers[notNA]), max(ev@centers[notNA]))
        ym <- c(min(ev@emp.vario[notNA]), max(ev@emp.vario[notNA])*1.1)
        plot(ev@centers[!is.nan(ev@emp.vario)],
             ev@emp.vario[!is.nan(ev@emp.vario)], pch=19) 
        return(0)
      } 
      plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
      return(0)
    }

    #baseModel <- get("baseModel",envir=ENVIR)
    tclvalue(entryScaleValue) <-
      round(exp(as.numeric(tclvalue(slScaleValue))), digits=2)
    tclvalue(entryVarianceValue) <-
      round(exp(as.numeric(tclvalue(slVarianceValue))), digits=2)
    tclvalue(entryNuggetValue) <-
      round(as.numeric(tclvalue(slNuggetValue)), digits=2)
    tclvalue(entryScaleAValue) <-
      round(exp(as.numeric(tclvalue(slScaleAValue))), digits=2)
    tclvalue(entryScaleBValue) <-
      round(exp(as.numeric(tclvalue(slScaleBValue))), digits=2)
    tclvalue(entryRotationValue) <-
      round(as.numeric(tclvalue(slRotationValue)), digits=2)


    newmodel <- GetGuiModel()
    assign("RFgui.model", newmodel, envir=parent.ev)
    
    if(as.numeric(tclvalue(showAniso))) {
      x1 <- rep(xcov, each=length(ycov))
      x2 <- rep(ycov, times=length(xcov))
      cv <- RFvariogram(x=as.matrix(expand.grid(xcov, ycov)),
                        model=newmodel, 
                        practicalrange = tclvalue(cbPracRangeVal) != "0")
      dim(cv) <- c(length(ycov),length(xcov))

      cv00 <- cv[1,1]
      if (xcov[1] == 0 && ycov[1] == 0) {
        zlim <- c(0, 1.1*max(cv))
        cv[1,1] <- NA
      }
      tranMatrix <- persp(x=xcov, y=ycov, z=cv,
                          theta = as.numeric(tclvalue(slTurnPlotValue)),
                          zlim = zlim, phi = 0, xlab = "x", ylab = "y",
                          zlab = as.character(tclvalue(plotVarCov)),
                          col = "lightblue", ltheta = 120, shade = 0.75,
                          ticktype = "detailed")
      if (xcov[1] == 0 && ycov[1] == 0)
        points(trans3d(xcov[1], ycov[1], cv00, pmat = tranMatrix), pch =16)


    #  Print(newmodel, "hier")
      assign("model", newmodel, envir = ENVIR)
      return(0)
    }
    
    cv <- xcov
    if(as.character(tclvalue(plotVarCov)) == "Covariance") {
      cv <- RFcov(x=xcov, model=newmodel,
                  practicalrange = tclvalue(cbPracRangeVal) != "0")
    }
    if(as.character(tclvalue(plotVarCov)) == "Variogram") {
      pr.dummy <- tclvalue(cbPracRangeVal) != "0"
      ##    save(file="~/xxxx", xcov, newmodel, pr.dummy)
      ##    Print(x=xcov, model=newmodel, practicalrange = pr.dummy )
      cv <- RFvariogram(x=xcov, model=newmodel,
                        practicalrange = pr.dummy)
    }

    xm <- c(min(xcov),max(xcov))
    ym <- c(min(cv), max(cv)*1.1)
    if(!is.null(ev) && plotev) {
      notNA <- !is.nan(ev@emp.vario)
      xm <- c(min(ev@centers[notNA]), max(ev@centers[notNA]))
      ym <- c(min(ev@emp.vario[notNA]), max(ev@emp.vario[notNA])*1.1)
    }
    plot(xcov[2:length(xcov)], cv[2:length(xcov)], type="l",
         xlab="", ylab="", xlim=xm, ylim=ym)
    points(xcov[1], cv[1])

    # plot emp.vario
    if(!is.null(ev) && plotev)    
      points(ev@centers[!is.nan(ev@emp.vario)],
             ev@emp.vario[!is.nan(ev@emp.vario)], pch=19)

    # Print(newmodel, "hierx")
     assign("model", newmodel, envir = ENVIR)
  }

  plotSimulation <- function(...)
  {

    par(cex=0.6, bg="lightgrey", mar=c(3,3,1,1))
    if(!exists("model", envir=ENVIR)) {
      plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=FALSE, xlab="", ylab="")
      return(0)
    }

    simu.model <- get("model", envir=ENVIR)
    if (!is.null(simu.model)) {
      
      ##  assign("cx", x, envir=.GlobalEnv)
      ##  assign("cmodel", model,  envir=.GlobalEnv)
      
      set.seed(fixed.rs)
      par(cex=0.6, bg="lightgrey")
      
      if (guiOpt$simu_method != "any method")
        simu.model <- list(guiOpt$simu_method, simu.model)
      
      ## Print(newmodel)
      yy <- (if (get("simDim", envir = ENVIR) =="sim1Dim") NULL else
             if (is.null(y)) x else y)
      pr <-  tclvalue(cbPracRangeVal) != "0"
      save(file="model", x, simu.model, yy, guiReg, pr)

      z <- try(RFsimulate(x=x, grid=TRUE, model=simu.model,
                          y=if (get("simDim", envir = ENVIR) =="sim1Dim") NULL
                          else if (is.null(y)) x else y,
                          register=guiReg, spConform=TRUE,
                          practicalrange = tclvalue(cbPracRangeVal) != "0"),
               silent=!TRUE)

      if (class(z) == "RFspatialGridDataFrame") plot(z, xlab=NULL, cex=.5)
      else {
        plot(Inf, Inf, xlim=c(0,1), ylim=c(0,1), axes=!FALSE, xlab="", ylab="",
             cex.main=2, col.main="brown",
              main=paste("\n\n\n\n\n\n\n\n",
               if (guiOpt$simu_method == "any method")
               "No suitable simulation algorithm found"
               else paste("Simulation method '", guiOpt$simu_method,
                          "'\ndoes not work for this specification", sep=""),
               ".\n\nSet 'simu_method = \"any method\"'",
               "\nor set 'same.algorithm=TRUE'",
              "\nor see RFoptions() for controlling parameters",
               if (guiOpt$simu_method != "any method")
               paste("\nof '", guiOpt$simu_method, "'"),
               sep=""))
        return(0)
      }
    }
  }

  
  plotDensity <- function(...)
  {
    #tkconfigure(labelOccupancy,textvariable=tclVar("Busy"))
    tkrreplot(imgVar)
    if(as.numeric(tclvalue(simAlways)))
      tkrreplot(imgSim)
    #tkconfigure(labelOccupancy,textvariable=tclVar("Free"))
  }

  OnChangeIsotropie <- function(...)
  {
    if(as.numeric(tclvalue(showAniso))) {
      tkgrid.remove(slScale)
      tkgrid.remove(entryScale)
      tkgrid.remove(labelScale)
    }else {
      tkgrid.remove(slScaleA)
      tkgrid.remove(entryScaleA)
      tkgrid.remove(labelScaleA)
      tkgrid.remove(slScaleB)
      tkgrid.remove(entryScaleB)
      tkgrid.remove(labelScaleB)
      tkgrid.remove(slRotation)
      tkgrid.remove(entryRotation)
      tkgrid.remove(labelRotation)
      tkgrid.remove(slTurnplot)
    }
    position()
    plotDensity()
  }

  OnTurnPlot <- function(...)
  {
    tkrreplot(imgVar)
  }

  OnNewSimu <- function(...) 
  {
    assign("fixed.rs", round(runif(1,1,100000)), envir=ENVIR)  
    tkrreplot(imgSim)
  }

  OnSimDimChanged <- function(...)
  {
    if(as.numeric(tclvalue(showAniso))) 
    {
      tclvalue(rb2DimValue) <-"sim2Dim"
      return(0)
    } 

    if(!sim_only1dim) 
    {
      assign("simDim", tclvalue(rb2DimValue), envir = ENVIR) 
      tkrreplot(imgSim)
    }
  }

  OnReturn <- function(...)
  {
    # hier muss eine rückgabe stehen, emp vario und model mit parametern
   #   Print(GetGuiModel())
    RFoptions(LIST=get("RFopt.old", envir=ENVIR))
    ##remove("RFopt.old", envir=ENVIR)
    assign(".RFgui.exit", TRUE, envir=parent.ev)
    tkdestroy(tt)    
  }

  position <- function(...)
  {  
    #--- DropDown-ComboBox for model selection -------------------------
    tkgrid.configure(labModelSelect, column=col.sl, row=row.sl)
    row.sl <- row.sl+1
    tkgrid.configure(comboBox, column=col.sl, row=row.sl, sticky = "e")
    row.sl <- row.sl+1

    #--- PLOT  ---------------------------------------------------------
   tkgrid.configure(imgVar, rowspan=image.rowspan, columnspan=image.colspan,
                     column=col.var, row=1,sticky="w") 
    tkgrid.configure(imgSim, rowspan=image.rowspan, columnspan=image.colspan,
                     column=col.sim, row=1,sticky="w")

    #--- Radiobutton zur Frage Variogram oder Covarianzfunktion --------
    tkgrid.configure(rbCovariance, column=col.var+image.colspan-1,
                     row=image.rowspan+1, sticky="w")
    tkgrid.configure(labelCovariance, column=col.var+image.colspan-1,
                     row=image.rowspan+1, sticky="e")
    tkgrid.configure(rbVariogram, column=col.var+image.colspan-1,
                     row=image.rowspan+2, sticky="w") 
    tkgrid.configure(labelVariogram, column=col.var+image.colspan-1,
                     row=image.rowspan+2, sticky="e")

    #--- Checkbox show the empirical variogram --------------------------
    tkgrid.configure(cbPlotEV, column=col.var, row=image.rowspan+1, sticky="w")
    tkgrid.configure(labelPlotEV, column=col.var, row=image.rowspan+1,
                     sticky="e")
    
    #--- Radiobutton: select dimension for simulation -------------------
    tkgrid.configure(rbSim1Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+1, sticky="w")
    tkgrid.configure(labelSim1Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+1, sticky="e")
    tkgrid.configure(rbSim2Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+2, sticky="w") 
    tkgrid.configure(labelSim2Dim, column=col.sim+image.colspan-1,
                     row=image.rowspan+2, sticky="e")

    #--- Checkbox simulate on slider movement --------------------------
    tkgrid.configure(cbSimAlways, column=col.sim, row=image.rowspan+1,
                     sticky="w")
    tkgrid.configure(labelSimAlways, column=col.sim, row=image.rowspan+1,
                     sticky="e")
  
    #--- Checkboxes Practical Range  and anisotropy option -------------
    if(!is.null(y)) {
      tkgrid.configure(cbAnisotropy, column=col.sl, row=row.sl, sticky="w")
      tkgrid.configure(labelAniso, column=col.sl, row=row.sl)
      row.sl=row.sl+1
    }
  
    tkgrid.configure(cbPracRange, column=col.sl, row=row.sl, sticky="w")
    tkgrid.configure(labelPracRange, column=col.sl, row=row.sl)
    row.sl=row.sl+1

    #--- Parameterwaehler ----------------------------------------------
   
    if(!as.numeric(tclvalue(showAniso))) {
      #Slider Scale
      tkgrid.configure(labelScale,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkgrid.configure(slScale,  column=col.sl, row=row.sl, sticky="w")
      tkgrid.configure(entryScale, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1
    }else {
      #Slider Rotation
      tkgrid.configure(labelRotation,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkgrid.configure(slRotation,  column=col.sl, row=row.sl, sticky="w")
      tkgrid.configure(entryRotation, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1

      #Slider Radius
      tkgrid.configure(labelScaleA,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkgrid.configure(slScaleA,  column=col.sl, row=row.sl, sticky="w")
      tkgrid.configure(entryScaleA, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1

      tkgrid.configure(labelScaleB,  column=col.sl, row=row.sl)
      row.sl <- row.sl+1
      tkgrid.configure(slScaleB,  column=col.sl, row=row.sl, sticky="w")
      tkgrid.configure(entryScaleB, column=col.sl, row=row.sl, sticky="e")
      row.sl <- row.sl+1

      #Slider turn the now 2dim covarianz plot
      tkgrid.configure(slTurnplot, column=col.var,columnspan=image.colspan,
                       row=image.rowspan+3)
    }

    #Slider variance
    tkgrid.configure(labelVariance,  column=col.sl, row=row.sl)
    row.sl <- row.sl+1
    tkgrid.configure(slVariance,  column=col.sl, row=row.sl, sticky="w")
    tkgrid.configure(entryVariance, column=col.sl, row=row.sl, sticky="e")
    row.sl <- row.sl+1

    #Slider nugget
    tkgrid.configure(labelNugget,  column=col.sl, row=row.sl)
    row.sl <- row.sl+1
    tkgrid.configure(slNugget,  column=col.sl, row=row.sl, sticky="w")
    tkgrid.configure(entryNugget, column=col.sl, row=row.sl, sticky="e")
    row.sl <- row.sl+1

    if(exists("baseModel", envir=ENVIR)) {
      baseModel <- get("baseModel",envir=ENVIR)
      baseParam <- baseModel$k
      if(length(baseModel$k) > 0)
        for (i in 1:length(baseModel$k)) { 
          tkgrid.configure(get(paste("slParam",i,"Name",sep=""),
                               get("slParam", envir=ENVIR),
                               envir=ENVIR), column=col.sl, row=row.sl)
          row.sl <- row.sl+1 
          tkgrid.configure(get(paste("slParam", i, sep=""),
                               get("slParam", envir=ENVIR),
                               envir=ENVIR),
                           column=col.sl, row=row.sl, sticky="w")  
          tkgrid.configure(get(paste("entryParam", i, sep=""),
                               get("entryParam", envir=ENVIR),
                               envir=ENVIR),
                           column=col.sl, row=row.sl, sticky="e")
          row.sl <- row.sl+1
        }
    }

    #--- Buttons - new simulation (new seed), return ---------------
    tkgrid.configure(buttonNewSimu, column=col.sl,
                     row=max(row.sl,image.rowspan+1), sticky="e")
    row.sl=row.sl+1
    tkgrid.configure(buttonReturn, column=col.sl,
                     row=max(row.sl,image.rowspan+2), sticky="e")
    row.sl=row.sl+1
    #--- Beschäftigungsindikator ---------------------------------------
#    tkgrid.configure(labelOccupancy, row=row.last, column=col.sl, sticky="e")
  } ## end fct position

  PRINT <- FALSE

  if (missing(x)) {
    stopifnot(is.null(y))    
    if (missing(data)) {
      x <- seq(1, 5, len=guiOpt$size[if (sim_only1dim) 1 else 2] )
    } else {      
      stopifnot(is(data, "RFsp"))
      if (!is.null(y)) stop("'x' and 'y' may not be given if 'data' is given.")

      str(data)
      
      r <- apply(data@coords, 2, range)
      len <- guiOpt$size[if (length(r) == 2) 1 else 2] 
      x <- seq(r[1], r[2], len=len)
      if (length(r) > 2) {
         y <- seq(r[3], r[4], len=len)
      }
    }
  }

  if(!missing(data) && !is.null(data)) {
    if (!is.null(ev)) stop("if 'data' is given, 'ev' may not be given.")
    ev <- RFempiricalvariogram(data=data, phi=1)
    ## RFfit(data) ???
  }
  
  if (any(diff(x) <= 0)) 
    stop("x should be a sequence of increasing numbers")

  if (!is.null(y) && any(diff(y) <= 0))
    stop("y should be a sequence of increasing numbers")
  
  if(is.null(xcov)) {
    if(is.null(ev))
      xcov <- seq(0,15,0.1)
    else
      xcov <-seq(min(0,0.9*min(ev@centers)), max(ev@centers),
                 by=diff(range(ev@centers))/100)
  }
  
  if(is.null(ycov) && !is.null(y)) {
    if(is.null(ev))
      ycov <- seq(0,15,by=0.1)
    else
      ycov <-seq(min(0, 0.9*min(ev@centers)), max(ev@centers),
                 by=diff(range(ev@centers))/100)
  }


  if (exists("baseModel", where=ENVIR)) remove("baseModel", envir=ENVIR)

  ##  Print("removing model")
  ##  if (exists("model", where=ENVIR)) remove("model", envir=ENVIR)

             
  if(is.null(fixed.rs))
    fixed.rs <- ifelse(exists(".Random.seed"), .Random.seed, runif(1,0,100000)) 


  # get all model names
  models <-
    RFgetModelNames(type=ZF_TYPE[c(TcfType, PosDefType, UndefinedType) + 1],
                    domain=ZF_DOMAIN[TRANS_INV + 1],
                    isotropy=ZF_ISOTROPY[ISOTROPIC + 1],
                    operator=FALSE, 
                    vdim=1)#, multivariate = 1)  
  
  #-------------------------------------------------------------------
  # Start Values and ranges
  #-------------------------------------------------------------------
  cbPracRangeVal <- tclVar(RFopt$general$practicalrange)
  simAlways <- tclVar(as.integer(guiOpt$alwaysSimulate))
  plotVarCov <- tclVar("Variogram")
  plotEV <- tclVar(ifelse(is.null(ev) && tclvalue(plotVarCov)=="Variogram",
                          "0", "1"))
  showAniso <- tclVar("0")
  slTurnPlotValue <- tclVar("0")
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
  slNuggetValue <- tclVar(nugget)
   ## Variance
  varianceMin <- round(log(0.01), digits=2)
  varianceMax <- log(max(1e-10, nuggetMax))
  slVarianceValue <- tclVar(log(variance))
  ## Scale
  scaleMin <- round(log(0.1*scale), digits=2)  
  scaleMax <- round(log(10*scale), digits=2)   
  slScaleValue <- tclVar(log(scale))
  
  ## die direkte eingabe muss als variable getrennt von den schiebern laufen
  entryScaleValue <- tclVar(scale)
  entryVarianceValue <- tclVar(variance)
  entryNuggetValue <- tclVar(nugget) 
  ## die direkte eingabe muss als variable getrennt von den schiebern laufen
  entryScaleValue <- tclVar(scale)
  entryVarianceValue <- tclVar(variance)
  entryNuggetValue <- tclVar(nugget) 

  # die anisotropie tierchen
  slRotationValue <- tclVar("0")
  entryRotationValue <- tclVar("0")
  anisoScale = "1"
  slScaleAValue <- tclVar(anisoScale)
  slScaleBValue <- tclVar(anisoScale)
  entryScaleAValue <- tclVar(anisoScale)
  entryScaleBValue <- tclVar(anisoScale)
  radiusMax <- 2

  #------------------------------------------------------------------
  # GUI
  #------------------------------------------------------------------
  tt <- tktoplevel()#title="RFgui")
  tktitle(tt) <- "RFgui"
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

  tkgrid(tklabel(tt, text="", width=1, height=0), column=0, row=0)
  tkgrid(tklabel(tt, text="", width=1), column=col.sim+image.colspan, row=1)
  tkgrid(tklabel(tt, text="", width=1), column=col.var+image.colspan, row=1)
  tkgrid(tklabel(tt, text="", width=1), column=col.sl+image.colspan, row=1)

  #--- DropDown-ComboBox for model selection -------------------------
  labModelSelect <- tklabel(tt,text="Adjust Model:")
  textModell <- tclVar("Please select a model...")
  #if(!is.null(model))
  #  textModel <- 
  comboBox <- ttkcombobox(tt,textvariable=textModell, state="readonly",
                          values=models)
  tkbind(comboBox, "<<ComboboxSelected>>", OnModelSelected)

  #--- PLOT  ---------------------------------------------------------
  imgVar <- tkrplot(tt,fun=plotFunction,hscale=plothscale,vscale=plotvscale)
  imgSim <- tkrplot(tt,fun=plotSimulation,hscale=plothscale,vscale=plotvscale)
    
  #--- Beschäftigungsindikator -------------------------------------
  labelOccText <- tclVar("Free")
  labelOccupancy <- tklabel(tt,text=tclvalue(labelOccText))

  #--- anistropy -----------------------------------------------------
  cbAnisotropy <- tkcheckbutton(tt, variable=showAniso,
                                command=OnChangeIsotropie)
  labelAniso <-  tklabel(tt,text="Anisotropy")

  #--- Radiobutton zur Frage Variogram oder Covarianzfunktion --------
  rbVariogram <- tkradiobutton(tt, command=OnPlotVarCovChanged)
  rbCovariance <- tkradiobutton(tt, command=OnPlotVarCovChanged)
  tkconfigure(rbVariogram,variable=plotVarCov,value="Variogram")
  tkconfigure(rbCovariance,variable=plotVarCov,value="Covariance")
  labelCovariance <- tklabel(tt,text="Covariance")
  labelVariogram <- tklabel(tt,text="Variogram")

  #--- Checkbox plot empirical variogram    --------------------------
  cbPlotEV <- tkcheckbutton(tt, variable=plotEV, command=OnplotEVChanged)
  labelPlotEV <- tklabel(tt,text="Plot empirical variogram")

  #--- Radiobutton: select dimension for simulation ------------------
  rbSim1Dim <- tkradiobutton(tt, command=OnSimDimChanged)
  rbSim2Dim <- tkradiobutton(tt, command=OnSimDimChanged)
  rb2DimValue <- tclVar(ifelse((sim_only1dim), "sim1Dim", "sim2Dim"))
  tkconfigure(rbSim1Dim,variable=rb2DimValue, value="sim1Dim")
  tkconfigure(rbSim2Dim,variable=rb2DimValue, value="sim2Dim")
  labelSim1Dim <- tklabel(tt,text="1 dim")
  labelSim2Dim <- tklabel(tt,text="2 dim")
  assign("simDim", as.character(tclvalue(plotVarCov)), envir=ENVIR)

  #--- Button - new simulation (new seed) ----------------------------
  buttonNewSimu <- tkbutton(tt,text="New Simulation",command=OnNewSimu)
  
  #--- Button - Return -----------------------------------------------
  buttonReturn <- tkbutton(tt,text="      Return       ",command=OnReturn)

  #--- Checkbox simulate on slider movement --------------------------
  cbSimAlways <- tkcheckbutton(tt, variable=simAlways)
  labelSimAlways <- tklabel(tt,text="Simulate always")

  #--- Checkbox Practical Range --------------------------------------
  cbPracRange <- tkcheckbutton(tt, command=plotDensity)
  tkconfigure(cbPracRange,variable=cbPracRangeVal)
  labelPracRange <- tklabel(tt,text="Practical Range")

  #--- Slider turn covariance plot in anisotropic case ----------------
  slTurnplot <- tkscale(tt, command = OnTurnPlot, from=0, to=360,
                        showvalue=TRUE, variable=slTurnPlotValue,
                        resolution=360/numberSteps, orient="horizontal",
                        length=2*length.slider, width=18)
 # labelRotation <- tklabel(tt,text="Rotation")

  #--- Parameterwaehler ----------------------------------------------
  # Aniso Rotation
  slRotation <- tkscale(tt, command = plotDensity, from=0, to=0.5*pi,
                        showvalue=FALSE, variable=slRotationValue,
                        resolution=0.5*pi/numberSteps, orient="horizontal",
                        length=length.slider, width=18)
  labelRotation <- tklabel(tt,text="Rotation")
  entryRotation <- tkentry(tt,width=size.entry, textvariable=entryRotationValue)
  tkbind(entryRotation, "<Return>", OnRotationEntryChanged)

  # Aniso Scales in first and second axis
  slScaleA <- tkscale(tt, command = plotDensity, from=scaleMin, to=scaleMax,
                      showvalue=FALSE, variable=slScaleAValue,
                      resolution=radiusMax/numberSteps, orient="horizontal",
                      length=length.slider, width=18)
  labelScaleA <- tklabel(tt,text="first axis scale")
  entryScaleA <- tkentry(tt,width=size.entry, textvariable=entryScaleAValue)
  tkbind(entryScaleA, "<Return>", OnRadiusEntryChanged)

  slScaleB <- tkscale(tt, command = plotDensity, from=scaleMin, to=scaleMax,
                      showvalue=FALSE, variable=slScaleBValue,
                      resolution=radiusMax/numberSteps, orient="horizontal",
                      length=length.slider, width=18)
  labelScaleB <- tklabel(tt,text="second axis scale")
  entryScaleB <- tkentry(tt,width=size.entry, textvariable=entryScaleBValue)
  tkbind(entryScaleB, "<Return>", OnRadiusEntryChanged)

  # Scale
  slScale <- tkscale(tt, command = plotDensity, from=scaleMin, to=scaleMax,
                     showvalue=FALSE, variable=slScaleValue,
                     resolution=(scaleMax-scaleMin)/numberSteps,
                     orient="horizontal", length=length.slider, width=18)
  labelScale <- tklabel(tt,text="Scale")
  entryScale <- tkentry(tt,width=size.entry, textvariable=entryScaleValue)
  tkbind(entryScale, "<Return>", OnScaleEntryChanged)

  #Slider variance
  slVariance <- tkscale(tt, command = plotDensity, from=varianceMin,
                        to=varianceMax,
                        showvalue=FALSE, variable=slVarianceValue,
                        resolution=(varianceMax-varianceMin)/numberSteps,
                        orient="horizontal", length=length.slider, width=18)
  labelVariance <- tklabel(tt,text="Variance")
  entryVariance <- tkentry(tt,width=size.entry, textvariable=entryVarianceValue)
  tkbind(entryVariance, "<Return>", OnVarEntryChanged)

  #Slider nugget
  slNugget <- tkscale(tt, command = plotDensity, from=nuggetMin, to=nuggetMax,
                      showvalue=FALSE, variable=slNuggetValue,
                      resolution=(nuggetMax-nuggetMin)/numberSteps,
                      orient="horizontal", length=length.slider, width=18)
  labelNugget <- tklabel(tt,text="Nugget")
  entryNugget <- tkentry(tt,width=size.entry, textvariable=entryNuggetValue)
  tkbind(entryNugget, "<Return>", OnNuggetEntryChanged)

  position()
}
