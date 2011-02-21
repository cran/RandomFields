

# options(warn=0);library(RandomFields);source("~/R/RF/RandomFields/R/getNset.R");source("~/R/RF/RandomFields/R/rf.R");source("RandomFields/R/ShowModels.R"); ShowModels(x=1:100, y=1:100, model=list(list(model="exp", var=1, aniso=c(1,0,0,1))), method="ci");



ShowModels <- function(x, y=NULL,
                       covx=ifelse(is.null(empirical),diff(range(x))/5,
                         max(empirical$c)),
                       fixed.rs=TRUE,
                       method=NULL,
                       empirical=NULL,
                       model=NULL,
                       param=NULL,
                       anisotropy=FALSE,
                       all.param=NULL,## var, nugg, scale
                       legends = TRUE,
                       register=0,
                       Mean=NULL,
                       erase=TRUE,
                       x.fraction=0.60,
                       cex.names=1,
                       covx.default = 100,
                       link.fct=NULL,
                       Zlim=NULL,
                       Col.rect="red", Col.bg="blue", Col.sep="grey",
                       Col.line="red", Col.txt="black", Col.flash="red",
                       Col.vario="blue", Col.main="black",
                       Col.model=c("red", "black"), ## rf, transformed,
                       ##                                   vario
                       vario.lty=c(1,2), ## main axis, sec. axis
                       cex.leg =0.7 * cex.names,
                       cex.eval=0.8 * cex.names,
                       update=TRUE,
                       screen.new=TRUE,
                       use.outer.RFparameters=FALSE,
                       debug=FALSE,
                       ...){
#  stop("unchecked")

  ## here: only non-operator models and  kappasmodel == kappasize1,
  ## getNset.cc !

  cat("I am afraid. ShowModels is currently not available.\nPlease use RandomFields Version 1.3.35, instead.\n")
  return(invisible(NULL))
} 
