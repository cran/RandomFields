\name{jss13}
\alias{jss13}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{\R package \pkg{RandomFields}: Analysis and simulation of
 multivariate random fields and more}
\description{
 Here the code of the paper on \sQuote{\R package \pkg{RandomFields}:
 Analysis and simulation of 
 multivariate random fields and more} is given. 
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de};
}
\references{
 \itemize{ 
 \item
 Schlather, M., Malinowski, A., Menck, P.J., Oesting, M. and Strokorb, K. (2013)
 \R package \pkg{RandomFields}: Analysis and simulation of
 multivariate random fields and more. \emph{Submitted}
 }
}
 

\examples{
set.seed(0)

# the full code within the paper is the following -- see also examples{weather}

\dontshow{\dontrun{
## Fig. 1
model <- RMmatrix(M = c(0.9, 0.43), RMwhittle(nu = 0.3)) + 
 RMmatrix(M = c(0.6, 0.8), RMwhittle(nu = 2))
x <- y <- seq(-10, 10, 0.2)
simu <- RFsimulate(model, x, y, grid=TRUE)
dev.new(height=5, width=8)
plot(simu)
dev.copy2pdf(file="pdf/lmc.pdf")


## Fig. 2
model <- RMdelay(RMstable(alpha=1.9, scale=2), s=c(4, 4))
simu <- RFsimulate(model, x, y, grid=TRUE)
dev.new(height=5, width=8)
plot(simu, zlim='joint')
dev.copy2pdf(file="pdf/delay0.pdf")


## Fig. 3
model <- RMdelay(RMstable(alpha=1.9, scale=2), s=c(0, 4)) + 
 RMdelay(RMstable(alpha=1.9, scale=2), s=c(4, 0))
simu <- RFsimulate(model, x, y, grid=TRUE)
dev.new(height=5, width=8)
plot(simu, zlim='joint')
dev.copy2pdf(file="pdf/delay.pdf")
dev.new(height=9, width=8)
plot(model, dim=2, xlim=c(-5, 5), main="Covariance function", 
 cex=1.25, labcex=1, xlab=" ")
dev.copy2pdf(file="pdf/delay_cov.pdf")


## Fig. 4
model <- RMgencauchy(alpha=1.5, beta=3)
simu <- RFsimulate(model, x, y, z=c(0,1), grid=TRUE)
dev.new(height=5, width=8)
plot(simu, MARGIN.slices=3, n.slices=2)
dev.copy2pdf(file="pdf/latent.pdf")


## Fig. 5
model <- RMbiwm(nudiag=c(1, 2), nured=1, rhored=1, cdiag=c(1, 5), 
 s=c(1, 1, 2))
simu <- RFsimulate(model, x, y, grid=TRUE)
dev.new(height=5, width=8)
plot(simu)
dev.copy2pdf(file="pdf/biwm.pdf")


## Fig. 6 & 7
model <- RMcurlfree(RMmatern(nu=5), scale=4)
simu <- RFsimulate(model, x, y, grid=TRUE)
dev.new(height=5, width=11)
plot(simu, select=list(1, 1:3, 4))
dev.copy2pdf(file="pdf/divfree.pdf")
dev.new(height=8.5, width=8)
plot(model, dim=2, xlim=c(-3, 3), main="")
dev.copy2pdf(file="pdf/divfree_cov.pdf")


## Fig. 8
x <- y <- seq(-2, 2, len=20)
model <- RMkolmogorov()
simu <- RFsimulate(model, x, y, z=0, grid=TRUE)
dev.new(height=5, width=4.5)
plot(simu, select=list(1:2), col=c("red"))
dev.copy2pdf(file="pdf/kolmogorov.pdf")
dev.new(height=6.5, width=6)
plot(model, dim=3, xlim=c(-3, 3), MARGIN=1:2, fixed.MARGIN=1, main="")
dev.copy2pdf(file="pdf/kolmogorov_cov.pdf")


## Section 5 : Inference

#################################
## output ##
#################################

PaperOutput <- function(m, sdP, sdT) {
 ## correct for miles2l=km and for sdP and sdT
 if (pars <- !is.null(m$C1$scale)) { ## ! parsimonious
 m$C1$scale <- m$C1$scale * miles2km 
 m$C1$phi$cdiag <- m$C1$phi$cdiag * c(sdP, sdT)^2
 m$C1$phi$c <- m$C1$phi$c * c(sdP^2, sdP * sdT, sdT^2)
 biwm <- m$C1$phi
 } else {
 m$C1$cdiag <- m$C1$cdiag * c(sdP, sdT)^2
 m$C1$c <- m$C1$c * c(sdP^2, sdP * sdT, sdT^2)
 m$C1$s <- m$C1$s * miles2km
 biwm <- m$C1
 }
 m$C0$M <- m$C0$M * c(sdP, 0, 0, sdT)
 
 ml <- RFfit(distances=Dist.mat * miles2km,
             dim = 2, data=t(t(PT) * c(sdP, sdT)),
             model=m, grid=FALSE, meth="ml",
             spC=FALSE)$ml$ml
 
 sigmaP <- sqrt(biwm$c[1])
 sigmaT <- sqrt(biwm$c[3])
 
 return(list(#model = m,
 sigmaP = sigmaP,
 sigmaT = sigmaT,
 nuP = biwm$nu[1],
 nuT = biwm$nu[3],
 nuPT = if(biwm$c[2]==0) NA else biwm$nu[2],
 inv.aP = if (pars) m$C1$s else biwm$s[1],
 inv.aT = if (pars) m$C1$s else biwm$s[3] ,
 inv.aPT= if (pars) m$C1$s else if(biwm$c[2]==0) NA else biwm$s[2],
 rhoPT = biwm$c[2] / (sigmaP * sigmaT), 
 tauP = m$C0$M[1],
 tauT = m$C0$M[4],
 ml = ml
 ))
}


#################################
## get the data ##
#################################
library(fields)
miles2km <- 1.608

data(weather)
sdP <- sd(weather[, 1])
sdT <- sd(weather[, 2])
PT <- cbind( weather[, 1] / sdP, weather[, 2] / sdT )

Dist.mat <- rdist.earth(weather[, 3:4])
Dist.mat <- Dist.mat[lower.tri(Dist.mat)] ## in miles

 
nug <- RMmatrix(M=matrix(nc=2, c(NA, 0, 0, NA)), RMnugget())

modes <- c("sloppy", "easy", "normal")

RFoptions(spConform = FALSE)

#################################
## MLE ##
#################################
for (mode in modes) { ## takes about 2 hours !!
 cat("mode\n")
 file <- paste(mode, ".rda", sep="")
 #% if (FileExists(file)) next;

 #################################
 ## first: independent model ##
 #################################
 indep.model <- nug + RMbiwm(nudiag=c(NA,NA), s=c(NA, 1, NA), c=c(NA, 0, NA))
 indep <- RFfit(indep.model, distances=Dist.mat, dim=2, data=PT,
                meth="ml", modus_op=mode)
 
 #################################
 ## second: parsimoninous model ## 
 #################################
 pars.model <- nug + RMbiwm(nudiag=c(NA,NA), scale=NA, c=rep(NA, 3))
 pars <- RFfit(pars.model, distances=Dist.mat, dim = 2, data=PT,
               meth="ml", modus_op=mode)

 #################################
 ## third: full biwm model ##
 #################################
 full.model <- nug + RMbiwm(nu=rep(NA, 3), s=rep(NA, 3), c=rep(NA, 3))
 full <- RFfit(full.model, distances=Dist.mat, dim = 2, data=PT,
               meth="ml", modus_op=mode)


 RFoptions(spConform = TRUE)

 #################################
 ## output ##
 #################################

 i <- PaperOutput(indep$ml$model, sdP, sdT)#% ever best: -1276.75
 p <- PaperOutput(pars$ml$model, sdP, sdT) #% ever best: -1265.73
 f <- PaperOutput(full$ml$model, sdP, sdT) #% ever best: -1265.30
 table <- rbind(indep=unlist(i), pars=unlist(p), full=unlist(f))
 print(table, digits=3)
 print(table[, ncol(table)], digits=8)

 save(file=file, indep, pars,full, full.model, pars.model, indep.model,
 sdP, sdT, table, modes)


 file <- paste(mode, ".tex", sep="")
 model.names <- c("indep.{}", "pars.{}", "full")
 write(file=file,
 paste("\\begin{table}[p]\\small\n\\centering\n\\caption{Comparison ",
 "between the bivariate Whittle-Mat\'ern models 'independent',",
 " 'parsimonious' and 'full'. Here $\\sigma_P=\\sqrt{c_{pp}}$ and",
 " $\\nu_P$ instead of $\\nu_{PP}$ etc., for short.", 
 " Optimization mode is '", mode, "'. Here, \\label{tab:", mode,
 "}}\n\\begin{tabular}{l",
 paste(rep("c", ncol(table)), collapse=""), "}\n",
 " & $\\sigma_P$ & $\\sigma_T$ & $\\nu_P$ & $\\nu_T$ & ",
 "$\\nu_{PT}$ & $s_P$ & $s_T$ & $s_{PT}$ & $\\rho_{PT}$ &",
 " $\\tau_P$ & $\\tau_T$ & ML \\\\\\hline",
 sep="")
 )
 for (nr in 1:length(modes)) {
 write(file=file, append=TRUE,
 paste(model.names[nr], " & ",
 paste(signif(table[nr, (0:1) - ncol(table)], 3),
 collapse=" & "),
 " & ", signif(table[nr, ncol(table)-1], 1),
 " & ", signif(table[nr, ncol(table)], 7), "\\\\"))
 }
 write(file=file, append=TRUE, "\\end{tabular}\n\\end{table}\n", sep="")
}

 }
}
}

\seealso{
 \link{weather}, \link{SS12}, \link{S10}
}

\keyword{spatial}

