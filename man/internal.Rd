\name{Low level functions}
\alias{rfGenerateModels}
\alias{checkExamples}
\alias{plotWithCircles}
\title{Low level functions}
\description{
 These functions are internal and should not be used.
}
\usage{
rfGenerateModels(assigning,
 RFpath = "~/R/RF/svn/randomfields_2",
 RMmodels.file = paste(RFpath, "R/RMmodels.R", sep="/")
 )
 
checkExamples(exclude = NULL, include=1:length(.fct.list),
               ask=FALSE, echo=ask, halt=FALSE, ignore.all = FALSE,
               path="randomfields_2", package = "RandomFields",
                          read.rd.files=TRUE)

plotWithCircles(data, factor=1.0,
                xlim=range(data[,1])+c(-maxr,maxr),
                ylim=range(data[,2])+c(-maxr,maxr),
                col=1, fill=0, ...)
}
\arguments{
  \item{assigning, RFpath, RMmodels.file}{internal}
  \item{exclude, include, ask, echo, halt, ignore.all, path, package,
    read.rd.files}{internal; 
    ignore.all referes to the \sQuote{all} export statement in the
    namespace -- whether this should be ignored.
    if \code{read.rf.files} is \code{TRUE} or a path to the Rd files, then
    the man pages are analysed to get all examples; \code{ignore.all} is
    then ignored. If \code{FALSE} only examples of functions (which are
                search in the environments) are run.
  `}
  \item{data, factor, xlim, ylim, col, fill, ...}{internal}
 % \item{REGISTER, COVREGISTER, RELAX}{internal}
}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}
}
\keyword{spatial}

\examples{
set.seed(0)
\dontshow{\dontrun{
## check own examples
checkExamples(ignore.all=TRUE, halt=!TRUE)


### Check the examples of the other packages:
dep.packages <- c( #"DSpat","lgcp", "constrainedKriging", "MarkedPointProcess",
                  "Geneland", "glmmBUGS", 
                  "ProbForecastGOP","geoR",
                  "CompRandFld", ## RFsim does not work in version 2.1.18
                  "fractaldim", "rpanel", "spatstat")
#for (i in dep.packages) library(i, character.only=TRUE)
#for (i in dep.packages) install.packages(i)


not.working <- list()
for (.i in 1:length(dep.packages)) {
  not.working[[.i]] <-
    checkExamples(path=paste("~/TMP/dep.packages", dep.packages[.i], sep="/"),
                   package=dep.packages[.i])
  Print(.i, not.working); repeat{ if (readline()=="e") break}
}

Print(not.working)
}}

# for (i in dep.packages) cat(maintainer(i), "\n") 
}
