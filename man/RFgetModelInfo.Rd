\name{RFgetModelInfo}
\alias{RFgetModelInfo}
\alias{RFgetModelInfo_register}
\alias{RFgetModelInfo_model}
%\alias{GetModelInfo}
%\alias{GetModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Information on RMmodels}
\description{
  The function returns information about a \link{RMmodel},
  either internal information when used in simulations, for instance,
  or general information
}
\usage{
RFgetModelInfo(...) 

RFgetModelInfo_register(register, level = 1, spConform =
                 RFoptions()$general$spConform, which.submodels =
                 c("user", "internal", "call+user", "call+internal",
                 "user.but.once", "internal.but.once",
                 "user.but.once+jump", "internal.but.once+jump",
                 "all"), modelname = NULL)

RFgetModelInfo_model(model, dim = 1, Time = FALSE,
                     kernel = FALSE, exclude_trend = TRUE, ...)
}

\arguments{
  \item{...}{any of the arguments below}
  \item{register}{\eqn{0,...,21} or an evaluating function,
    e.g. \command{\link{RFsimulate}}.
    Place where intermediate calculations are stored.
    See also section \code{Registers} in \command{\link{RFoptions}}.
  }
  \item{level}{integer [0..5]; level of details, i.e. the higher the
    number the more details are given. 
  }

  % \item{max.elements}{integer; since \command{\link{RFgetModelInfo}}
  % might be a dump of simulation that needs a huge amount of
  % memory and since all entries are copied, the maximal amount
  % of available memory might be easily exceeded. Therefore,
  % only the size of the critical parts are returned and not the
  % vector or matrix itself, if the number of elements exceeds
  % \code{max.elements}. 
% }
\item{spConform}{see \command{\link{RFoptions}}
}
\item{which.submodels}{
  Internally, the sub-models are represented in two different ways
  \sQuote{internal} and \sQuote{user}. The latter is very close to
  the model defined by the user.
   
  Most models have a leading internal model.
  The values \code{"call+user"}, \code{"call+internal"} return also
  this leading model if existent.
  
  The values \code{"user.but.once"}, \code{"internal.but.once"}
  \code{"user.but.once"} returns the user path of the internal model
  following the leading model.  \code{"internal.but.once"} would return the
  internal path of the user model following the leading model, but this
  path should never exist. So as all the other options if a certain
  direction does not exist, the alternative path is taken.
  
  The values \code{"user.but.once+jump"},
  \code{"internal.but.once+jump"}
  same as \code{"user.but.once"} and \code{"internal.but.once"},
  except that the first submodel below the leading model is not given.

  The value \code{"all"} returns the whole tree of models (very advanced).
}
\item{modelname}{ string. If \code{modelname} is given then
  it returns the first appearance of the
  covariance model with name \code{modelname}.
  If \code{meth} is given then the model within the method is returned.
}
\item{model}{an \command{\link{RMmodel} with \code{NA}s where
    information on the parameters is requested}.
}
\item{dim}{positive integer. Spatial dimension.
}
\item{Time}{logical. Should time be considered too?
}
\item{kernel}{logical. Should the model be considered as a kernel?
}
\item{exclude_trend}{logical. Currently, only \code{TRUE} is available.
}
}
\details{
  \command{RFgetModelInfo} branches either into
  \command{RFgetModelInfo_register} or \command{RFgetModelInfo_model},
  depending on the type of the \emph{first} argument. The latter two are
  usually not called by the user.

  \command{RFgetModelInfo} has three standard usages:

  \itemize{
    \item \code{RFgetModelInfo()} returns internal information on the
    last call of an \code{\link{RF}} function.
    \item \code{RFgetModelInfo(RFfunction)} returns internal information on the
    last call of \code{\link{RF}function}.
    \item \code{RFgetModelInfo(RMmodel)} returns general information on
    \code{\link{RMmodel}}
  }
    
  Whereas \command{RFgetModelInfo()} 
  can return detailed internal information,
  \command{\link{RFgetModel}} returns a model that can be re-used
  by the user.
}
\note{
 Put \code{Storing=TRUE}, see \command{\link{RFoptions}}
 if you like to have more
 internal information in case of failure of an initialisation of
 a random field simulation.
}
\value{
  If \code{RFgetModelInfo(model)} is called a list returned with the
  following elements:
  \itemize{
    \item \code{trans.inv} : logical. Whether the model is translation
    invariant (stationary)
    \item \code{isotropic} : logical. Whether the model is rotation
    invariant (stationary)
    \item \code{NAs} : is case of an additive model it gives the number
    of NAs in each submodel
    \item \code{minmax} : a data frame containing information on all
    arguments set to \code{NA}s
    \itemize{
      \item \code{pmin}, \code{pmax} : lower and upper endpoint of the
      parameter values usually found in practice
      \item \code{type} : integer; recognized particularities of a parameter;
      an explanation of the values is given after the table, if printed.
      \item \code{NAN} : the number of \code{NAN}s found
      \item \code{min}, \code{max} : mathetically valid lower and upper
      endpoints of the parameter values
      \item \code{omin}, \code{omax} : logical. If \code{FALSE} the
      respective mathematical endpoint is included
      \item \code{col}, \code{row} : the dimension of the parameter.
      If the parameter is a scalar then \code{col = row = 1}. If it is a
      vector then \code{col = 1}.
      \item \code{bayes} : currently not used (always \code{FALSE})
    }
  }
  
  Else a list of internal structure is returned.
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\seealso{command{\link{RFgetModel}}, \command{\link{RFsimulate}}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
model <- RMexp(scale=4, var=2) + RMnugget(var=3) + RMtrend(mean=1)
z <- RFsimulate(model, 1:4, storing=TRUE)
RFgetModelInfo()

model <-  RMwhittle(scale=NA, var=NA, nu=NA) + RMnugget(var=NA)
RFgetModelInfo(model)
\dontshow{FinalizeExample()}

}
\keyword{spatial}
