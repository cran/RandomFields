\name{RFsimulate}
\alias{RFsimulate}
\title{Simulation of Random Fields}

\description{
  This function simulates \bold{unconditional} random fields:
  \itemize{
    \item univariate and multivariate,
    spatial and spatio-temporal \link[=Gaussian]{Gaussian random fields}
    \item fields based on Gaussian fields such as \link[=RPchi2]{Chi2 fields}
    or \link[=RPbernoulli]{Binary fields}, see \link{RP}.
    \item \link[=RPpoisson]{stationary Poisson fields}
    \item \link[=RPmaxstable]{stationary max-stable random fields}.
  }

  It also simulates \bold{conditional} random fields for
  \itemize{
    \item univariate and multivariate,
    spatial and spatio-temporal Gaussian random fields
  }
  
  Here, only the simulation of Gaussian random fields is described.
  For other kinds of random fields (binary, max-stable, etc.) or
  more sophisticated approaches see \link{RFsimulateAdvanced}.
}

\usage{
RFsimulate(model, x, y=NULL, z=NULL, T=NULL, grid=NULL,
           distances, dim, data, given=NULL, err.model, params,
           err.params, n=1, ...)
}

\arguments{
 \item{model,params}{\argModel
%   \itemize{
%     \item if of class \code{\link[=RMmodel-class]{RMmodel}}, \code{model}
%     specifies a  covariance or variogram model of a Gaussian random     field;
%     type \code{\link{RFgetModelNames}(type="variogram")} for a list of
%     available models; see also \command{\link{RMmodel}}.
%     
%     \item if of class \code{\link{RFformula}} or
%     \code{\link[methods:formula-class]{formula}} ,
%     \code{submodel} specifies a linear mixed model where random
%     effects can be modelled by Gaussian random fields;
%     see \command{\link{RFformula}} for details on model
%     specification.
%     
%     \item for (many) more options see \link{RFsimulateAdvanced}.
%   }
 }
 \item{x}{\argX}
 \item{y,z}{\argYz}
 \item{T}{\argT}
 \item{grid}{\argGrid}
 \item{distances,dim}{\argDistances}
 \item{data}{For conditional simulation and random imputing only. If
   \code{data} is missing, unconditional 
   simulation is performed.\cr
   
   \argData \argDataGiven

   If the argument
 \code{x} is missing,
 \code{data} may contain \code{NA}s, which are then replaced by
 conditionally simulated values (random imputing);
 }
 \item{given}{\argGiven }
 
 \item{err.model,err.params}{For conditional simulation and random
 imputing only.\cr\argErrmodel.}

 \item{n}{number of realizations to generate.
   For a very advanced feature, see the notes in \link{RFsimulateAdvanced}.
 }
  \item{...}{\argDots}
}


\details{
 By default, all Gaussian random fields have zero mean. 
 Simulating with trend can be done by including \command{\link{RMtrend}}
 in the model, see the examples below.


 If \code{data} is passed, conditional simulation based on
 simple kriging is performed: 
 \itemize{
 
 \item
 If of class \code{\link[=RFsp-class]{RFsp}},
 \code{ncol(data@coords)} must equal the dimension of the index
 space. If \code{data@data} contains only a single variable,
 variable names are optional. If \code{data@data} contains
 more than one variable, variables must be named and \code{model}
 must be given in the tilde notation \code{resp ~ ... } (see
 \command{\link{RFformula}}) and \code{"resp"} must be contained
 in \code{names(data@data)}.
 \item
 % Beschreibung hier stimmt nicht so ganz mit Examples unten ueberein
 If \code{data} is a matrix or a data.frame, either \code{ncol(data)}
 equals \eqn{(dimension of index space + 1)} and the order of the
 columns is (x, y, z, T, response) or, if \code{data} contains
 more than one 
 response variable (i.e. \code{ncol(data) > (dimension of index
 space + 1)}), \code{colnames(data)} must contain
 \code{colnames(x)} or those of \code{"x", "y", "z", "T"} that
 are not missing. The response variable name is matched with
 \code{model}, which must be given in the tilde notation. If
 \code{"x", "y", "z", "T"} are missing and \code{data} contains
 \code{NA}s, \code{colnames(data)} must contain an element which starts
 with \sQuote{data}; the corresponding column and those behind it are
 interpreted as the given data and those before the corresponding
 column are interpreted as the coordinates.
 \item
 If \code{x} is missing, \command{\link{RFsimulate}} searches for
 \code{NA}s in the data and performs a conditional simulation
 for them.
 }

 Specification of \code{err.model}:
 In geostatistics we have two different interpretations of a nugget
 effect: small scale variability and measurement error.
 The result of conditional simulation usually does not include the
 measurement error. Hence the measurement error \code{err.model}
 must be given separately. For sake of generality, any model (and not
 only the nugget effect) is allowed.
 Consequently, \code{err.model} is ignored
 when unconditional simulation is performed.

 
}

\value{By default, 
 an object of the virtual class \command{\link[=RFsp-class]{RFsp}};
 result is of class \code{\link[=RMmodel-class]{RMmodel}}.
 
 \itemize{
 \item
 \command{\link[=RFspatialGridDataFrame]{RFspatialGridDataFrame}}
 if the space-time dimension is greater than 1
 and the coordinates are on a grid,
 \item
 \command{\link[=RFgridDataFrame]{RFgridDataFrame}}
 if the space-time dimension equals 1 and the coordinates are on a grid,
 \item
 \command{\link[=RFspatialPointsDataFrame]{RFspatialPointsDataFrame}}
 if the space-time dimension is greater than 1 and the coordinates are not on a grid,
 \item
 \command{\link[=RFpointsDataFrame]{RFpointsDataFrame}}
 if the space-time dimension equals 1 and the coordinates are not on a
 grid.
 }

 In case of a multivariate 
 
 If \code{n > 1} the repetitions make the last dimension.

 See \link{RFsimulateAdvanced} for additional options.
 
}

\references{
 \litLantue
 
 \litIntro

 See \link{RFsimulateAdvanced} for more specific literature.
}

\note{Several advanced options can be found in sections \sQuote{General
    options} and \sQuote{coords} of \command{\link{RFoptions}}.
  In particular, option \code{spConform=FALSE} leads to a simpler
  (and faster!) output, see  \command{\link{RFoptions}} for details.
}

\me

\seealso{
 \command{\link{RFvariogram}},
 \command{\link{RFfit}},
 \command{\link{RFgetModelInfo}},
 \command{\link{RFgui}},
 \command{\link{RMmodel}},
 \command{\link{RFoptions}},
 \command{\link{RFsimulateAdvanced}},
 \command{\link{RFsimulate.more.examples}}
}

\keyword{spatial}


\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

#############################################################
## ##
## ONLY TWO VERY BASIC EXAMPLES ARE GIVEN HERE ##
## see ##
## ?RMsimulate.more.examples ##
## and ##
## ?RFsimulateAdvanced ##
## for more examples ##
## ##
#############################################################

#############################################################
## ##
## Unconditional simulation ## 
## ##
#############################################################

## first let us look at the list of implemented models
RFgetModelNames(type="positive definite", domain="single variable",
                iso="isotropic") 

## our choice is the exponential model;
## the model includes nugget effect and the mean:
model <- RMexp(var=5, scale=10) + # with variance 4 and scale 10
 RMnugget(var=1) + # nugget
 RMtrend(mean=0.5) # and mean
 
## define the locations:
from <- 0
to <- 20
x.seq <- seq(from, to, length=200) 
y.seq <- seq(from, to, length=200)

simu <- RFsimulate(model, x=x.seq, y=y.seq)
plot(simu)



#############################################################
## ##
## Conditional simulation ## 
## ##
#############################################################

# first we simulate some random values at 
# 100 random locations:
n <- 100
x <- runif(n=n, min=-1, max=1)
y <- runif(n=n, min=-1, max=1)
dta <- RFsimulate(model = RMexp(), x=x, y=y, grid=FALSE)
plot(dta)

# let simulate a field conditional on the above data
L <- if (interactive()) 100 else 5
x.seq.cond <- y.seq.cond <- seq(-1.5, 1.5, length=L)
model <- RMexp()
cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=dta)
plot(cond, dta)
\dontshow{FinalizeExample()}}
