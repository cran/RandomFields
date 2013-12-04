\name{RFgetModelNames}
\alias{RFgetModelNames}
%\alias{PrintModelList}
\title{Names of implemented covariance and variogram models}
\description{Displays the names of covariance and
 variogram models (see \command{\link{RMmodel}}) and returns them as a
 list. The user may
 specify and group the models according to the following properties:
 \itemize{
 \item types of stationarity and isotropy
 \item whether the model is an operator
 \item whether the model is a normal scale mixture
 \item whether the model has a finite range covariance
 \item validity in certain dimensions of the coordinate space 
 \item maximal possible dimension of the coordinate space
 \item uni- or multivariety
 }
 See \command{Details} for an explanation and
 \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
 for possible states (values) of these properties.
}

\usage{RFgetModelNames(type = ZF_TYPE, domain = ZF_DOMAIN,
                       isotropy = ZF_ISOTROPY, operator = c(TRUE, FALSE),
                       normalmix = c(TRUE, FALSE),
                       finiterange = c(TRUE, FALSE),
                       valid.in.dim = c(1, Inf), 
                       vdim = c(1, 5),
                       group.by=NULL,
                       internal, newnames
)
}

\arguments{
  \item{type, domain, isotropy, operator, normalmix, finiterange, vdim}{
    see \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}}.
  }
  \item{valid.in.dim}{an optional integer indicating the dimension of
  the space where the model is valid}
  \item{group.by}{an optional character string; must be one of
    \code{'stationarity'}, \code{'isotropy'}, \code{'operator'},
    \code{'normalmix'},
    \code{'finiterange'},\code{'maxdim'},\code{'vdim'}}
  \item{internal, newnames}{both logical;
    \code{internal} might be also integer valued.
    If any of them are given,
    \command{\link{RFgetModelNames}} behaves very differently.
    See the Notes below.
  }
}

\note{
  In case \code{internal} or \code{newnames} is given,
  only the values of \code{internal},
  \code{newnames} and \code{operator} are considered.
  All the other parameters are ignored and
  \command{\link{RFgetModelNames}} prints a table of the currently
  implemented covariance functions and the matching methods:
  \itemize{
    \item \code{internal}:\cr
    if \code{TRUE} also \code{\link{RMmodels}} are listed that are
    internal, hence invisible to the user. Default: \code{FALSE}.
    \item \code{newnames}:\cr
    The model names of version 2 of \pkg{RandomFields} and earlier
    can still be used in the model definitions. Namely when the
    list notation is chosen; see  \link{Advanced RMmodels} for the
    latter. If the \code{internal} or \code{newnames} is given, then
    these old names are shown; if \code{newnames=TRUE} then also the
    usual names are shown.
    Default: \code{FALSE}.

    In fact, both internal and public
    models can have different variants implemented.
    These variants are also shown if
    \code{internal} has a value greater than or equal to\code{2},
    \item \code{operator}:\cr
    see above.
  }
  Here, also an indication is given, which method for simulating
  Gaussian random fields matches the model.
}

\details{ 
 The plain call \code{\link{RFgetModelNames}()} simply gives back a vector of
 the names of all implemented covariance and variogram models and operators,
 i.e. members of the class
 \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}.

 The following arguments can be specified:
 
 \describe{
   \item{\code{type}}{specifies the class of functions; for
     the meaning of the possible values see
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
   }
   
   \item{\code{stationarity}}{specifies the type of stationarity; for
     the meaning of the possible values see
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
   }
   
   \item{\code{isotropy}}{specifies the type of isotropy; for
     the meaning of the possible values see
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
   }
   
   \item{\code{operator}}{indicates whether the model is an operator, 
     i.e. it requires at least one submodel, 
     e.g. \command{\link{RMplus}} or
     \command{\link{RMdelay}} are operators; see
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
   }
   
   \item{\code{normalmix}}{indicates whether the model is a normal
     scale mixture, e.g. \command{\link{RMexp}} or
     \command{\link{RMcauchy}} are normal scale mixtures; see
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
   }
   
   \item{\code{finiterange}}{indicates whether the covariance of the
     model has finite range, e.g. \command{\link{RMcircular}} or
     \command{\link{RMnugget}} have covariances with finite range; see
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
   }
   
   \item{\code{valid.in.dim}}{If \code{valid.in.dim=n} is
     passed, all models which are valid in dimension \eqn{n} are
     displayed. Otherwise \code{valid.in.dim} should be bivariate vector
     giving the range of requested dimensions.
   }
   
   \item{\code{maxdim}}{if a positive integer, it specifies the maximal
     possible dimension of the coordinate space;
     note that a model which is valid in dimension
     \eqn{n} is also valid in dimension \eqn{n-1};
     \code{maxdim=-1} means that the maximal possible dimension depends
     on the parameters of the \command{\link{RMmodel}} object;
     \code{vdim=-2} means that the maximal possible dimension is
     adopted from the called submodels;
     see also
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}
   }
 
   \item{\code{vdim}}{if a positive integer, vdim specifies, whether
     the model is \eqn{vdim}-variate; 
     \code{vdim=-1} means that being multivariate
     in a certain dimension depends on the parameters of the
     \command{\link{RMmodel}} object;
     \code{vdim=-2} means that being multivariate
     in a certain dimension is adopted from the called submodels;
     see also
     \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}

     If \code{vdim} is bivariate then a range is given.
   }
   
   \item{\code{group.by}}{If \code{group.by="propertyname"} is
     passed, the displayed models are grouped according to
     \code{propertyname}.
   }
 }

 All parameters allow also for vectors of values. In case of
 \code{valid.in.dim} the smallest value is taken.
 The interpretation is canonical.
 
 
 Note that the arguments \command{stationarity}, \command{isotropy},
 \command{operator}, \command{normalmix}, \command{finiterange},
 \command{maxdim}, \command{vdim}
 are also slots (attributes) of the SP4-class
 \command{\link[=RMmodelgenerator-class]{RMmodelgenerator}}.
}

\value{
  Either a vector of model names if the argument
 \command{group.by} is not used;
 or a list of vectors of model names if the argument \command{group.by} is
 used
 (with list elements specified by the categories of the grouping
 argument).

 In case \code{internal} or \code{newnames} is given,
 \command{\link{RFgetModelNames}} prints a table of the currently
 implemented covariance functions and the matching methods.
 \command{\link{RFgetModelNames}} returns \code{NULL}.
}

%\references{
%}

\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\seealso{
 \code{\link[=RMmodelgenerator-class]{RMmodelgenerator}},
 \command{\link{RMmodel}},
 \code{\link[=RandomFields-package]{RandomFields}}.
}

\examples{
set.seed(0)
# get vector of names of all objects of class RMmodelgenerator
RFgetModelNames()

# get vector of names of all stationary objects of class RMmodelgenerator
RFgetModelNames(type="positive definite", domain="single variable")

# get list of models grouped by the stationarity attribute
RFgetModelNames(group.by=c("type"))


\dontrun{
# get list of all univariate stationary models
# additionally grouped by the isotropy attribute
str(RFgetModelNames(type="positive definite", domain="single variable",
                    vdim=1, group.by="isotropy"))

# get vector of all models which are operators
# and valid in the two-dimensional coordinate space
RFgetModelNames(type=c("tail correlation function",
                       "positive definite",
                       "negative definite",
                       "undefined"),
                operator=TRUE, valid.in.dim=2)

# processes and covariance function grouped by the stationarity
# argument and subsequently grouped by the isotropy argument
str(RFgetModelNames(type=c("positive definite", "negative definite", "process"),
                    group.by=c("type", "domain", "isotropy")))
}


}
\keyword{spatial}
