\name{RMmodelgenerator-class}
\docType{class}
\alias{RMmodelgenerator-class}
\alias{show,RMmodelgenerator-method}
\alias{[,RMmodelgenerator-method}
\alias{[,RMmodelgenerator,ANY,ANY-method}
\alias{[,RMmodelgenerator,ANY,ANY,ANY-method}
\alias{[<-,RMmodelgenerator,ANY,ANY-method}
\alias{[<-,RMmodelgenerator,ANY,ANY,ANY-method}
\alias{[<-,RMmodelgenerator-method}
\alias{print.RMmodelgenerator}

\title{Class \code{RMmodelgenerator} }
\description{ Class for all functions of this package with prefix \code{RM},
  i.e. all functions that generate objects of class
  \code{\link[=RMmodel-class]{RMmodel}}; direct extension of
  class \code{\link[methods:function-class]{function}}
}

 
\section{Creating Objects}{
 Objects should not be created by the user!
}

\section{Slots}{
 \describe{
 \item{\code{.Data}:}{function; the genuine funtion that generates an
 objects of class
 \command{\link[=RMmodel-class]{RMmodel}} }
 \item{\code{type}:}{character string; specifies the cathegory of
 RMmodel-function, see Details} 
 \item{\code{domain}:}{character string; specifies whether the 
   corresponding function(s) depend on 1 or 2 variables, see Details} 
 \item{\code{isotropy}:}{character string; specifies the type of
 isotropy of the corresponding covariance model, see Details} 
 \item{\code{operator}:}{logical; specifies whether the underlying
 covariance model is an operator, see Details} 
 \item{\code{monotone}:}{character string; specifies the kind of monotonicity
 of the model}
 \item{\code{finiterange}:}{logical; specifies whether the underlying
   covariance model has finite range, see Details}
 % finiterange waere wuenschenswert an Abhaengigkeit von parametern
 % und submodellen anzupassen, siehe maxdim und vdim
 \item{\code{simpleArguments}:}{logical. If \code{TRUE} than all the
   parameters are real valued (or integer valued).
 }
 \item{\code{maxdim}:}{numeric; the maximal dimension, in which the
 corresponding model is a valid covariance model, see Details} 
 \item{\code{vdim}:}{numeric; dimension of the value of the random
 fiels at a single fixed location, equals 1 in most cases, see
 Details} 
 }
}

\section{Extends}{
 Class \code{\link[methods:function-class]{function}}, directly.
}

\section{Methods}{
 \describe{
 \item{show}{\code{signature(x = "RMmodel")}: returns the structure
 of \code{x}}
 \item{print}{\code{signature(x = "RMmodel")}: identical with
 \command{show}-method}
 \item{[}{\code{signature(x = "RMmodelgenerator")}: enables accessing
 the slots via the "["-operator, e.g. x["maxdim"]}
 \item{[<-}{\code{signature(x = "RMmodelgenerator")}: enables replacing
 the slots via the "["-operator}
 }
}


\section{Details}{
 \describe{
 \item{\code{type}:}{can be one of the following strings:
   \describe{
     \item{\code{'tail correlation function'}:}{
       indicates that the function returns a
       tail correlation function (a subclass of the set of
       positive definite functions)
       }
     \item{\code{'positive definite'}:}{indicates that the function returns a
       covariance function (positive definite function)}
     \item{\code{'negative definite'}:}{indicates that the function returns a
       variogram model (negative definite function)}
     \item{\code{'process'}:}{functions of that type determine
       the class of 
       processes to be simulated} 
     \item{\code{'method for Gauss processes'}:}{methods to simulate
       Gaussian random fields}
     \item{\code{'method for Brown-Resnick processes'}:}{methods to
       simulate Brown-Resnick fields}
     \item{\code{'point-shape function'}:}{functions of that type
       determine the distribution of points in space}
     \item{\code{'distribution family'}:}{
       e.g. (multivariate) uniform distribtion, normal distribution,
       etc.,
       defined in \pkg{RandomFields}.
       See \link{RR} for a complete list.
     }
     \item{\code{'shape function'}:}{functions used in, e.g., M3 processes (\link{RPsmith})}
     \item{\code{'trend'}:}{\link{RMtrend} or a \link[=RFformula]{mixed model} }
     \item{\code{'interface'}:}{indicates internal models which are usually
       not visible for the users. These functions are the internal
       representations of \command{\link{RFsimulate}},
       \command{\link{RFcov}}, etc.. See \link{RF} for a complete list.
     }%\item{\code{'undefinded'}:}{}
     \item{\code{'undefined'}:}{some models can take different types,
 depending on the parameter values and/or the submodels }
    \item{\code{'other type'}:}{very very special internal functions,
      not belonging to any of the above types.
    }
   }
 }
 \item{\code{domain}:}{can be one of the following strings:
   \describe{
     \item{\code{'single variable'}:}{Function depending on a single variable}
     \item{\code{'kernel'}:}{model refers to a kernel, e.g., an
       non-stationary convariance function} 
     \item{\code{'framework dependent'}:}{domain depends on the calling model}
     \item{\code{'mismatch'}:}{this option is used only internaly and should never appear}
   }
 }
 \item{\code{isotropy}:}{can be one of the following strings:
   \describe{
     \item{\code{'isotropic'}:}{indicates that the model is isotropic}
     \item{\code{'space-isotropic'}:}{indicates that the spatial part of a
       spatio-temporal model is isotropic} 
     \item{\code{'zero-space-isotropic'}:}{this property refers to space-time
       models; the model is called zerospaceisotropic if it is
       isotropic as soon as the time-component is zero}
     \item{\code{'vector-isotropic'}:}{multivariate vector model (flow
 fields) have a different notion of isotropy}
     \item{\code{'symmetric'}:}{the most basic property of any
       covariance function or variogram model}
     \item{\code{'cartesian system'}, \code{'earth system'},
	 \code{'spherical system'}, \code{'cylinder system'}:}{
	 different coordinate systems
     }
     \item{\code{'non-dimension-reducing'}:}{the property \eqn{f(x)
	 = f(-x)^\top} does not hold
     }
     \item{\code{'parameter dependent'}:}{indicates that the type of
       isotropy of the model depends on the parameters passed to the
       model; in particular parameters may be submodels if an operator model
       is considered}
    \item{\code{'<mismatch>'}:}{this option is used only internaly and should never appear}
   }
 }
 \item{\code{operator}:}{if \code{TRUE}, the model requires at least
 one submodel}
 \item{\code{monotone}:}{
   \describe{
     \item{\code{'mismatch in monotonicity'}:}{used if a statement on
       the monotonocity does not make sense, e.g. for
       \code{\link{RRmodels}}
       }
     \item{\code{'submodel dependent monotonicity'}:}{only for operators,
 e.g. \code{\link{RMS}}}
     \item{\code{'previous model dependent monotonicity'}:}{internal;
       should not be used}
     \item{\code{'parameter dependent monotonicity'}:}{some models change
 their properties according to the parameters}
     \item{\code{'not monotone'}:}{none of the above cathegories; either
 the function is not monotone or properties are not known}
     \item{\code{'monotone'}:}{isotone or antitone}
     \item{\code{'Gneiting-Schaback class'}:}{function belonging to
 Euclid's hat in Gneiting's 1999 paper} 
     \item{\code{'normal mixture'}:}{scale mixture of the Gaussian model} 
     \item{\code{'completely monotone'}:}{completely monotone function}
     \item{\code{'Bernstein'}:}{Bernstein function}
   }
   Note that
   \itemize{
     \item \code{'not monotone'} includes \code{'monotone'}
     and \code{'Bernstein'}
     \item \code{'monotone'} includes \code{'Gneiting-Schaback class'} 
     \item \code{'Gneiting-Schaback class'} includes \code{'normal mixture'} 
     \item \code{'normal mixture'} includes \code{'completely monotone'}
   }
 }
 \item{\code{finiterange}:}{if \code{TRUE}, the covariance of the
 model has finite range}
 \item{\code{maxdim}:}{if a positive integer, \code{maxdim} gives the
 maximum dimension in which the model is a valid covariance model,
 can be \code{Inf};
 \code{vdim=-1} means that the actual maxdim depends on the
 parameters; \code{vdim=-2} means that the actual maxdim depends on
 the submodel(s)} 
 \item{\code{vdim}:}{if a positive integer, \code{vdim} gives the
 dimension of the random field, i.e. univariate, bi-variate, ...;
 \code{vdim=-1} means that the actual vdim depends on the
 parameters; \code{vdim=-2} means that the actual vdim depends on
 the submodel(s)}
 }
}

\author{Alexander Malinowski \email{malinows@math.uni-goettingen.de};
  Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de/de/publications/software}}


\references{
  \itemize{
    \item Gneiting, T. (1999) Radial positive definite functions
    generated by Euclid's hat, \emph{J.  Multivariate Anal.},
    \bold{69}, 88-119.
    }
}

\seealso{
 \code{\link[=RMmodel-class]{RMmodel}},
 \code{\link{RFgetModelNames}}
}

\keyword{classes}
\keyword{hplot}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
RFgetModelNames()
\dontshow{FinalizeExample()}
}
