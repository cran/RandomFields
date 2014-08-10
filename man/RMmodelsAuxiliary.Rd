\name{Others}
\alias{AuxiliaryModels}
\alias{Auxiliary Models}
\alias{Auxiliary RMmodels}
\alias{Other models}
\alias{RMmodelsAuxiliary}
\title{Auxiliary and other Models}
\description{
 Here, auxiliary models are given that are not covariance functions or
 variograms, but which might be used within the definition of a model.
}

\section{Implemented models}{
  \bold{Distribution families}
  See \link{RR}.

  \bold{Evaluation operators}
  See \link{RF}.

  \bold{Random Fields / Random Processes}
  See \link{RP}.  



  \bold{Shape functions}
  
  Besides any of the covariance functions the following functions can be
  used as shape functions.
 \tabular{ll}{
    \command{\link{RMangle}} \tab defines an anisotropy matrix by angle
    and a diagonal matrix.
    \cr
    \command{\link{RMball}} \tab Indicator of a ball of radius
    \eqn{1/2}.
    \cr
    \command{\link{RMm2r}} \tab spectral function belonging to a
    tail correlation function of the Gneiting class \eqn{H_n}
    \cr 
    \command{\link{RMm3b}}  \tab spectral function belonging to a tail
    correlation function of the Gneiting class \eqn{H_n}
    \cr 
    \command{\link{RMmppplus}} \tab operator to define mixed moving
    maxima (M3) processes
    \cr 
    \command{\link{RMmps}} \tab spectral functions belonging to a tail
    correlation function of the Gneiting class \eqn{H_n}
    \cr 
    \command{\link{RMpolygon}} \tab Indicator of a typical Poisson
    polygon. 
    \cr
    \command{\link{RMrational}} \tab shape function used in the Bernoulli
    paper given in the references
    \cr
    \command{\link{RMrotat}} \tab shape function used in the Bernoulli
    paper given in the references
    \cr
    \command{\link{RMsign}} \tab random sign
    \cr
    \command{\link{RMtruncsupport}} \tab truncates the support of a shape in
    a Poisson based model
    \cr
  }

  \bold{Special transformations}
  \tabular{ll}{
    \command{\link{RMeaxxa}} \tab shape function used in the Bernoulli
    paper given in the references
    \cr
    \command{\link{RMetaxxa}} \tab shape function used in the Bernoulli
    paper given in the references
    \cr
    \command{\link{RMid}} \tab identity \cr
  \command{\link{RMrotation}} \tab shape function used in the Bernoulli
   paper given in the references
   \cr
}
 
  \bold{Other models}
  \tabular{ll}{
    \command{\link{RMcoord}} \tab passing new coordinates in a mixed
    model 
    \cr
    \command{\link{RMuser}} \tab User defined covariance model
    \cr
  }



}


%\section{Methods}{
% \describe{
% \item{[}{\code{signature(x = "RFgridDataFrame")}: selects %]
% slot by name}
% \item{[<-}{\code{signature(x = "RFgridDataFrame")}: replaces
% slot by name}
% \item{as}{\code{signature(x = "RFgridDataFrame")}:
% converts into other formats, only implemented for target class
% \command{\link[=RFpointsDataFrame-class]{RFpointsDataFrame}} } 
% \item{cbind}{\code{signature(...)}: if arguments have identical
% topology, combine their attribute values}
% }
%}


\seealso{\command{\link{RFformula}}, \command{\link{RMmodels}},
 \command{\link{RMmodelsAdvanced}}
}

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
(type <- RC_TYPE[(0:1) -length(RC_TYPE)])
RFgetModelNames(type=type, group.by="type")
\dontshow{FinalizeExample()}
}

\author{
 Alexander Malinowski, \email{malinowski@math.uni-mannheim.de}
 
 Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}
