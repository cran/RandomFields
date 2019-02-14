\name{Independent Variables}
\alias{RPnugget}
\alias{Nugget}
\title{Method to simulate the Nugget effect}



\description{
 Method to simulate the Nugget effect. (Only for advanced users)
}

\usage{
RPnugget(phi, boxcox, tol, vdim)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 specifies the covariance model to be simulated. The only possible
 model for \code{phi} is \command{\link{RMnugget}}.}
 \item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  See \command{\link{RFboxcox}} for details.
 }
\item{tol}{
   points at a distance less than or equal to \code{nugget.tol}
 are considered as being identical. This strategy applies to
 the simulation method and the covariance function itself.
 Hence, the covariance function is only positive definite
 if \code{nugget.tol=0.0}. However, if the anisotropy matrix
 does not have full rank and \code{nugget.tol=0.0}, then
 the simulations are likely to be odd. 
 The value of \code{nugget.tol}
 should be of order \eqn{10^{-15}}{1e-15}.
 
 Default: \code{0.0}
 }
 \item{vdim}{positive integer; the model is treated
   \code{vdim}-variate, \code{vdim=1} (default) corresponds to a
   univariate random field.
   Mostly, the value of \code{vdim} is set automatically.

   Default is that it takes the value of the submodel \code{phi}.}
}
 

\details{
  \describe{
    \item{General}{
      This method only allows \command{\link{RMnugget}} as a submodel.
    }
    
    \item{Anisotropy}{
      The method also allows for zonal nugget effects. Only there the
      argument \code{tol} becomes important.
      For the zonal nugget effect, the anisotropy matrix \code{Aniso}
      should be given in \command{\link{RMnugget}}. There, only the
      kernel of the 
      matrix is important.
    }

    \item{Points close together}{
      The 
      locations at a distance less than or equal to the \link{RFoptions}
      \code{nugget.tol}
      are considered as being identical. This strategy applies to
      the simulation method and the covariance function itself.
      Hence, the covariance function is only positive definite
      if \code{nugget.tol=0.0}. However, if the anisotropy matrix
      does not have full rank and \code{nugget.tol=0.0}, then the
      simulations are likely to be odd.  
      The value of \code{nugget.tol}
      should be of order \eqn{10^{-15}}{1e-15}.
    }
    
    \item{Repeated measurements}{
      Measurement errors are mathematically not distinguishable from spatial
      nugget effects as long as measurements are not repeated at the very
      same space-time
      location. So there is no need to distinguish the spatial nugget
      effect from a measurement error.
      This is the default, see
      \code{allow_duplicated_locations} in \link{RFoptions}.
      
      In case several measurement have been taken in single space-time
      locations,
      measurement errors can be separated from spatial noise.
      In this case \code{RMnugget()} models the measurement error (which
      corresponds to a non-stationary model in an abstract space) by
      default and the
      measurement error model cannot be extended beyond the given
      locations.
      On the other hand \code{RMnugget(Ansio=something)} and
      \code{RMnugget(proj=something)} model the spatial nugget effect
      (with and without zonal anisotropy in case \code{Aniso} has low and
      full rank respectively).
    }

    \item{Role of \command{RPnugget}}{
      Even for advanced users, there is no need to call
      \command{RPnugget} directly, as this is done internally when
      the \link{RMnugget} is involved in the covariance model.
    }
  }
}

\value{
 \code{RPnugget} returns an object of class
 \code{\link[=RMmodel-class]{RMmodel}}. 

}

\references{
 \itemize{
 \item
 Schlather, M. (1999) \emph{An introduction to positive definite
 functions and to unconditional simulation of random fields.}
 Technical report ST 99-10, Dept. of Maths and Statistics,
 Lancaster University. 
 }
}
\me

\seealso{ \link{Gaussian},
  \link{RP},
 \command{\link{RPcoins}},
 \command{\link{RPhyperplane}},
 \command{\link{RPspectral}},
 \command{\link{RPtbm}}.
}

\keyword{methods}

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


x <- y <- 1:2
xy <- as.matrix(expand.grid(x, y)) ## we get 4 locations
## Standard use of the nugget effect
model <- RMnugget(var = 100)
RFcovmatrix(model, x=xy)
as.vector(RFsimulate(model, x=x, y=x, tol=1e-10))

## zonal nugget effect, which is not along the axes
model <- RMnugget(Aniso=matrix(1, nr=2, nc=2))
RFcovmatrix(model, x=xy)
as.vector(RFsimulate(model, x=x, y=x, tol=1e-10))


## All the following examples refer to repeated measurements
RFoptions(allow_duplicated_locations = TRUE) 
(xy <- rbind(xy, xy)) ## now, the 4 locations are repeated twice 

## standard situation: the nugget is interpreted as measurement error:
model <- RMnugget()
RFcovmatrix(model, x=xy)
as.matrix(RFsimulate(model, x=xy)) 

## any anisotropy matrix with full rank: spatial nugget effect
model <- RMnugget(Aniso=diag(2))
RFcovmatrix(model, x=xy)
as.matrix(RFsimulate(model, x=xy))

## anisotropy matrix with lower rank: zonal nugget effect
model <- RMnugget(Aniso=matrix(c(1, 0, 0, 0), nc=2))
RFcovmatrix(model, x=xy)
as.matrix(RFsimulate(model, x=xy))

## same as before: zonal nugget effect
model <- RMnugget(Aniso=t(c(1,0)))
RFcovmatrix(model, x=xy)
as.matrix(RFsimulate(model, x=xy))
\dontshow{FinalizeExample(); RFoptions(allow_duplicated_locations = FALSE) }}