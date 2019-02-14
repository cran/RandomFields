\name{Specific}
\alias{Specific}
\alias{RPspecific}
\title{Methods that are specific to certain covariance models}

\description{
 This model determines that the (Gaussian) random field should
 be modelled by a particular method that is specific to the given
 covariance model.
}

\usage{
RPspecific(phi, boxcox)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
   specifies the covariance model to be simulated.}
\item{boxcox}{the one or two parameters of the box cox transformation.
  If not given, the globally defined parameters are used.
  See \command{\link{RFboxcox}} for details.
 }
% \item{loggauss}{optional arguments; same meaning as for 
% \command{\link{RPgauss}}.}
}

\details{
  \code{RPspecific} is used for specific algorithms or specific features
  for simulating certain covariance functions.
  
  \itemize{ % i.W. alle Modele mit struct und do Funktion
    \item{\command{\link{RMplus}}}{
      is able to simulate separately
      the fields given by its summands. This is necessary, e.g., when
      a trend model \command{\link{RMtrend}} is involved.
    }

    \item{\command{\link{RMmult}}} {
      for Gaussian random fields only.
      \command{RMmult} simulates the random fields
      of all the components and multiplies them. This is repeated
      several times and averaged.
    }
    
    \item{\command{\link{RMS}}}{
      Then, for instance,
      \code{sqrt(var)} is multiplied onto the (Gaussian) random
      field after the field has been simulated.
      Hence, when \code{var} is random, then for each realization
      of the Gaussian field (for \code{n>1} in \command{\link{RFsimulate}})
      a new realization of \code{var} is used.
      
      Further, new coordinates are created where the old coordinates
      have been divided by the \code{scale} and/or multiplied with the
      \code{Aniso} matrix or a \code{proj}ection has been performed.

      \code{\link{RPspecific}(\link{RMS}())} is called internally when
      the user wants to simulate \code{Aniso}tropic fields with
      isotropic methods, e.g. \command{\link{RPtbm}}.          
    }
    \item{\command{\link{RMmppplus}}}{
    }
    \item{\command{\link{RMtrend}}}{
    }
    % \item{\command{\link{RM}}}{}
  }

  Note that \code{RPspecific} applies only to the first model or
  operator in the argument \code{phi}.
}

\value{
 \command{RPspecific} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

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

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## example for implicit use
model <- RMgauss(var=10, s=10) + RMnugget(var=0.1)
plot(model)
plot(RFsimulate(model=model, 0:10, 0:10, n=4))
## The following function shows the internal structure of the model.
## In particular, it can be seen that RPspecific is applied to RMplus.
RFgetModelInfo(level=0, which="internal")

## example for explicit use: every simulation has a different variance
model <- RPspecific(RMS(var=unif(min=0, max=100), RMgauss()))
x <- seq(0,50,0.02)
plot(RFsimulate(model, x=x, n=4), ylim=c(-15,15))

\dontshow{FinalizeExample()}}

\seealso{ \link{Gaussian},
  \link{RP}.
}

\keyword{methods}
