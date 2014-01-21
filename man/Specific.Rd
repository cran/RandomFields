
\name{Specific}
\alias{Specific}
\alias{RPspecific}
\title{Methods that are specific to certain covariance models}

\description{
 This model defines that the (Gaussian) random field should
 be modelled by a particular method that is specific to the given
 covariance model.

 \bold{These is a very advanced model of \pkg{RandomFields}.}
}

\usage{
RPspecific(phi, loggauss)
}

\arguments{
 \item{phi}{object of class \code{\link[=RMmodel-class]{RMmodel}};
 specifies the covariance model to be simulated.}
 \item{loggauss}{optional parameters; same meaning as for 
 \command{\link{RPgauss}}.}
}

\details{
  \code{RPspecific} is used for specific algorithms or specific features
  for simulating certain covariance functions
  
  \itemize{ % i.W. alle Modele mit struct und do Funktion
    \item{\command{\link{RMplus}}}{
      \command{\link{RMplus}} is able to simulate separately
      the fields given by its summands. This is necessary, e.g., when
      a \command{\link{RMtrend}} is involved.
    }
    \item{\command{\link{RMS}}}{
      Then the \code{sqrt(var)} is multiplied onto the Gaussian random
      fields after the field has been simulated.
      Hence, when \code{var} is random, then, for each realisation
      of the Gaussian field (for \code{n>1} in \command{\link{RFsimulate}})
      a new realisation of \code{var} is used.
      
      Further, new coordinates are created where the old coordinates
      have been devided by the \code{scale} and/or multiplied with the
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
}

\value{
 \command{RPspecific} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}

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
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}

\examples{
RFoptions(seed=0)
## example for implicite use
model <- RMgauss(var=10, s=10) + RMnugget(var=0.1)
plot(model)
plot(RFsimulate(model=model, 0:10, 0:10, grid=TRUE, n=4))
## The following function shows the internal structure of the model.
## In particular, it can be seen that RPspecific is applied to RMplus.
str(RFgetModelInfo(level=0, which="keys"))

## example for explicite used
model <- RPspecific(RMS(var=unif(min=0, max=10), RMgauss()))
n <- if (interactive()) 6 else 1
k <- if (interactive()) 10 else 1
x <- seq(0,10,0.02)
\dontshow{RFoptions(seed=NA)}
for (i in 1:k) {
  readline(paste("Simulation no.", i, ": press return", sep=""))
  plot(RFsimulate(model, x=x, n=n), ylim=c(-5,5))
}

\dontshow{RFoptions(seed=NA)}
}

\seealso{
  \link{RP},
 \command{\link{RPcoins}},
 \command{\link{RPhyperplane}},
 \command{\link{RPspectral}},
 \command{\link{RPtbm}}.
}

\keyword{methods}
