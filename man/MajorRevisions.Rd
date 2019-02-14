\name{Major Revisions}
\alias{MajorRevisions}
\title{Documentation of major changings} 

\description{
  This man page documents some major changings in RandomFields.
}



\section{Changes done in 3.1.0 (Summer 2015)}{
  \itemize{
    \item full (univariate) trend modelling
    \item error in particular in \command{\link{RFfit}} corrected
    \item \command{\link{RFfit}} runs much faster now
    \item effects of \code{modus operandi} changed for estimating
  }
}


\section{Corrections done in 3.0.56 (Jan 2015)}{
  \itemize{
    \item log Gauss field corrected (has been a log log Gauss field)
    \item \command{RMconstant} is now called \command{\link{RMfixcov}}
  }
}


\section{Corrections done in 3.0.55 (Jan 2015)}{
  \itemize{
    \item Conditional simulation: several severe typos corrected.
  }
}

\section{Major Revision: changings from Version 2 to Version 3 (Jan 2014)}{
  \itemize{
    \item \bold{S4 objects}
    \itemize{
      \item
      \pkg{RandomFields} is now based on S4 objects using the package
      \pkg{sp}. The functions accept both \pkg{sp} objects and simple
      objects as used in version 2. See also above.      
    }
    
    \item\bold{Documentation}
    \itemize{
      \item each model has now its own man page;
      \item classes of models and functions are bundled in several pages:
      Covariance models start with \code{\link{RM}}, distribution families
      with \code{\link{RR}}, processes with \code{\link{RP}}, user
      functions with \code{\link{RF}}
      \item the man pages of several functions are split into two parts:
      
      (i) a beginners man page which includes a link to
      
      (ii) man pages for advanced users
    }

   
    \item\bold{Interfaces}
    \itemize{
      \item The interfaces become simpler, at the same time more powerful
      than the functions in version 2. E.g., \code{RFsimulate} can perform
      unconditional simulation, conditional simulation and random
      imputing.
      
      \item Only those arguments are kept in the functions that are considered
      as being absolutely necessary. All the other arguments can be
      included as \link[=RFoptions]{options}.
      
      \item \command{\link{RFgui}} is an instructive interface based on tcl/tk,
      replacing the former \code{ShowModels}
    }
    
    
    \item \bold{Inference for Gaussian random fields}    
    \itemize{
      \item \command{\link{RFfit}} has undergone a major revision. E.g.:
    
      (i) estimation of random effect models with spatial covariance
      structure

      (ii) automatic estimation of 10 and more arguments in
      multivariate and/or space-time models

      \item \command{\link{RFvariogram}} is now based on an fft algorithm
      if the data are on a grid, even allowing for missing values.

      \item \command{\link{RFratiotest}} has been added.
    }
    
    
    \item \bold{Processes}
   \itemize{
     \item Modelling of \link[=maxstable]{maxstable
       processes} has been enhanced, including
     
     (i) the simulation of Brown-Resnick processes
     
     (ii) initial support of
     \link[=RMmodelsAdvanced]{tail correlation functions};
     \item{Further processes}
     \link[=RPchi2]{chi2 processes}, \link[=RPpoisson]{compound
       Poisson processes}, \link[=RPbernoulli]{binary processes} added.
   }
   
   \item \bold{Models}
   \itemize{
     \item the \link[=RFformula]{formula notation} for
     linear models may now be defined
     \item{Novel, user friendly definition of the covariance models}
     \item \link[=RMmodelsAdvanced]{Multivariate and vector-valued random
       fields} are now fully included 
     \item The \link[=RMuser]{user}
       may now define his own functions, to some extend.
     \item The \link[=RMtrend]{trend} allows for much more flexibility
     \item \link[=RR]{Distributions} may now be included which will be
     extended to \link{Bayesian} modelling in future.
   }
 }
}

\me

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## S4 vs S3
x <- seq(0, 10, 0.1)
model <- RMexp()
plot(RFsimulate(model, x)) ## S4
plot(RFsimulate(model, x, spConform=FALSE)) ## no class

\dontshow{FinalizeExample()}}
