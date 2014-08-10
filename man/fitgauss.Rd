\name{fitgauss}
\alias{fitgauss}
%\alias{AIC}
%\alias{AICc}
%\alias{logLik}
\title{Details on fitting Gaussian random fields}
\description{
 Here some details of \command{\link{RFfit}} are given
 concerning the fitting of models for Gaussian random fields

\bold{  This documentation is far from being complete. }
}
% to do : merge with privat/rffitAdvanced.Rd

%\usage{}
%\arguments{}
%\value{}
\section{Maximum likelihood}{
  The application of the usual maximum likelihood method and
  reporting the result is the default.
}

%\section{Restricted maximum likelihood}{
% to do
%}

%\section{Composed likelihood}{
% to do
%}

\section{Least squares}{
 
 The weighted least squares methods minimise
 \deqn{
 \sum_{i} w_i (\hat \gamma(h_i) - \gamma(h_i))^2
 }
 over all parametrised models of \eqn{\gamma}.
 Here, \eqn{i} runs over all \eqn{N} bins of the binned variogram
 \eqn{\hat \gamma} 
 and \eqn{h_i} is the centre of bin \eqn{i}.

 The following variants of the least squares methods,
 passed as \code{sub.methods} in \command{\link{RFfit}} are implemented.
 \describe{
 \item{\code{'self'}}{
 \eqn{w_i = (\gamma(h_i))^{-2}}
 }
 \item{\code{'plain'}}{\eqn{w_i = 1} for all \eqn{i}.
 }
 \item{\code{'sqrt.nr'}}{\eqn{w_i^2} equals the number of points
 \eqn{n_i} in bin \eqn{i}.
 }
 \item{\code{'sd.inv'}}{\eqn{1 / w_i} equals the standard deviation
 of the variogram cloud within bin \eqn{i}.
 }
 \item{\code{'internal'}}{
 Three subvariants are implemented:
 \describe{
	\item{\code{'internal1'}}{
	 \eqn{w_i^2 = (N-i+1) n_i } 
	}
	\item{\code{'internal2'}}{
	 \eqn{w_i = N-i+1}
	}
	\item{\code{'internal3'}}{
	 \eqn{w_i^2 = N-i+1}
	}
      }
    }
  }
}


\seealso{
 \code{\link{RFfit}}, \code{RFfit-class}
}
\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
## see 'RFfit'
\dontshow{FinalizeExample()}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
  \url{http://ms.math.uni-mannheim.de/de/publications/software}}

\keyword{spatial}
\keyword{optimize}
