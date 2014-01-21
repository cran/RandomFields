\name{fitgauss}
\alias{fitgauss}
\title{Details on fitting Gaussian random fields}
\description{
 Here some details of \command{\link{RFfit}} are given
 concerning the fitting of models for Gaussian random fields

 This documentation is far from being complete.
}
% to do : merge with privat/rffitAdvanced.Rd

%\usage{}
%\arguments{}
%\value{}
\section{Maximum likelihood}{
 to do
}

%\section{Restricted maximum likelihood}{
% to do
%}

%\section{Composed likelihood}{
% to do
%}

\section{Least squares}{
  The following least square estimation methods are available:

 The weighted least squares methods minimises
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
%\references{}
\seealso{
 \code{\link{RFfit}}
}
\examples{
RFoptions(seed=0)
## to do
\dontshow{RFoptions(seed=NA)}
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}}
\keyword{spatial}
