\name{SBS14}
\alias{SBS14}
%\alias{Papers}
%\alias{schlather}
%\alias{Schlather}
\title{Systematic co-occurrence of tail correlation functions
among max-stable processes
}
\description{
 Here the code of the paper on \sQuote{On some covariance models based on normal scale mixtures} is given.
}
 
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de};
}
\references{
 \itemize{ 
 \item
 Strokorb, K., Ballani, F. and Schlather, M. (2014)
 Systematic co-occurrence of tail correlation functions
 among max-stable processes. Work in progress.
 }

}

\examples{
set.seed(0)
\dontshow{if (!.C("isAuthor", a=integer(1))$a) {}}

\dontshow{\dontrun{
%  library(RandomFields, lib="~/TMP")
RFoptions(every=1000, xi=0)
model <- RMfbm(var=4, alpha=1) ## not 8 as in the paper
x <- seq(0, 5, 0.01)

z <- RFsimulate(RPbrownresnick(model), x, x)
plot(z)



}}
}


\keyword{spatial}

