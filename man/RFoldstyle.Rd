\name{RFoldstyle}
\alias{RFoldstyle}
\title{RFoldstyle}
\description{
 This function is written only for package writers who have based
 their code on RandomFields version 2.

 It avoids warnings if the old style is used, and sets
 \code{spConform = FALSE}.
 
 }
\usage{
RFoldstyle(old=TRUE) 
}
\arguments{
 \item{old}{logical}
}
\value{
 \code{NULL}
}

%\note{
%  see for partial documentation of the former
%  functions. Note that the documentation is not maintained.
%}
\seealso{
 See \sQuote{\href{../doc/version.pdf}{version2}}
 for details on the commands of version 2.
}

\me

\examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

GaussRF(x=1:10, model="exp", param=c(0,1,0,1), grid=TRUE)

RFoldstyle()
GaussRF(x=1:10, model="exp", param=c(0,1,0,1), grid=TRUE)

\dontshow{RFoldstyle(FALSE)}
\dontshow{FinalizeExample()}}

\keyword{spatial}


