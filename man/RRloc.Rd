\name{RRloc}
\alias{RRloc}
\title{Location and Scale Modification of A Distribution}
\description{

 \command{RRloc} modifies location and scale of a distribution.
 
 }
 \usage{
RRloc(phi, mu, scale, pow) 
}

\arguments{
  \item{phi}{distribution \command{\link{RMmodel}}}
  \item{mu}{location shift}
  \item{scale}{scale modification}
  \item{pow}{argument for internal use only} % 
}

\value{
 \command{\link{RRloc}} returns an object of class \code{\link[=RMmodel-class]{RMmodel}}.

}

\details{
  It has the same effect as
  \code{\link{RRdistr}(\link{norm}(mu=mu, sd=sd, log=log))}
}

\me

\seealso{
 \command{\link{RMmodel}},
 \command{\link{RRdistr}},
\command{\link{RRgauss}}.
 }

 \examples{\dontshow{StartExample()}
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again

## empirical density of the distribution 'RRspheric'
model <- RRspheric(balldim=2)
hist(RFrdistr(model, n=1000), 50)

## empirical density of the distribution 'RRspheric', shifted by 3
model <- RRloc(mu=3, RRspheric(balldim=2))
hist(RFrdistr(model, n=1000), 50)
\dontshow{FinalizeExample()}}

\keyword{spatial}
\keyword{models}




