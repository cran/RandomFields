\name{Trend Modelling}
\alias{RMmodelsTrend}
\alias{trend modelling}
\title{Trend Modelling}
\description{
  The coding of trend, in particular multivariate trends,
  will be described here
}

\details{
  See \link{RFcalc}, \link{RMtrend} and also the examples below
  for some insight on the possibilities of trend modelling.
}

%\references{Ribeiro}
\seealso{
  \link{RFcalc}, \link{RM}, \link{RMmodels}, \link{RMtrend}, \link{RMmodelsMultivariate}
}

\examples{
%RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
%##                   RFoptions(seed=NA) to make them all random again
\dontshow{StartExample();if (!interactive()) RFoptions(modus="careless")} %ok

require(geoR)

data(ca20) ## data set from geoR
ca20.df <- as.data.frame(ca20)
head(ca20.df)
RFoptions(coordnames=c("east", "north"), varnames="data") 

## covariance model with variance, scale and nugget to be estimated;
## just to abbreviate later on
M <- RMexp(var=NA, scale=NA) + RMnugget(var=NA)
\dontshow{if (!interactive()) M <- RMexp(var=NA, scale=NA)} % ok

## short definition of a trend using the fact that ca20.df is a
## data.frame
ca20.RFmod02 <- ~ 1 + altitude + M
(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df))


## long definition, which allows also for more general constructions
ca20.RFmod02 <- NA + NA*RMcovariate(ca20.df$altitude) + M 
(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df))


## Note that the following also works.
## Here, the covariance model must be the first summand
ca20.RFmod02 <- M + NA + ca20.df$altitude 
print(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df))

### The following does NOT work, as R assumes (NA + ca20.df$altitude) + M
(ca20.RFmod02 <- NA + ca20.df$altitude + M) 
try(ca20.fit02.RF <- RFfit(ca20.RFmod02, data=ca20.df)) ### error ...


## factors:
ca20.RFmod03 <- ~ 1 + area + M ### 
(ca20.fit03.RF <- RFfit(ca20.RFmod03, data=ca20.df))

\dontshow{FinalizeExample(); RFoptions(coordnames=NA, varnames=NA, modus="normal")}

}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
}
\keyword{spatial}
