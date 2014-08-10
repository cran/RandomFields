\name{RFfractaldim}
\alias{RFfractaldim}
\title{RFfractaldimension}
\description{
 The function estimates the fractal dimension of a process
}
\usage{
RFfractaldim(x, y = NULL, z = NULL, data, grid, 
 bin= seq(min(ct$x[2, ]) / 2, 
          min(ct$x[2,] * ct$x[3,] / 4, vario.n * min(ct$x[2,]) + 1),
          min(ct$x[2,])),
 vario.n=5,
 sort=TRUE,
 fft.m = c(65, 86), ## in % of range of l.lambda
 fft.max.length=Inf,%
 fft.max.regr=150000,
 fft.shift = 50, # in %; 50:WOSA; 100: no overlapping
 method=c("variogram", "fft"), % "box", "range", 
 mode=c("plot", "interactive"),
 pch=16, cex=0.2, cex.main=0.85,
 printlevel = RFoptions()$general$printlevel,
 height=3.5,
 ...)
}
\arguments{
 \item{x}{matrix of coordinates, or vector of x coordinates; if
 \code{x} is not given a grid with unit grid length is assumed}
 \item{y}{vector of y coordinates}
 \item{z}{vector of z coordinates}
% \item{T}{vector of time coordinates. 
% \code{T} must always be given in the \code{gridtriple} format,
% independently how the spatial part is defined.
% }
 \item{data}{the values measured.}
 \item{grid}{determines whether the vectors \code{x},
 \code{y}, and \code{z} should be
 interpreted as a grid definition, see Details. \code{grid}
 does not apply for \code{T}.}
 \item{bin}{sequence of bin boundaries for the empirical variogram}
 \item{vario.n}{first \code{vario.n} value of the empirical variogram
 are used for the regression fit that are not \code{NA}.}
 \item{sort}{If \code{TRUE} then the coordinates are permuted
 such that the largest grid length is in \code{x}-direction; this is
 of interest for algorithms that slice higher dimensional fields
 into one-dimensional sections.}
% \item{box.sequ}{sequence of box length for which the number of
% covering boxes is calculated.}
% \item{box.enlarge.y}{The grid of coordinates is always standardised to
% the marginal interval \eqn{[0, 1]}. The \code{data} are standardised
% to the interval \eqn{[0, \code{box.enlarge.y}]}.}
% \item{range.sequ}{sequence of box lengths; for each box length and a
% covering set of boxes the sum of the ranges for the values in each
% box is calculated.}
 \item{fft.m}{numeric vector of two components; interval of frequencies
 for which the regression should be calculated; the interval is given
 in percent of the range of the frequencies in log scale.}
 \item{fft.max.length}{The first dimension of the data is cut into pieces
 of length \code{fft.max.length}. For each piece the FFT is
 calculated and then the average for all pieces is taken. The pieces
 may overlap, see the argument \code{fft.shift}.}
 \item{fft.max.regr}{If the \code{fft.m} is too large, parts of the
 regression fit will take a very long time.
 Therefore, the regression fit is calculated only if the number points
 given by \code{fft.m} is less than \code{fft.max.regr}.
 }
 \item{fft.shift}{This argument is given in percent [of
 \code{fft.max.length}] and defines the overlap of the pieces defined
 by \code{fft.max.length}. If \code{fft.shift=50} the WOSA estimator is
 given; if \code{fft.shift=100} no overlap exist.}
 \item{method}{
 list of implemented methods to calculate the fractal dimension; see Details
 }
 \item{mode}{character. A vector with components
 \code{'nographics'}, \code{'plot'}, or \code{'interactive'}: 
 \describe{
 \item{\code{'nographics'}}{no graphical output}
 \item{\code{'plot'}}{the regression line is plotted}
 \item{\code{'interactive'}}{the regression domain can be chosen interactively}
 }
 Usually only one mode is given. Two modes may make sense
 in the combination \code{c("plot", "interactive")} in which case all the
 results are plotted first, and then the interactive mode is called. 
 In the interactive mode, the regression domain is chosen by
 two mouse clicks with the left
 mouse; a right mouse click leaves the plot.
 }
 \item{pch}{vector or scalar; sign by which data are plotted.}
 \item{cex}{vector or scalar; size of \code{pch}.}
 \item{cex.main}{The size of the title in the regression plots.}
 \item{printlevel}{integer. If \code{printlevel} is 0 nothing is
 printed. If \code{printlevel=1} error messages are printed. 
 If \code{printlevel=2} warnings and the regression results
 are given. If \code{printlevel>2} tracing information is given.
 }
 \item{height}{height of the grahics window}
 \item{...}{graphical arguments}
}
\details{ 
 The function calculates the fractal dimension by various methods:
 \itemize{
 \item variogram method
% \item box counting
% \item min / max method
 \item Fourier transform
 }
}
\value{
 The function returns a list with elements
 \code{vario}, %\code{box}, \code{range},
 \code{fft} corresponding to
 the 2 methods given in the Details.

 Each of the elements is itself a list that contains the
 following elements.
 \item{x}{the x-coordinates used for the regression fit}
 \item{y}{the y-coordinates used for the regression fit}
 \item{regr}{the return list of the \command{\link[stats]{lm}}.}
 \item{sm}{smoothed curve through the (x,y) points}
 \item{x.u}{\code{NULL} or the restricted x-coordinates given
 by the user in the interactive plot}
 \item{y.u}{\code{NULL} or y-coordinates according to \code{x.u}}
 \item{regr.u}{\code{NULL} or the return list of
 \command{\link[stats]{lm}} for \code{x.u} and \code{y.u}}
 \item{D}{the fractal dimension}
 \item{D.u}{\code{NULL} or the fractal dimension corresponding to the
 user's regression line}
 }
\references{
% Overviews:
% \itemize{
% \item{-}{-}
% }

 variogram method
 \itemize{
 \item Constantine, A.G. and Hall, P. (1994)
 Characterizing surface smoothness via estimation of effective
 fractal dimension. \emph{J. R. Statist. Soc. Ser. B} \bold{56}, 97-113.
 }

% box counting
% \itemize{
% \item{-}{-}
% }

% min/max
% \itemize{
% \item{-}{-}
% }

 fft
 \itemize{
 \item Chan, Hall and Poskitt (1995)
 }


}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de/de/publications/software}}
\seealso{
 \command{\link{RMmodel}}, \command{\link{RFhurst}}
}

\keyword{ spatial }%-- one or more ...

\examples{
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again
x <- seq(0, 10, if (interactive()) 0.001 else 1)
z <- RFsimulate(RMexp(), x)
if (interactive()) str(RFfractaldim(data=z))
\dontshow{FinalizeExample()}
}