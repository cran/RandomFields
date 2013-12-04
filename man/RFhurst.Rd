\name{RFhurst}
\alias{RFhurst}
%\alias{Hurst}
%\alias{long memory dependence}
\title{Hurst coefficient}
\description{
The function estimates the Hurst coefficient of a process
}
\usage{
RFhurst(x, y = NULL, z = NULL, data, sort = TRUE,
 block.sequ = unique(round(exp(seq(log(min(3000, dimen[1]/5)),
 log(dimen[1]),
 len = min(100, dimen[1]))))),
 fft.m = c(1, min(1000, (fft.len - 1)/10)),
 fft.max.length = Inf, method = c("dfa", "fft", "var"),
 mode = c("plot", "interactive"),
 pch = 16, cex = 0.2, cex.main = 0.85,
 printlevel = RFoptions()$general$printlevel, height = 3.5,
 ...)
}
\arguments{
 \item{x}{matrix of coordinates, or vector of x coordinates}
 \item{y}{vector of y coordinates}
 \item{z}{vector of z coordinates}
 \item{data}{the data}
% \item{grid}{logical; determines whether the vectors \code{x},
% \code{y}, and \code{z} should be
% interpreted as a grid definition, see Details. \code{grid}
% does not apply for \code{T}.}
 \item{sort}{logical. If \code{TRUE} then the coordinates are permuted
 such that the largest grid length is in \code{x}-direction; this is
 of interest for algorithms that slice higher dimensional fields
 into one-dimensional sections.
 }
 \item{block.sequ}{ascending sequences of block lengths for which the
 detrended fluctuation analysis and the variance method is
 performed.}
 \item{fft.m}{vector of 2 integers; lower and upper endpoint of
 indices for the frequency which are used in the calculation of the
 regression line for the periodogram near the origin.}
 \item{fft.max.length}{if the number of points in \code{x}-direction is
 larger than \code{fft.max.length} then the segments of length
 \code{fft.max.length} are considered, shifted by
 \code{fft.max.length/2} (WOSA-estimator).}
 \item{method}{list of implemented methods to calculate the Hurst parameter; see Details}
 \item{mode}{character. A vector with components
 \code{'nographics'}, \code{'plot'}, or \code{'interactive'}: %; see the Details.
 \describe{
 \item{\code{'nographics'}}{no graphical output}
 \item{\code{'plot'}}{the regression line is plotted}
 \item{\code{'interactive'}}{the regression domain can be chosen
 interactively}
 }
 Usually only one mode is given. Two modes may make sense
 in the combination c("plot", "interactive") in which case all the
 results are plotted first, and then the interactive mode is called. 
 In the interactive mode, the regression domain is chosen by
 two mouse clicks with the left
 mouse; a right mouse click leaves the plot.
 }
 \item{pch}{vector or scalar; sign by which data are plotted.}
 \item{cex}{vector or scalar; size of \code{pch}.}
 \item{cex.main}{font size for title in regression plot;
 only used if mode includes \code{'plot'}
 or \code{'interactive'}}
 \item{printlevel}{integer. If \code{printlevel} is 0 or 1
 nothing is printed. 
 If \code{printlevel=2} warnings and the regression results
 are given. If \code{printlevel>2} tracing information is given.
 }
 \item{height}{height of the graphics window}
 \item{...}{graphical parameters}
}
\details{
 The function is still in development. Several functionalities do
 not exist - see the code itself for the current stage.

 The function calculates the Hurst coefficient by various methods:
 \itemize{
 \item detrended fluctuation analysis (dfa)
 \item aggregated variation (var)
 \item periodogram or WOSA estimator (fft)
 }
}
\value{
 The function returns a list with elements
 \code{dfa}, \code{varmeth}, \code{fft} corresponding to
 the three methods given in the Details.
 
 Each of the elements is itself a list that contains the
 following elements.
 
 \item{x}{the x-coordinates used for the regression fit}
 
 \item{y}{the y-coordinates used for the regression fit}
 
 \item{regr}{the coefficients of the \command{\link[stats]{lm}}.}
 
 \item{sm}{smoothed curve through the (x,y) points}
 
 \item{x.u}{\code{NULL} or the restricted x-coordinates given
 by the user in the interactive plot}

 \item{y.u}{\code{NULL} or y-coordinates according to \code{x.u}}

 \item{regr.u}{\code{NULL} or the coefficients of 
 \command{\link[stats]{lm}} for \code{x.u} and \code{y.u}}

 \item{H}{the Hurst coefficient}

 \item{H.u}{\code{NULL} or the Hurst coefficient corresponding to the
 user's regression line}
}
\references{
% Overviews:
% \itemize{
% \item 
% }

 detrended fluctuation analysis
 \itemize{
 \item Peng, C.K., Buldyrev, S.V., Havlin, S., Simons, M., Stanley,
 H.E. and Goldberger, A.L. (1994)
 Mosaic organization of DNA nucleotides
 \emph{Phys. Rev. E} \bold{49}, 1685-1689
 }

 aggregated variation
 \itemize{
 \item Taqqu, M.S. and Teverovsky, V. (1998)
 On estimating the intensity of long range dependence in finite and
 infinite variance time series. In: Adler, R.J., Feldman, R.E., and
 Taqqu, M.S. \emph{A Practical Guide to Heavy Tails, Statistical
 Techniques an Applications.} Boston: Birkhaeuser

 \item
 Taqqu, M.S. and Teverovsky, V. and Willinger, W. (1995)
 Estimators for long-range dependence: an empirical study.
 \emph{Fractals} \bold{3}, 785-798
 }

 periodogram
 \itemize{
 \item Percival, D.B. and Walden, A.T. (1993)
 \emph{Spectral Analysis for Physical Applications: Multitaper and
 Conventional Univariate Techniques}, Cambridge: Cambridge
 University Press.
 
 \item Welch, P.D. (1967) The use of {F}ast {F}ourier {T}ransform for
 the estimation of power spectra: a method based on time averaging
 over short, modified periodograms \emph{IEEE Trans. Audio
	Electroacoustics} \bold{15}, 70-73. 
 }
}
\author{Martin Schlather, \email{schlather@math.uni-mannheim.de}
 \url{http://ms.math.uni-mannheim.de}}
\seealso{
 \command{\link{RMmodel}}, \command{\link{RFfractaldim}}
}
\examples{
set.seed(0)
\dontshow{.randomfields.options = options()$warn; options(warn=0)}
x <- runif(1000)
if (interactive()) {
 h <- RFhurst(1:length(x), data=x)
} else {
 h <- RFhurst(1:length(x), data=x, mode = "nographics")
}
\dontshow{options(warn = .randomfields.options)}

}
\keyword{ spatial }%-- one or more ...






% LocalWords: hurst gridtriple sequ exp len fft Inf dfa var pch cex WOSA regr
% LocalWords: printlevel RFoptions periodogram nographics itemize varmeth
% LocalWords: lm sm Schlather url ealso RMmodel
