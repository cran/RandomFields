\name{regression}
\alias{regression}
\title{Regression plot}
\description{
Interactive regression plot
}
\usage{
regression(x, y, main, scr, value=function(regr) regr$coeff[2], 
           variable="var", mode=c("nographics", "plot", "interactive"),
           regr.select=c("coefficients"), averaging=FALSE, cex.main=0.85,
           col.passive = "grey", col.active = "green",
           col.main.passive = "black", col.main.active = "red",
           col.points.chosen = "lightblue", col.abline = "yellow",
           col.abline.chosen = "darkblue", col.smooth = "red", ...)
}
\arguments{
  \item{x}{x-coordinates of the regression plot}
  \item{y}{y-coordinates}
  \item{main}{main title}
  \item{scr}{\code{split.screen} must have been called beforehand;
    \code{scr} gives the screen number}
  \item{value}{function; the characteristic of interest is a
    function of the regression output.}
  \item{variable}{name of the characteristic of interest, plotted in the graphs}
  \item{mode}{A scalar character of one of the following values
    'silent', 'plot', or 'interactive': %; see the Details.
    \describe{
      \item{'silent':}{no output except the return of the result}
      \item{'plot':}{the regression line is plotted}
      \item{'interactive':}{the regression domain can be chosen
	interactively}
    }
    Any mode includes all the functionalities of the modes that precede in
    the above enumeration. 
    In the interactive mode, the regression domain is chosen by
    two mouse clicks with the left
    mouse; a right mouse click leaves the plot.
  }
  \item{regr.select}{vector of characters that select elements from the
    \command{\link[stats]{lm}} return list; these elements are returned by
    \code{regression}.}
  \item{averaging}{logical. If \code{TRUE} and \code{y} is matrix then
    the values of \code{y} are averaged row-wise beforehand}
  \item{cex.main}{font size of the title}
  \item{col.passive}{colour of the data points if plot is not active}
  \item{col.active}{colour of the data points if plot is active}
  \item{col.main.passive}{colour of the title if plot is not active}
  \item{col.main.active}{colour of the title if plot is active}
  \item{col.points.chosen}{colour of the selected points}
  \item{col.abline}{colour of the regression line for all points}
  \item{col.abline.chosen}{colour of the regression line for the selected points}
  \item{col.smooth}{colour for the smoothed data}
  \item{...}{further graphical parameters}
}
%\details{}
\value{
  The function returns a list with the following components
  \item{\code{regr}}{part of the return list of \command{\link[stats]{lm}}}
  \item{\code{val}}{the calculated value of the variable of interest}
  \item{\code{regr.u}}{\code{NULL} or part of the return list of
    \command{\link[stats]{lm}} for \code{x.u} and \code{y.u}}
  \item{\code{val.u}}{\code{NULL} or the calculated value of the
    variable of interest for the restricted domain}
  \item{\code{x.u}}{\code{NULL} or the restricted x-coordinates given
    by the user in the interactive plot}
  \item{\code{y.u}}{\code{NULL} or y-coordinates according to \code{x.u}} 
  \item{sm}{smoothed curve through the (x,y) points}
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/institute}}
%\seealso{}
%\examples{}
\keyword{ spatial }%-- one or more ...

