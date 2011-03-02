%%% stimmt die Beschreibung noch??!!

\name{eval.parameters}
\alias{eval.parameters}

\title{Interactive menu}
\description{
  eval.parameters provides an interactive menu on a X11 graphical
  device of R
}
\usage{
eval.parameters(variable, entry, update, simulate, dev, create = TRUE, 
		col.rect = "red", col.bg = "blue", col.sep = "grey", 
                col.left = "red", col.mid = "white", col.right = "white", 
                col.line = "red", col.txt = "black",
                cex=0.8, cex.i=cex, sep="-----", ...)
}
\arguments{
  \item{variable}{string.  The name of the variable to be changed.
    The name consist of \code{$} and  \code{[[]]} expression pointing
    to sublists of a list. The complete list must be given by
    name in \code{\dots}}
  \item{entry}{A list of lists.  See Details.}
  \item{update}{logical.  If \code{TRUE} then \code{simulate} is called
    after each interactive input.}
  \item{simulate}{function that is called if simulations are to
    updated. The parameters must equal the variables given by
    \code{...}; the function must return the complete list indicated by
    \code{variable}.}
  \item{dev}{Before calling \code{eval.parameters()} \code{split.screen}
    must have been called.  \code{dev} gives the screen on which the
    interactive menu should be plotted.}
  \item{create}{logical.  If \code{TRUE} missing list elements of
    \code{variable} are created.}
  \item{col.rect}{colour of the button for free input.}
  \item{col.bg}{colour of a interactive bar}
  \item{col.sep}{colour of the separating line}
  \item{col.left}{colour of preceding element}
  \item{col.mid}{colour for the message}
  \item{col.right}{colour of subsequent element}
  \item{col.line}{colour of the marking line in interactive bars of
    absolute choice.}
  \item{col.txt}{colour of headers}
  \item{cex}{font height of headers}
  \item{cex.i}{font height of text for elements}
  \item{sep}{style of added characters in separating line.}
  \item{\dots}{The input variables given by name;
    the names may not start with a dot;
    There the complete list to which \code{variable} refers must be
    given. Further additional parameters for the function \code{simulate}.
  }
}
\note{
  To the list given by \code{variable} the element \code{.history}
  is added.  \code{.history} is a list that contains the history
  of the user input. Each element is a list where the first entry
  is the number of the menu, the second and the third entries are the
  former and the new value. Exception: for entries with character
  \code{val}, the value of \code{val} is returned as second entry.
  Consequently, the name \code{.history} should not be used for other
  purposes in \code{variable}.

  Further, any variable name given in \code{...} must start with a letter.
}

  \details{
  \code{eval.parameters} shows a menu list on X11.  Depending on the
  mode of the variables the menu bars have a different appearance and
  behave differently if the user clicks on the bar.
  Most of the menu bars have a small rectangle on
  the right hand side.  If this rectangle is pressed the input of a
  variable is expected in the xterm where R is run.

  
  \code{entry} is a list of lists.  Each list may contain the following
  elements:
  \itemize{
    \item \code{name} : header for menu button if \code{var} is not
    \code{NULL}; otherwise printed instead of a menu button
 
    \item \code{var} :
    \itemize{
      \item \code{NULL} : if \code{val=NULL} then
      it is a separating line in colour \code{col.sep};
      \code{name} is surrounded by \code{sep}; all other elements of
      the list are ignored. 
      If \code{val} is a string then \code{val} is interpreted as a
      function; a special string is "simulate", which entails the call of
      the function \code{simulate} with the appropriate parameters.
      \item string : selected element of the list that is given
      by \code{variable}. \emph{A special string "undo" will be
	installed to undo things}.
    }
    \item \code{val} :
    \itemize{
      \item \code{function(d, v)} gives the update for
      \code{var}.  If \code{v} is missing, a starting value (for
      d=1/2) is expected.  Otherwise, \code{v} is the current value
      of \code{var} and \code{d} is the choice of the user,
      a value in \eqn{[0,1]}
      \item \code{TRUE} : logical bar. 
      \item \code{FALSE} : logical bar.  The value for \code{var} is
              negated before shown.
      In the menu, the negative value for \code{var} is shown.
      \item \code{NULL} : a string is read from the terminal inot \code{var}.
      \item character (vector): if \code{var} is given then
      this vector of strings is interpreted as belonging to a 
      categorical variable
      \eqn{1,\ldots,\code{length(val)}}{1,...,\code{length(val)}}
      and \code{var} gives the number of the selected elements. 
      If \code{var=NULL} then \code{val} is interpreted as a function; a
      special string is "simulate", which entails the call of 
      the function \code{simulate} with the appropriate parameters.
    }
    \item \code{delta} : logical.
    In the menu bar absolute values are plotted if \code{delta=FALSE} 
    and increments otherwise.  Only considered if \code{val} is a function.
    \item \code{cond} : The menue points is shown only if the given
    condition is satisfied. The condition must be expressed by named elements
    of the list variable, see example.
    \item \code{col} : colour that overwrites the standard colour for the
    rectangle or the separating text.
    \item \code{update} : logical. If not missing, its value overwrites
    locally the value of the global parameter \code{update}.
    \item \code{...} : additional parameters for \code{simulate} that overwrite
    the values given in \code{...} in the call of \code{eval.parameters}.
  }
}

\value{
  The first variable given in \dQuote{\code{...}}, which is a list.
  To this list the entry \code{.history} is added.
  
  If the users enters \sQuote{exit immediately} at any point, the program
  stops with an error message.
}

\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}

%\seealso{\command{\link[SoPhy]{define.horizons}}, \command{\link[SoPhy]{SoPhy}}}
\seealso{\command{\link{useraction}}}

\examples{
\dontrun{
  ## the following lines define a menu that does not make
  ## too much sense, but shows the various kinds of buttons

  quadratic <- function(d, v, a, mini=0, maxi=Inf) {
    d <- pmin(1, pmax(0, d)) - 0.5
    d <- ((d>0) * 2 - 1) * d^2 * a * 4
    if (missing(v)) d else pmax(mini, pmin(maxi, v + d))
  }

  simulate <- function(H, par) { ## not a serious example
    Print(c(H$x$var, par, runif(1)))
    return(H)  ## the function must return the first parameter
  }

  entry <- list(
     list(name="Nonsense Menu"),
     list(name="Simulate!", val="simulate", col="blue"),
     list(name="show H", val="str(H)", col="blue"),
     list(name="colx", var="colour",
          val=c("red", "green", "blue", "brown")),
     list(name="open", var="closed", val=FALSE, par=4.5),
     list(name="modifying", var="modify", val=TRUE, par=5),
     list(name="probability", var="probab", delta=FALSE,
          val=function(d, v) pmin(1, pmax(0, d))),
     list(name="variance", var="var", delta=TRUE,
          val=function(d, v) quadratic(d, v, 10)),
     list(name="name", var="name", par=3, cond="modify")
    )

  scr <- split.screen(rbind(c(0, 0.45, 0, 1), c(0.5, 1, 0, 1)))
  ## before proceeding make sure that both the screen and the xterm
  ## are completely visible

  H <- list(modify=5, x=list()) # note that in this example eval.parameters 
  ##          will be called by by H$x, hence modify=5 will be left
  ##          unchanged.
  options(locatorBell=FALSE)
%   source("/home/schlather/R/RF/RandomFields/R/auxiliary.R"); source("/home/schlather/R/RF/RandomFields/R/evalpar.R");
  useraction("start.register")  ## registring the user's input
  Print(eval.parameters("H$x", entry, simulate, update=TRUE, dev=scr[2],
     H=H, par=17)) # do not forget to call by name
  getactions()

  ## replay the user's input
  useraction("replay")
  Print(eval.parameters("H$x", entry, simulate, update=TRUE, dev=scr[2],
     H=H, par=17))
}
}
\keyword{spatial}







%  LocalWords:  val cond dQuote SoPhy
