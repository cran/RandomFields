\name{Dev}
\alias{Dev}

\title{Choosing the device}
\description{
  \code{Dev} chooses between the graphical devices
  screen, postscript and pdf
}
\usage{
Dev(on, dev, ps=NULL, cur.cex=TRUE, paper="special", width=5, height=5,
    quiet=FALSE, ...)
}
\arguments{
  \item{on}{logical. Indicates whether \code{dev} should be switched on
    or off}
  \item{dev}{see Details}
  \item{ps}{name of the pdf or postscript file}
  \item{cur.cex}{logical.  If \code{TRUE} the par parameters of the current
    device are used, not the standard parameters}
  \item{paper}{kind of paper. Postscript parameter}
  \item{width}{width of figure. Postscript and pdf parameter}
  \item{height}{height of figure. Postscript and pdf parameter}
  \item{quiet}{logical.  If \code{FALSE} additional information is given.}
  \item{...}{further parameters for pdf or postscript}
}
\details{
  The parameter \code{dev} might be
  \describe{
    \item{logical}{If suffix of \code{ps} is either "eps", "ps" or "pdf"
      the respective file is created. Otherwise, a postscript file ending
      with suffix "eps" is created if \code{dev=TRUE} and a pdf file with
      suffix "pdf" if \code{dev=FALSE}}
    \item{character}{A function with name \code{dev} is called, and the
      suffix \code{dev} is added to \code{ps}}
    \item{numeric}{\code{Dev} switches to the device with number \code{dev};
      if such a device does not exist, a new X11 device is created.
      If \code{dev<2} the device opened by \code{Dev(on=TRUE,...)} is
      closed by \code{Dev(on=FALSE)}; otherwise only \code{par(new=F)}
      is applied to the device.
    }
    The parameters \code{dev}, \code{ps}, \code{cur.cex}, \code{paper},
    \code{width}, \code{height} and \code{...} are ignored if \code{on=FALSE}.
  }
}
\value{
  \code{NULL}. Side effect is that the global variable \code{.dev.orig}
  is created when \code{Dev(on=TRUE,...)} is called.
}

\author{Martin Schlather, \email{martin.schlather@cu.lu}
  \url{http://www.cu.lu/~schlathe}}

\examples{
  ## first an eps-file test.eps is created, then a jpeg-file,
  ## finally the figure is plotted on the screen
  dir(pattern="test*")
  for (dev in list(TRUE, "jpeg", 1)) {
    print(dev)
    size <- if (dev=="jpeg") 450 else 5
    Dev(TRUE, dev, ps="test", height=size, width=size)
    plot(0, 0, main=paste("dev=", dev))
    # readline("press return")
    Dev(FALSE)
  }
  dir(pattern="test*")
}

\keyword{device}
\keyword{utilities}









