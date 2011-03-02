\name{useraction}
\alias{useraction}
\alias{putactions}
\title{Set input behaviour}
\description{
  The functions store values that influence the behaviour
  of \command{\link{Locator}} and \command{\link{Readline}}
}
\usage{
useraction(action=c("none", "start.register", "continue.register",
                    "replay", "endless", "delete"),
		    sleep, wait, PrintLevel, actionlist)
putactions(l)
}
\arguments{
  \item{action}{string of the following values
    \describe{
      \item{\sQuote{"none"}}{the behaviour of
	\command{\link{Readline}} and \command{\link{Locator}}
	is that of \command{\link[base]{readline}} and
	\command{\link[graphics]{locator}}, respectively.
	The values are \bold{not} stored.
      } 
      \item{\sQuote{"start.register"}}{the behaviour of
	\command{\link{Readline}} and \command{\link{Locator}}
	is that of \command{\link[base]{readline}} and
	\command{\link[graphics]{locator}}, respectively.
	The current list of stored values is
	deleted and new inputs by \command{\link{Readline}} and
	\command{\link{Locator}} are \bold{stored}.}
      \item{\sQuote{"continue.register"}}{the behaviour of
	\command{\link{Readline}} and \command{\link{Locator}}
	is that of \command{\link[base]{readline}} and
	\command{\link[graphics]{locator}}, respectively.
	The current list of stored values is continued with
	inputs from \command{\link{Readline}} and
	\command{\link{Locator}}.}
      \item{\sQuote{"replay"}}{
	\command{\link{Locator}} and \command{\link{Readline}}
	do not get the data from the mouse or terminal, but
	read them from the stored list.
	If the end of storage list is reached, the user is informed and
	the \code{action} is switched to \dQuote{none}.
      }
      \item{\sQuote{"endless"}}{similar to \sQuote{"replay"}, but
	storage list is looped.
      }
      \item{\sQuote{"delete"}}{
	The current list of stored values is
	deleted	and \code{action} is set to \dQuote{none}
      }
    }
  }
  \item{actionlist}{list of stored values by \command{\link{Readline}} and
    \command{\link{Locator}}} 
  \item{sleep}{sleeping time before the blinking starts (replay mode of
    \command{\link{Locator}}); value is kept until explicit changing;
    starting value: 2.5
  } 
  \item{wait}{regulated the blinking speed of the characters
    (\command{\link{Locator}}) and of the 
    speed of appearance of characters (\command{\link{Readline}}); value
    is kept until explicit changing; 
    starting value: 0.1}
  \item{PrintLevel}{
    if \code{Printlevel>1} and \command{\link{Readline}} or
    \command{\link{Locator}} is in replay mode then
    the current pointer position within the storage list
    and the tracing information (if any) is shown. The value
    is kept until explicit changing; 
    starting value: 0
  }
  \item{l}{list of all parameters relevant to the behaviour of
    \command{\link{Readline}} and 
    \command{\link{Locator}}.
  }
}
\value{
  \command{useraction} returns \code{NULL}.

  \command{putactions} returns invisibly the state before the new state is set.
}
\details{
  The function \command{putactions} restores the state
  taken by \command{\link{getactions}}.
  
  Concerning the function \command{useraction}, the parameter
  \code{action} must be given. If any of the other 
  parameters is not given, its current value is kept.
  If \code{action="start.register"} or \code{action="delete"},
  the value of \code{actionlist} is ignored.
  
  If replaying, \command{\link{Readline}} shows the input
  values on the screen and \command{\link{Locator}} shows
  blinking characters as side effects.
}
\seealso{
  \command{\link{eval.parameters}},
  \command{\link{getactions}},
  \command{\link{Locator}},
  \command{\link{Readline}},
 % \command{\link{userinput}},
}
\author{Martin Schlather, \email{martin.schlather@math.uni-goettingen.de}
  \url{http://www.stochastik.math.uni-goettingen.de/~schlather}}
\examples{
%   library(RandomFields, lib="~/TMP"); source("/home/schlather/R/RF/RandomFields/R/auxiliary.R")
sample.interaction <- function(read.r2=TRUE, type="p", pch=16) {
  r1 <- Readline(prompt="first input: ")
  r2 <- if (read.r2) Readline(prompt="second input: ", info="2nd input")
  else " -- no input  -- "
  cat("Mark some points with the right mouse key, then leave with the left mouse key\n")
  do.call(getOption("device"), list())
  plot(Inf, Inf, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="")
  l <- Locator(100, info="locator input", type=type, pch=pch)
  r3 <- Readline(prompt="third input: ", info="last input")
  return(list(r1=r1, r2=r2, l=l, r3=r3))
}

%   source("/home/schlather/R/RF/RandomFields/R/auxiliary.R")

\dontrun{

#########################
## user terminal input ##
#########################
{
  useraction("start.register")
  str(sample.interaction())
  str(l <- getactionlist())
}

#################
## just replay ##
#################
%   source("/home/schlather/R/RF/RandomFields/R/auxiliary.R")
useraction("replay", sleep=2, wait=interactive() * 0.2 * 5,
           PrintLevel=2, actionlist=l)
str(sample.interaction(type="l"))


##############################################
## modify first the input list, then replay ##
##############################################
l2 <- l[-2]
l2[[1]] <- list("some other words", info="changed input")
str(l2)
useraction("replay", sleep=1, wait=interactive() * 0.05,
           PrintLevel=0, l2) 
str(sample.interaction(read.r2=FALSE, type="o")) # input now from l2
Readline(prompt="new input: ", info="?!") # switch to terminal
##                                                since end of stored list
# str(getactionlist())    # new input has not been not stored ...


####################################################
## use of the two lists, l and l2, in a mixed way ##
####################################################
useraction("replay", sleep=2, wait=interactive() * 0.05,
           PrintLevel=0, actionlist=l)
Readline(prompt="first input of 1: ")
dump <- getactions()
useraction("replay", sleep=0.5, wait=interactive() * 0.05,
           PrintLevel=0, actionlist=l2)
Readline(prompt="first input of 2: ")
dump2 <- putactions(dump)
Readline(prompt="second input of 1: ")
## locator call reading from list 1:
plot(Inf, Inf, xlim=c(-1,1), ylim=c(-1,1), xlab="", ylab="")
Locator(100, info="locator input", type="p")
dump <- putactions(dump2)
##locator call reading from list 2:
Locator(100, info="locator input", type="p", pch=20, col="blue")
Readline(prompt="last input of 2: ")
putactions(dump)
Readline(prompt="last input of 1: ")


#####################################################
## see help("eval.parameters") for another example ##
#####################################################
}

}
\keyword{utilities}

