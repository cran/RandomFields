\name{Changings}
\alias{Changings}
\alias{changings}

\alias{tbmdim}
\alias{coord_units}
\alias{new_coord_units}
\alias{variab_units}

\title{Documentation of some further changings}


\description{
  \itemize{
    \item{Version 3.3}
     \itemize{
       \item\code{RFempiricalvariogram}, \code{RFempiricalcovariance}
       and \code{RFempiricalmadogram} became obsolete. Use 
       \command{\link{RFvariogram}},
       \command{\link{RFcov}},
       \command{\link{RFmadogram}}
       instead.
       \item \code{RFoptions(grDefault=FALSE)} returns to the old style
       of graphical device handling. Otherwise there is no handling.
       \item C code is started to be parallelized.
       \item Some new \link{Multivariate RMmodels}
       \item New way of passing models, see \link{RFformula}, which
       allows connections (formulae) between parameters, e.g. one
       parameter value might be twice as large as another parameter
       value. Also dummy variables can be \command{\link{RMdeclare}}d.
     }
       
    \item{Options getting obsolete (Version 3 and older)}
    \itemize{
      \item\code{oldstyle} is becoming \code{warn_oldstyle} 
      \item\code{newstyle} is becoming \code{warn_newstyle} 
      \item\code{newAniso} is becoming \code{warn_newAniso} 
      \item\code{ambiguous} is becoming \code{warn_ambiguous} 
      \item\code{normal_mode} is becoming \code{warn_normal_mode} 
      \item\code{colour_palette} is becoming \code{warn_colour_palette} 
    }
    
    
    \item \bold{Changings in option names}
    \itemize{
      \item several changes in \code{RFoptions()$graphics} in version 3.1.11
      \item\code{pdfnumber} became in version 3.0.42 \code{filenumber}      
      \item\code{pdfonefile} became in version 3.0.42 \code{onefile}      
      \item\code{pdffile} became in version 3.0.42 \code{file}      
      \item\code{tbmdim} became in version 3.0.41 \code{reduceddim} 
      \item\code{coord_units} became in version 3.0.39 \code{coordunits} 
      \item\code{new_coord_units} became in version 3.0.39
      \code{new_coordunits} 
      \item\code{variab_units} became in version 3.0.39 \code{varunits}      
    }
  } 
}

\seealso{
  \link{MajorRevisions},
  \link{RandomFields}.
}

\me

\examples{\dontshow{StartExample()}
## no examples given
\dontshow{FinalizeExample()}}
