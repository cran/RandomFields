
## sudo apt-get install tcl8.5-dev
## sudo apt-get install tk8.5-dev
## see R Installation and Administration
##./configure --with-tcl-config=/usr/lib/tcl8.5/tclConfig.sh --with-tk-config=/usr/lib/tk8.5/tkConfig.sh
# sudo tar xvf ~/TMP/bwidget-1.9.5.tar 


# readline <- function(...) return("") 

.onAttach <- function (lib, pkg) {
  #packageStartupMessage("\n")
}

#.onUnload <- function(lib, pkg){
#  }
#Implementierung von Cox & Isham's non-separable model

#individuelle Praeferenzliste:
#  1) erste praeferenz
#  2) if # pkte < 10 / 50 -> Gauss??
#  4) if nugget/variance > 0.01 -> CE (and not spectral, e.g.)
#  3) if C(maxh)/C(0) > 0.05  -> TBM else CE


.onUnload <- function(lib, pkg){
  RFoptions(storing=FALSE) ## delete everything
}
 
