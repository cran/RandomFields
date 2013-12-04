
print.RFopt <- function(x, ..., max.level = 3) {
  str(x, max.level = max.level) #
}




## definition of class 'RFempVariog'
#setClass("RFopt", representation(opt="ANY"))

## rules for validity checking of 'RFempVariog' objects
#setValidity("RFopt", 
#            function(object){
#              if (!is.list(object@opt) || !all(sapply(object@opt, is.list))) {
#                Print(object@opt, is.list(object@opt))
#                return("object must be a list of lists")
#              }
#              return(TRUE)
#            })
#
#setMethod(f="show", signature="RFopt", #, useSource="missing"),
#          definition=function(object) {str(object@opt)})
