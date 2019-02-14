# This file has been created automatically by 'rfGenerateMaths'
setMethod("asin", signature = CLASS_CLIST, definition=function(x) R.asin(x))
setMethod("atan", signature = CLASS_CLIST, definition=function(x) R.atan(x))
setMethod("atan2", signature = c(CLASS_CLIST,'ANY'), definition=function(y, x) R.atan2(y, x))
setMethod("atan2", signature = c('ANY',CLASS_CLIST), definition=function(y, x) R.atan2(y, x))
setMethod("cos", signature = CLASS_CLIST, definition=function(x) R.cos(x))
setMethod("sin", signature = CLASS_CLIST, definition=function(x) R.sin(x))
setMethod("tan", signature = CLASS_CLIST, definition=function(x) R.tan(x))
setMethod("acosh", signature = CLASS_CLIST, definition=function(x) R.acosh(x))
setMethod("asinh", signature = CLASS_CLIST, definition=function(x) R.asinh(x))
setMethod("atanh", signature = CLASS_CLIST, definition=function(x) R.atanh(x))
setMethod("cosh", signature = CLASS_CLIST, definition=function(x) R.cosh(x))
setMethod("sinh", signature = CLASS_CLIST, definition=function(x) R.sinh(x))
setMethod("tanh", signature = CLASS_CLIST, definition=function(x) R.tanh(x))
setMethod("exp", signature = CLASS_CLIST, definition=function(x) R.exp(x))
setMethod("log", signature = CLASS_CLIST, definition=function(x) R.log(x))
setMethod("expm1", signature = CLASS_CLIST, definition=function(x) R.expm1(x))
setMethod("log1p", signature = CLASS_CLIST, definition=function(x) R.log1p(x))
exp2 <- R.exp2
setMethod("log2", signature = CLASS_CLIST, definition=function(x) R.log2(x))
setMethod("^", signature = c(CLASS_CLIST,'ANY'), definition=function(e1, e2) R.pow(e1, e2))
setMethod("^", signature = c('ANY',CLASS_CLIST), definition=function(e1, e2) R.pow(e1, e2))
setMethod("sqrt", signature = CLASS_CLIST, definition=function(x) R.sqrt(x))
hypot <- R.hypot
cbrt <- R.cbrt
setMethod("ceiling", signature = CLASS_CLIST, definition=function(x) R.ceil(x))
setMethod("abs", signature = CLASS_CLIST, definition=function(x) R.fabs(x))
setMethod("floor", signature = CLASS_CLIST, definition=function(x) R.floor(x))
setMethod("%%", signature = c(CLASS_CLIST,'ANY'), definition=function(e1, e2) R.fmod(e1, e2))
setMethod("%%", signature = c('ANY',CLASS_CLIST), definition=function(e1, e2) R.fmod(e1, e2))
setMethod("round", signature = c(CLASS_CLIST, 'missing'), definition=function(x, digits=0) R.round(x))
setMethod("trunc", signature = CLASS_CLIST, definition=function(x) R.trunc(x))
erf <- R.erf
erfc <- R.erfc
setMethod("lgamma", signature = CLASS_CLIST, definition=function(x) R.lgamma(x))


