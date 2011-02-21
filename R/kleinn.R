
"Forest" <- function(x, model, percent, ref.edge, coarse, gausspercent,
                      debug=FALSE) {
  len <- c(length(percent), length(ref.edge), length(coarse), 
           length(gausspercent))
 
  idx <- len > 1
  names <- list(paste("perc", percent, sep=""),
                paste("area", ref.edge, sep=""),
                paste("coar", coarse, sep=""),
                paste("gaus", gausspercent, sep=""))[idx]
  len<- len[idx]
  
  result <- double(prod(len))
  GaussRF(x, x, model=model, gridtriple=TRUE, n=0, Stor=TRUE,
          Print=1 + debug * 4, method="TBM3")
  nx <- length(seq(x[1], x[2], x[3]))
  
  if (debug) {    
    gausspercent <- gausspercent[length(gausspercent)]
    coarse <- coarse[length(coarse)]
    ref.edge <- ref.edge[length(ref.edge)]
    percent <- percent[length(percent)]

    gauss <- double(nx * nx)
    binary <- integer(nx * nx)
    decreased <- integer(nx * nx)
    refarea <- integer(nx * nx)
    nrdecr <- integer(1)
    ncdecr <- integer(1)
    areathreshold <- integer(1)
    
    if (nx < 100)
      Print("analyseForstImages",
       as.integer(0), # keynr
       as.double(gausspercent), as.integer(length(gausspercent)), # gausspercent
       as.integer(coarse), as.integer(length(coarse)), # coarse
       as.integer(ref.edge), as.integer(length(ref.edge)),# def of reference area
       as.double(percent), as.integer(length(percent)), # to be forest
       gauss,
       binary,
       decreased,
       nrdecr,
       ncdecr,
       refarea,
       result,
       DUP=FALSE, PACKAGE="RandomFields")
    
    .C("analyseForstImages",
       as.integer(0), # keynr
       as.double(gausspercent), as.integer(length(gausspercent)), # gausspercent
       as.integer(coarse), as.integer(length(coarse)), # coarse
       as.integer(ref.edge), as.integer(length(ref.edge)),# def of reference area
       as.double(percent), as.integer(length(percent)), # to be forest
       gauss,
       binary,
       decreased,
       nrdecr,
       ncdecr,
       refarea, areathreshold,
       result,
       DUP=FALSE, PACKAGE="RandomFields")

    dim(gauss) <- dim(binary) <- dim(decreased) <- dim(refarea) <- c(nx, nx)

    decreased <- decreased[1: (nrdecr * ncdecr)]
    refarea <- refarea[1: (nrdecr * ncdecr)]
    dim(decreased) <-  dim(refarea) <- c(nrdecr, ncdecr)

    if (nrdecr < 20)
      Print(gauss, binary,
            decreased, refarea,
            sum(refarea >= areathreshold),
            areathreshold,
            nrdecr, ncdecr,
            result
            )
 
    par(mfcol=c(2, 2), mar=c(2,2,0,0))
    bw <- c("white", "green")
    image(gauss, col=rainbow(100))
    image(binary, col=bw)
    image(decreased, col=bw)
    image(refarea >= areathreshold, col=bw)

   if (nrdecr < 50) { 
     cat("\n", apply(binary, 1, function(v) {paste(paste(v, collapse=""), "\n")}))
     cat("\n", apply(decreased, 1, function(v) {
       paste(paste(format(v, wi=2), collapse=""), "\n")}))
     cat("\n", apply(refarea, 1, function(v) {
      paste(paste(format(v, wi=2), collapse=""), "\n")}))
   }
  } else {
    .C("analyseForst",
       as.integer(0), # keynr
       as.double(gausspercent), as.integer(length(gausspercent)), # gausspercent
       as.integer(coarse), as.integer(length(coarse)), # coarse
       as.integer(ref.edge), as.integer(length(ref.edge)),# def of reference area
       as.double(percent), as.integer(length(percent)), # to be forest
       result,
       DUP=FALSE, PACKAGE="RandomFields")
    gauss <- binary <- decreased <- refarea <- NULL
  }
  if (length(result) > 1)  {
    dim(result) <- len
    dimnames(result) <- names
  }
  return(list(result=result, gauss=gauss, binary=binary, decreased=decreased, refarea=refarea))
}
