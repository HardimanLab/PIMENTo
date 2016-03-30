#' @title Visualize arrays to identify cutoff for background subtraction
#' @description Plot the histogram of max gene illumination across all arrays
#' to identify cutoff for background subtraction. Upper and lower limits can be
#' changed to narrow focus. Plots will be saved in analysis pipeline 
#' directory.
#' @usage backgroundCutoff(preprocessData.obj, method=c("mloess","quantile"), 
#' xlim.lo=0, xlim.hi)
#' @param preprocessData.obj Object returned from call to 
#' preprocessData
#' @param method Type of normalization to use: "quantile" for quantile 
#' normalization; "mloess" for MLOESS normalization
#' @param xlim.lo Lower bound on X for histogram plot (binary logarithm)
#' @param xlim.hi Upper bound on X for histogram plot (binary logarithm)
#' @export

backgroundCutoff <- function(preprocessData.obj, method,
                             xlim.lo=0,xlim.hi=0) {

  if (grepl("quantile",method))
    backgroundData <- preprocessData.obj$quantile[,preprocessData.obj$dataCol]
  else if (grepl("mloess",method))
    backgroundData <- preprocessData.obj$mloess[,preprocessData.obj$dataCol]
  else
    stop("Input argument 'method' must be either 'quantile' or 'mloess'")
  
  geneMax <- apply(backgroundData,1,max)
  logData <- log2(geneMax)
  
  if(xlim.hi != 0 | xlim.lo != 0) {
      ps.plotsFile <- paste0("./",preprocessData.obj$pipelineName,"_pipeline/",
                          preprocessData.obj$pipelineName,"_",method,
                          "_backgroundCutoff_",xlim.lo,"-",xlim.hi,".ps")
      pdf.plotsFile <- paste0("./",preprocessData.obj$pipelineName,"_pipeline/",
                          preprocessData.obj$pipelineName,"_",method,
                          "_backgroundCutoff_",xlim.lo,"-",xlim.hi,".pdf")

      postscript(file=ps.plotsFile,paper="letter")
      hist(logData,breaks=200,xlim=c(xlim.lo,xlim.hi), 
           main = "Maximum Illumination Across Arrays", 
           xlab=expression('Illumination ('*log[2]*')'))
      garbage <- dev.off()

      pdf(file=pdf.plotsFile,paper="letter")
      hist(logData,breaks=200,xlim=c(xlim.lo,xlim.hi), 
           main = "Maximum Illumination Across Arrays", 
           xlab=expression('Illumination ('*log[2]*')'))
      garbage <- dev.off()

      cat("Histogram plots saved at",ps.plotsFile,"\n")
  }
  else {
    ps.plotsFile <- paste0("./",preprocessData.obj$pipelineName,"_pipeline/",
                        preprocessData.obj$pipelineName,"_",method,
                        "_backgroundCutoff_full.ps")
    pdf.plotsFile <- paste0("./",preprocessData.obj$pipelineName,"_pipeline/",
                        preprocessData.obj$pipelineName,"_",method,
                        "_backgroundCutoff_full.pdf")

    postscript(file=ps.plotsFile,paper="letter")
    hist(logData,breaks=200,main = "Maximum Illumination Across Arrays",
         xlab=expression('Illumination ('*log[2]*')'))
    garbage <- dev.off()

    pdf(file=pdf.plotsFile,paper="letter")
    hist(logData,breaks=200,main = "Maximum Illumination Across Arrays",
         xlab=expression('Illumination ('*log[2]*')'))
    garbage <- dev.off()

    cat("Histogram plots saved at",ps.plotsFile,"\n")
  }
}
