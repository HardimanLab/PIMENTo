#' @title Visualize arrays to identify cutoff for background subtraction
#' @description Plot the histogram of max gene illumination across all arrays
#' to identify cutoff for background subtraction. Upper and lower limits can be
#' changed to narrow focus. Plots will be saved in analysis pipeline 
#' directory.
#' @usage backgroundCutoff(initializePipeline.obj, method=c("loess","quantile"), 
#' xlim.lo=0, xlim.hi)
#' @param initializePipeline.obj Object returned from call to 
#' initializePipeline
#' @param method Type of normalization to use: "quantile" or "q" for quantile 
#' normalization; "mloess" or "m" for MLOESS normalization
#' @param xlim.lo Lower bound on X for histogram plot (binary logarithm)
#' @param xlim.hi Upper bound on X for histogram plot (binary logarithm)
#' @export

backgroundCutoff <- function(initializePipeline.obj,method=c("mloess","quantile"),
                             xlim.lo=0,xlim.hi=0) {

  method <- match.arg(method)
  backgroundData <- switch(method,
        "mloess" = initializePipeline.obj$mloess[,initializePipeline.obj$dataCol],
        "quantile" = initializePipeline.obj$quantile[,initializePipeline.obj$dataCol]
  )
  geneMax <- apply(backgroundData,1,max)
  logData <- log2(geneMax)
  
  if(xlim.hi != 0 | xlim.lo != 0) {
      plotsFile <- paste0("./",initializePipeline.obj$pipelineName,"_pipeline/",
                          initializePipeline.obj$pipelineName,"_",method,
                          "_backgroundCutoff_",xlim.lo,"-",xlim.hi,".ps")
      postscript(file=plotsFile,paper="letter")
      hist(logData,breaks=200,xlim=c(xlim.lo,xlim.hi), 
           main = "Maximum Illumination Across Arrays", 
           xlab=expression('Illumination ('*log[2]*')'))
      invisible(dev.off())
      cat("Histogram plots saved at",plotsFile,"\n")
  }
  else {
    plotsFile <- paste0("./",initializePipeline.obj$pipelineName,"_pipeline/",
                        initializePipeline.obj$pipelineName,"_",method,
                        "_backgroundCutoff_full.ps")
    postscript(file=plotsFile,paper="letter")
    hist(logData,breaks=200,main = "Maximum Illumination Across Arrays",
         xlab=expression('Illumination ('*log[2]*')'))
    invisible(dev.off())
    cat("Histogram plots saved at",plotsFile,"\n")
  }
}
