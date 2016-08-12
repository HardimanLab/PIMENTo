#' @title Visualize arrays to identify cutoff for background subtraction
#' @description Plot the histogram of max gene illumination across all arrays
#' to identify cutoff for background subtraction. Upper and lower limits can be
#' changed to narrow focus. Plots will be saved in analysis pipeline 
#' directory.
#' @usage BackgroundCutoff(preprocess.data.obj, method=c("mloess","quantile"), 
#' xlim.lo=0, xlim.hi)
#' @param preprocess.data.obj Object returned from call to 
#' preprocess.data
#' @param method Type of normalization to use: "quantile" for quantile 
#' normalization; "mloess" for MLOESS normalization
#' @param xlim.lo Lower bound on X for histogram plot (binary logarithm)
#' @param xlim.hi Upper bound on X for histogram plot (binary logarithm)
#' @export

BackgroundCutoff <- function(preprocess.data.obj, method,
                             xlim.lo=0, xlim.hi=0) {

  if (grepl("quantile",method))
    background.data <- 
      preprocess.data.obj$quantile[, preprocess.data.obj$data.col]
  else if (grepl("mloess",method))
    background.data <- 
      preprocess.data.obj$mloess[, preprocess.data.obj$data.col]
  else
    stop("Input argument 'method' must be either 'quantile' or 'mloess'")
  
  if (xlim.hi < xlim.lo)
    stop("Input argument 'xlim.lo' must be less than 'xlim.hi'")

  gene.max <- apply(background.data, 1, max)
  gene.max[gene.max==0] <- 1
  log.data <- log2(gene.max)
  
  if (xlim.hi != 0 | xlim.lo != 0) {
      ps.plot <- paste0("./", preprocess.data.obj$pipeline.name, "_pipeline/",
                          preprocess.data.obj$pipeline.name, "_", method,
                          "_background_cutoff_", xlim.lo, "_", xlim.hi, ".ps")
      pdf.plot <- paste0("./", preprocess.data.obj$pipeline.name, "_pipeline/",
                          preprocess.data.obj$pipeline.name, "_", method,
                          "_background_cutoff_", xlim.lo, "_", xlim.hi, ".pdf")

      postscript(file=ps.plot, paper="letter")
      hist(log.data, breaks=200, xlim=c(xlim.lo, xlim.hi), 
           main = "Maximum Illumination Across Arrays", 
           xlab=expression('Illumination ('*log[2]*')'))
      garbage <- dev.off()

      pdf(file=pdf.plot, paper="letter")
      hist(log.data, breaks=200, xlim=c(xlim.lo, xlim.hi), 
           main = "Maximum Illumination Across Arrays", 
           xlab=expression('Illumination ('*log[2]*')'))
      garbage <- dev.off()

      cat("Histogram plots saved at", ps.plot, "\n")
  }
  else {
    ps.plot <- paste0("./", preprocess.data.obj$pipeline.name, "_pipeline/",
                        preprocess.data.obj$pipeline.name, "_", method,
                        "_background_cutoff_full.ps")
    pdf.plot <- paste0("./", preprocess.data.obj$pipeline.name, "_pipeline/",
                        preprocess.data.obj$pipeline.name, "_", method,
                        "_background_cutoff_full.pdf")

    postscript(file=ps.plot, paper="letter")
    hist(log.data, breaks=200, main = "Maximum Illumination Across Arrays",
         xlab=expression('Illumination ('*log[2]*')'))
    garbage <- dev.off()

    pdf(file=pdf.plot, paper="letter")
    hist(log.data, breaks=200, main = "Maximum Illumination Across Arrays",
         xlab=expression('Illumination ('*log[2]*')'))
    garbage <- dev.off()

    cat("Histogram plots saved at", ps.plot,"\n")
  }
}
