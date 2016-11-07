#' @title Perform background subtraction on normalized data
#' @description Remove all genes that have an maximum intensity level below that
#' of the provided cutoff.
#' @usage BackgroundSubtraction(preprocess.data.obj, 
#' method=c("mloess", "quantile"), cutoff)
#' @param preprocess.data.obj Object returned from call to 
#' preprocess.data
#' @param method Type of normalization to use: "quantile" for quantile 
#' normalization; "mloess" for MLOESS normalization; "raw" for raw data
#' @param cutoff The cutoff value (binary logarithm) identified through
#' background.cutoff
#' @return A list with components
#' \item{ntext}{Number of leading text columns}
#' \item{data.col}{Vector of column indices containing array data}
#' \item{id}{Vector containing gene ID's}
#' \item{symbol}{Vector containing gene symbols}
#' \item{desc.stats}{Vector of column indices containing descriptive statistics}
#' \item{pipeline.name}{Name of pipeline generated from input file name sans 
#' extension}
#' \item{normalized}{Data frame of descriptive stats and normalized data for 
#' the chosen method}
#' @export

BackgroundSubtraction <- function(preprocess.data.obj, method, cutoff) {

  if(missing(cutoff))
    stop("Must submit cutoff value for background subtraction.")
  
  if (grepl("quantile", method)) {
    background.data <- preprocess.data.obj$quantile[, preprocess.data.obj$data.col]
  }
  else if (grepl("mloess", method)) {
    background.data <- preprocess.data.obj$mloess[, preprocess.data.obj$data.col]
  }
  else if (grepl("raw", method)) {
    background.data <- preprocess.data.obj$raw[, preprocess.data.obj$data.col]
  }
  else
    stop("Input argument 'method' must be either 'raw', 'quantile' or 'mloess'")
  
  gene.max <- apply(background.data, 1, max)
  keep.indices <- !(gene.max < 2^cutoff)
  desc.stats <- preprocess.data.obj$desc.stats[keep.indices, ]
  subtracted.data <- background.data[keep.indices, ]
  
  preprocess.data.obj$symbol <- preprocess.data.obj$symbol[keep.indices]
  preprocess.data.obj$id <- preprocess.data.obj$id[keep.indices]
  preprocess.data.obj$normalized <- cbind(desc.stats, subtracted.data)
  preprocess.data.obj$desc.stats <- desc.stats
  preprocess.data.obj$background.cutoff <- cutoff
  preprocess.data.obj$normalize.method <- method
  preprocess.data.obj$mloess <- NULL
  preprocess.data.obj$quantile <- NULL

  pipeline.name <- preprocess.data.obj$pipeline.name

  write.csv(preprocess.data.obj$normalized, 
            paste0("./", pipeline.name, "_pipeline/", pipeline.name, 
                   "_backgroundSub_", method, "_data.csv"), row.names=F)

  return(preprocess.data.obj)
}
