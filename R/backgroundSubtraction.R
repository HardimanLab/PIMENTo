#' @title Perform background subtraction on normalized data
#' @description Remove all genes that have an maximum intensity level below that
#' of the provided cutoff.
#' @usage backgroundSubtraction(preprocessData.obj, 
#' method=c("mloess","quantile"), cutoff)
#' @param preprocessData.obj Object returned from call to 
#' preprocessData
#' @param method Type of normalization to use: "quantile" for quantile 
#' normalization; "mloess" for MLOESS normalization
#' @param cutoff The cutoff value (binary logarithm) identified through
#' backgroundCutoff
#' @return A list with components
#' \item{ntext}{Number of leading text columns}
#' \item{dataCol}{Vector of column indices containing array data}
#' \item{id}{Vector containing gene ID's}
#' \item{idInd}{Column index containing gene ID information}
#' \item{symbol}{Vector containing gene symbols}
#' \item{symbolIndex}{Column index containing gene symbol information}
#' \item{descStats}{Vector of column indices containing descriptive statistics}
#' \item{pipelineName}{Name of pipeline generated from input file name sans 
#' extension}
#' \item{data}{Data frame of descriptive stats and normalized data for chosen 
#' method}
#' @export

backgroundSubtraction <- function(preprocessData.obj, method, cutoff) {

  if(missing(cutoff))
    stop("Must submit cutoff value for background subtraction.")
  
  if (grepl("quantile",method))
    backgroundData <- preprocessData.obj$quantile[,preprocessData.obj$dataCol]
  else if (grepl("mloess",method))
    backgroundData <- preprocessData.obj$mloess[,preprocessData.obj$dataCol]
  else
    stop("Input argument 'method' must be either 'quantile' or 'mloess'")
  
  geneMax <- apply(backgroundData,1,max)
  keepIndices <- !(geneMax < 2^cutoff)
  descStats <- preprocessData.obj$descStats[keepIndices,]
  subtractedData <- backgroundData[keepIndices,]
  
  preprocessData.obj$symbol <- preprocessData.obj$symbol[keepIndices]
  preprocessData.obj$id <- preprocessData.obj$id[keepIndices]
  preprocessData.obj$data <- cbind(descStats,subtractedData)
  preprocessData.obj$descStats <- descStats
  preprocessData.obj$mloess <- NULL
  preprocessData.obj$quantile <- NULL
  

  return(preprocessData.obj)
}