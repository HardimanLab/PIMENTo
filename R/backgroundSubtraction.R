#' @title Perform background subtraction on normalized data
#' @description Remove all genes that have an maximum intensity level below that
#' of the provided cutoff.
#' @usage backgroundSubtraction(initializePipeline.obj, 
#' method=c("mloess","quantile"), cutoff)
#' @param initializePipeline.obj Object returned from call to 
#' initializePipeline
#' @param method Type of normalization to use: "quantile" or "q" for quantile 
#' normalization; "mloess" or "m" for mloess normalization
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

backgroundSubtraction <- function(initializePipeline.obj,
                                  method=c("mloess","quantile"),cutoff) {

  if(missing(cutoff))
    stop("Must submit cutoff value for background subtraction.")
  
  method <- tryCatch(match.arg(method), error=function(e){
    stop("The provided method of", method, "is not one of the avaiable methods.",
      " Use 'quantile' or 'mloess'",call.=FALSE)
  })
  if(method=="quantile")
    backgroundData <- initializePipeline.obj$quantile[,initializePipeline.obj$dataCol]
  else if (method=="mloess")
    backgroundData <- initializePipeline.obj$mloess[,initializePipeline.obj$dataCol]
  
  descStats <- initializePipeline.obj$descStats
  geneMax <- apply(backgroundData,1,max)
  keepIndices <- !(geneMax < 2^cutoff)
  descStats <- descStats[keepIndices,]
  subtractedData <- backgroundData[keepIndices,]
  initializePipeline.obj$data <- cbind(descStats,subtractedData)
  initializePipeline.obj$mloess <- NULL
  initializePipeline.obj$quantile <- NULL
  

  return(initializePipeline.obj)
}