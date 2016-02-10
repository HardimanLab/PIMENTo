#' @title Generate MA plots of raw and normalized data
#' @description Normalize input raw data using quantile and mloess methods. 
#' Plots of the normalized data along with a dendrogram clustering all samples 
#' will be stored in newly created pipeline directory.
#' @usage preprocessPlots(runSAM.obj, pathwaysDir, fileFormat=c("geneid",
#' "symbol"))
#' @param runSAM.obj Object returned from call to runSAM
#' @param pathwaysDir Directory containing files of genes output from pathway 
#' analysis
#' @param fileFormat Indicator of how genes are identified in each file, be it
#' "geneid" or "symbol"
#' @export

pathwayHeatmap <- function(runSAM.obj,pathwaysDir,
                              fileFormat=c("geneid","symbol")){
    
  outputDir <- paste0(getwd(),"/heatmap_output")
  
  if (!grepl(runSAM.obj$pipelineName,gsub("^.*\\/","",getwd()))) {
      stop(paste0("Function must be run from pipeline output directory"))
  }
  if(is.null(pathwaysDir) || !dir.exists(pathwaysDir) || 
       length(list.files(pathwaysDir)) == 0){
    stop("Provided input directory does not exist or is empty.")
  }
  if(is.null(fileFormat) || !(fileFormat %in% c("geneid","symbol"))){
    stop("No file format provided, use \'geneid\' or \'symbol\'")
  }
  if (!dir.exists(outputDir)) {
    cat("Creating output directory at ",outputDir,"\n",sep="")
    dir.create(outputDir, showWarnings = FALSE)
  }

  
  if(str_sub(pathwaysDir,-1)=="/"){
    pathwaysDir <- substr(pathwaysDir,1,nchar(pathwaysDir)-1)
  }
  

  for (inputFile in list.files(pathwaysDir)){
    
    baseFilename <- basename(file_path_sans_ext(inputFile))
    epsFull <- paste0(outputDir,"/",baseFilename,"-all.eps")
    pdfFull <- paste0(outputDir,"/",baseFilename,"-all.pdf")
    tiffFull <- paste0(outputDir,"/",baseFilename,"-all.tiff")
    xfigFull <- paste0(outputDir,"/",baseFilename,"-all.fig")
    csvFull <- paste0(outputDir,"/",baseFilename,"-filtered.csv")
    fcFull <- paste0(outputDir,"/",baseFilename,"-foldChange.csv")
    fileList <- list(eps=epsFull,pdf=pdfFull,fc=fcFull,tiff=tiffFull,
                     xfig=xfigFull,outputDir=outputDir,csv=csvFull,
                     pathways=pathwaysDir)
    cat("### Processing",inputFile,"###\n")
    preprocessedData <- heatmapPreprocess(runSAM.obj,inputFile,pathwaysDir,
                                          fileFormat,fileList)
    if (typeof(preprocessedData)=="logical") { break }
    heatmapReady <- heatmapFC(preprocessedData,fileList,inputFile)
    if (typeof(heatmapReady)=="logical") { break }
    heatmapMake(heatmapReady,preprocessedData$title,preprocessedData$cluster,
                fileList) 
  }
}