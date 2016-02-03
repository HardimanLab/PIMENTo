heatmapPreprocess <- function(inputData,inputFile,pathwaysDir,method,filenames){
  
  # Read in contents of current file
  geneList <- readLines(paste0(pathwaysDir,"/",inputFile))
  
  if (geneList[1] %in% c("gene","sample","both")) {
    cluster <- geneList[1]
    geneList <- geneList[-1]
  } else {
    warning(paste0("No tree status indicated for ",inputFile,", will skip this file until format corrected"))
    return(FALSE)
  }
  # Graph title is now first line
  if ((!geneList[1] %in% inputData$data[,1]) && !(geneList[1] %in% inputData$data[,2])) {
    graphTitle <- geneList[1]
    geneList <- geneList[-1]
  } else{
    warning(paste0("No graph title provided for ",inputFile,", will skip this file until format corrected"))
    return(FALSE)
  }
  if (!dir.exists(paste0(filenames$outputDir,"/",pathwaysDir))){
    print(paste0("Creating output directory at ",filenames$outputDir,"/",pathwaysDir))
    dir.create(paste0(filenames$outputDir,"/",pathwaysDir), showWarnings = FALSE)
  }

  colnames(inputData$data)[1] <- "Symbol"
  colnames(inputData$data)[2] <- "GeneID"

  # Extract gene name and ID
  if(method == "symbol") {
    heatmapID <- toupper(inputData$data$Symbol)
    geneList <- toupper(geneList)
  }
  else if (method == "geneid") {
    heatmapID <- toupper(inputData$data$GeneID)
  }
  
  # Filter down to desired genes
  matchingGenes <- na.omit(match(geneList,heatmapID))
  heatmapGenes <- inputData$data[matchingGenes,]
  heatmapValues <- inputData$data[matchingGenes,inputData$dataCol]
  heatmapValues[heatmapValues==0] <- 1
  
  median <- apply(heatmapValues,1,median)

  removeIndices <- which(median < 15)
  
  # Extract desired values, check if desired best p-values
  if (length(removeIndices) > 0) {
    heatmapGenes <- heatmapGenes[-removeIndices,]
    heatmapValues <- heatmapValues[-removeIndices,]
    if (dim(heatmapValues)[1] == 0) {
      warning(paste0("No genes from file ", inputFile, " have sufficient p-value or base mean. Moving to next file in subset."))
      return(FALSE)
    }
  }
  # Use available columns from the list desired in #71 for 
  # output csv of truncated data to just desired genes
  outputCSV <- cbind(heatmapGenes[,1:2],heatmapValues)
  write.csv(outputCSV,file=filenames$csv,row.names=F)
  return(list(genes=heatmapGenes,values=heatmapValues,title=graphTitle,cluster=cluster))
}