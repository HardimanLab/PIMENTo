heatmapFC <- function(data,filenames,inputFile){
  
  heatmapValues <- data$values
  print(head(heatmapValues))
  baseMedians <- apply(heatmapValues,1,median)
  sigGenesID <- make.names(data$genes[,data$methodIndex],unique=TRUE)
  print(head(sigGenesID))
  
  log2spaceDiff <- log2(heatmapValues) - log2(baseMedians)
  
  # Apply fold-change to boost +/-
  foldChange <- ifelse(log2spaceDiff > 0, 2^log2spaceDiff, (-1)*2^(-log2spaceDiff))
  
  # Write out to CSV the fold-change data
  processedData <- data.frame(sigGenesID,foldChange,row.names=1)
  write.csv(processedData,file=filenames$fc)
  outputDirStr <- regexpr("([^/]+/+[^/]+/+[^/]+)/*$",filenames$fc)
  outputDirStr <- substr(filenames$fc,outputDirStr,
                         attr(outputDirStr,"match.length")+outputDirStr)
  cat(paste0("CSV of fold-change values has been created at ./",outputDirStr,
             "\n"),sep="")
  
  if (dim(foldChange)[1] == 1){
    print(paste0("Only one gene found to be significant - No heatmap creation"))
    return(FALSE)
  }
  
  # Normalize fold-change for heatmap drawing
  posIndices <- which(foldChange >= 1)
  negIndices <- which(foldChange < -1)
  
  # LogBASE2 of the fold-changes
  foldChange[posIndices] <- log2(foldChange[posIndices])
  foldChange[negIndices] <- -log2(-foldChange[negIndices])
  
  # Global normalization, divide by max value
  foldChange.global <- foldChange/max(abs(foldChange))
  
  # Square-root scaling
  foldChange.global.norm <- ifelse(foldChange.global > 0, foldChange.global^0.5,
                                   -(-(foldChange.global))^0.5)
  
  # Create data frame combining labels with normalized fold-change values
  processedData.global.norm <- data.frame(sigGenesID,foldChange.global.norm,row.names=1)
  matProcessedData.global.norm <- data.matrix(processedData.global.norm)
  return(matProcessedData.global.norm)
}
