HeatmapFC <- function(data,filenames,input.file){
  
  heatmap.values <- data$values
  base.medians <- apply(heatmap.values, 1, median)
  sig.genes.ID <- make.names(data$genes[, data$method.index], unique=TRUE)
  log2space.diff <- log2(heatmap.values) - log2(base.medians)
  
  # Apply fold-change to boost +/-
  fold.change <- ifelse(log2space.diff > 0, 
                        2^log2space.diff, 
                        (-1)*2^(-log2space.diff))
  
  # Write out to CSV the fold-change data
  process.data <- data.frame(sig.genes.ID, fold.change, row.names=1)
  write.csv(process.data, file=filenames$fc)
  output.dir.str <- regexpr("([^/]+/+[^/]+/+[^/]+)/*$", filenames$fc)
  output.dir.str <- substr(filenames$fc, output.dir.str,
                         attr(output.dir.str, "match.length") + output.dir.str)
  cat(paste0("CSV of fold-change values has been created at ./", output.dir.str,
             "\n"), sep="")
  
  if (dim(fold.change)[1] == 1){
    print(paste0("Only one gene found to be significant - No heatmap creation"))
    return(FALSE)
  }
  
  # Normalize fold-change for heatmap drawing
  pos.indices <- which(fold.change >= 1)
  neg.indices <- which(fold.change < -1)
  
  # LogBASE2 of the fold-changes
  fold.change[pos.indices] <- log2(fold.change[pos.indices])
  fold.change[neg.indices] <- -log2(-fold.change[neg.indices])
  
  # Global normalization, divide by max value
  fold.change.global <- fold.change/max(abs(fold.change))
  
  # Square-root scaling
  fold.change.global.norm <- ifelse(fold.change.global > 0, 
                                    fold.change.global^0.5,
                                    -(-(fold.change.global))^0.5)
  
  # Create data frame combining labels with normalized fold-change values
  process.data.global.norm <- data.frame(sig.genes.ID, fold.change.global.norm, 
                                         row.names=1)
  matprocess.data.global.norm <- data.matrix(process.data.global.norm)
  return(matprocess.data.global.norm)
}