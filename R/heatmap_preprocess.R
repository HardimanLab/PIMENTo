HeatmapPreprocess <- function(input.data, input.file, subsets.dir, 
                              method, filenames){
  
  # Read in contents of current file
  gene.list <- readLines(paste0(subsets.dir, "/", input.file))
  
  if (grepl("gene|sample|both",tolower(gene.list[1]))) {
    cluster <- tolower(gene.list[1])
    gene.list <- gene.list[-1]
  } else {
    warning(paste0("No tree status indicated for ", input.file,
                  ", will skip this", " file until format corrected"))
    return(FALSE)
  }
  # Graph title is now first line
  if (!gene.list[1] %in% input.data$data[, 1:input.data$ntext]) {
    graph.title <- gene.list[1]
    gene.list <- gene.list[-1]
  } else{
    warning(paste0("No graph title provided for ", input.file, 
                   ", will skip this", " file until format corrected"))
    return(FALSE)
  }
  if (!dir.exists(paste0(filenames$output.dir, "/",subsets.dir))){
    print(paste0("Creating output directory at ",
                 filenames$output.dir, "/", subsets.dir))
    dir.create(paste0(filenames$output.dir, "/", subsets.dir), 
               showWarnings = FALSE)
  }
  colnames(input.data$data)[input.data$symbol.index] <- "Symbol"
  colnames(input.data$data)[input.data$id.index] <- "GeneID"

  # Extract gene name and ID
  if(method == "symbol") {
    heatmap.ID <- toupper(input.data$data$Symbol)
    gene.list <- toupper(gene.list)
    method.index <- input.data$symbol.index
  } else if (method == "geneid") {
    heatmap.ID <- toupper(input.data$data$GeneID)
    gene.list <- toupper(gene.list)
    method.index <- input.data$id.index
  }

  # Filter down to desired genes
  matching.genes <- na.omit(match(gene.list,heatmap.ID))
  heatmap.genes <- input.data$data[matching.genes, ]
  
  if ("class.compare.cols" %in% names(input.data)) {
    heatmap.values <- 
      input.data$data[matching.genes, input.data$class.compare.cols]
  } else {
    heatmap.values <- input.data$data[matching.genes, input.data$data.col]
  }
  
  temp <- as.matrix(heatmap.values)
  temp[temp == 0] <- 1
  heatmap.values <- as.data.frame(temp)
#  heatmap.values[heatmap.values==0] <- 1
    
  if (nrow(heatmap.values)[1] == 0) {
    warning(paste0("No genes from file ", input.file, 
                   " match those in dataset.", 
                   " Moving to next file in subset"))
    return(FALSE)
  }
  output.csv <- cbind(heatmap.genes[, 1:2], heatmap.values)
  write.csv(output.csv, file=filenames$csv, row.names=F)
  return(list(genes=heatmap.genes, values=heatmap.values, title=graph.title,
              cluster=cluster, method.index=method.index))
}
