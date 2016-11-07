HeatmapPreprocess <- function(input.data, input.file, subsets.dir, 
                              method, filenames){
  
  # Read in contents of current file
  gene.list <- readLines(paste0(subsets.dir, "/", input.file))
  if(length(gene.list) > 502) {
    warning("It is advised to use gene lists less than 500 for visible text.")
  }
  
  if (grepl("gene|sample|both",tolower(gene.list[1]))) {
    cluster <- tolower(gene.list[1])
    gene.list <- gene.list[-1]
  } else {
    warning(paste0("No tree status indicated for ", input.file,
                  ", will skip this", " file until format corrected"))
    return(FALSE)
  }
  # Graph title is now first line
  if (!gene.list[1] %in% input.data$normalized[, 1:input.data$ntext]) {
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
  
  colnames(input.data$normalized)[input.data$symbol.index] <- "Symbol"
  colnames(input.data$normalized)[input.data$id.index] <- "GeneID"

  # Extract gene name and ID
  if(method == "symbol") {
    heatmap.ID <- toupper(input.data$normalized$Symbol)
    gene.list <- toupper(gene.list)
    method.index <- input.data$symbol.index
  } else if (method == "geneid") {
    heatmap.ID <- toupper(input.data$normalized$GeneID)
    gene.list <- toupper(gene.list)
    method.index <- input.data$id.index
  }

  # Filter down to desired genes
  write.table(heatmap.ID, paste0("norm_genes_", input.file, ".txt"))
  matching.genes <- na.omit(match(gene.list,heatmap.ID))
  print(gene.list[attr(na.omit(matching.genes), "na.action")])
  heatmap.genes <- input.data$normalized[matching.genes, ]
  
  if ("class.compare.cols" %in% names(input.data)) {
    heatmap.values <- 
      input.data$normalized[matching.genes, input.data$class.compare.cols]
  } else {
    heatmap.values <- input.data$normalized[matching.genes, input.data$data.col]
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
