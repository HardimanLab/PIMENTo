#' @title Create heatmaps on desired subsets of genes
#' @description Normalize input raw data using quantile and mloess methods. 
#' Plots of the normalized data along with a dendrogram clustering all samples 
#' will be stored in newly created pipeline directory.
#' @usage CreateHeatmap(sig.genes.sam.obj, subsets.dir, file.format=c("geneid",
#' "symbol"))
#' @param sig.genes.sam.obj Object returned from call to SigGenesSAM
#' @param subsets.dir Directory containing files of genes output from pathway 
#' analysis or simply genes of interest
#' @param file.format Indicator of how genes are identified in each file, be it
#' "geneid" or "symbol"
#' @export

CreateHeatmap <- function(sig.genes.sam.obj, subsets.dir, file.format){
    
  output.dir <- paste0(getwd(),"/heatmap_output")
  
  if (!grepl(sig.genes.sam.obj$pipeline.name,gsub("^.*\\/","",getwd()))) {
      stop(paste0("Function must be run from pipeline output directory"))
  }
  if(is.null(subsets.dir) || !dir.exists(subsets.dir) || 
       length(list.files(subsets.dir)) == 0){
    stop("Provided input directory does not exist or is empty.")
  }
  if(is.null(file.format) || !(file.format %in% c("geneid","symbol"))){
    stop("No file format provided, use \'geneid\' or \'symbol\'")
  }
  if (!dir.exists(output.dir)) {
    cat("Creating output directory at ",output.dir,"\n",sep="")
    dir.create(output.dir, showWarnings = FALSE)
  } 
  if(stringr::str_sub(subsets.dir,-1)=="/"){
    subsets.dir <- substr(subsets.dir,1,nchar(subsets.dir)-1)
  }
  

  for (input.file in list.files(subsets.dir)){
    
    base.filename <- basename(tools::file_path_sans_ext(input.file))
    eps.full <- paste0(output.dir,"/",subsets.dir,"/",base.filename,"_all.eps")
    pdf.full <- paste0(output.dir,"/",subsets.dir,"/",base.filename,"_all.pdf")
    tiff.full <- paste0(output.dir,"/",subsets.dir,"/",base.filename,"_all.tiff")
    xfig.full <- paste0(output.dir,"/",subsets.dir,"/",base.filename,"_all.fig")
    csv.full <- paste0(output.dir,"/",subsets.dir,"/",base.filename,"_filtered.csv")
    fc.full <- paste0(output.dir,"/",subsets.dir,"/",base.filename,"_fold_change.csv")
    file.list <- list(eps=eps.full,pdf=pdf.full,fc=fc.full,tiff=tiff.full,
                     xfig=xfig.full,output.dir=output.dir,csv=csv.full,
                     subsets=subsets.dir)
    cat("### Processing",input.file,"###\n")
    
    preprocessed.data <- HeatmapPreprocess(sig.genes.sam.obj,input.file,subsets.dir,
                                          file.format,file.list)
    
    if (typeof(preprocessed.data)=="logical") {
      next
    }
    
    heatmap.ready <- HeatmapFC(preprocessed.data,file.list,input.file)
    
    if (typeof(heatmap.ready)=="logical") {
      next
    }
    
    HeatmapMake(heatmap.ready, preprocessed.data$title, 
                preprocessed.data$cluster, file.list, subsets.dir) 
  }
}
