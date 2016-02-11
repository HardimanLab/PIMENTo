#' @title Generate plot of sample identify matrix and PCA
#' @description Heatmap of sample similarity will be created along with a PCA
#' plot labelled with control and experimental groups.
#' @usage sampleSimilarity(runSAM.obj)
#' @param runSAM.obj Object returned from call to runSAM
#' @import ggplot2
#' @export

sampleSimilarity <- function(runSAM.obj) {

  groups = sapply(runSAM.obj$response,function(sample) 
    ifelse(sample == 1,"control","experiment"))
  
  similarityFile=paste0("./",runSAM.obj$pipelineName,"_pipeline/",
                        runSAM.obj$pipelineName,"_sampleSimilarity.ps")  
  postscript(file=similarityFile,paper="letter")
  
  ## Create heatmap plot of sample similarity
  
  normData <- as.matrix(runSAM.obj$data[,runSAM.obj$dataCol])
  storage.mode(normData) <- "integer"
  rld <- DESeq2::rlog(normData,fitType="local")
  distMatrix <- as.matrix(dist(t(rld)))
  rownames(distMatrix) = colnames(distMatrix) = colnames(rld)
  
  hclustfunc <- function(x) hclust(x, method="complete")
  distfunc <- function(x) dist(x, method="euclidean")
  cl.row <- hclustfunc(distfunc(distMatrix))
  cl.col <- hclustfunc(distfunc(t(distMatrix)))
  gr.row <- cutree(cl.row, 2)
  gr.col <- cutree(cl.col, 2)
  sideColours <- c("red","blue")
  colours = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  gplots::heatmap.2(distMatrix,symm=TRUE,trace="none",breaks=256,
            col=colours,margins=c(12,9),ColSideColors=sideColours[gr.col],
            RowSideColors=sideColours[gr.row])

  ## Create PCA plot of sample similarity
  
  pca <- prcomp(t(rld))
  percentVar = pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group=groups)
  
  print(ggplot(data = d,aes_string(x ="PC1",y="PC2",color="group"))+
          geom_point(size=3) + 
          xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
          ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
          coord_fixed())
  
  invisible(dev.off())
  
  cat("Plots of sample similarity created at ",similarityFile,sep="")
 
}