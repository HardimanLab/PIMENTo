#' @import ggplot2

sampleSimilarity <- function(runSAM.obj) {

  groups = sapply(runSAM.obj$response,function(sample) 
    ifelse(sample == 1,"control","experiment"))
  
  if ("classCompareName" %in% names(runSAM.obj)) {
    ps.similarityFile=paste0("./",runSAM.obj$pipelineName,"_pipeline/",
                          runSAM.obj$pipelineName,"_sampleSimilarity-",
                          runSAM.obj$classCompareName,".ps")  
    pdf.similarityFile=paste0("./",runSAM.obj$pipelineName,"_pipeline/",
                          runSAM.obj$pipelineName,"_sampleSimilarity-",
                          runSAM.obj$classCompareName,".pdf")  
  }
  else {
    ps.similarityFile=paste0("./",runSAM.obj$pipelineName,"_pipeline/",
                          runSAM.obj$pipelineName,"_sampleSimilarity.ps")  
    pdf.similarityFile=paste0("./",runSAM.obj$pipelineName,"_pipeline/",
                          runSAM.obj$pipelineName,"_sampleSimilarity.pdf")  
  }

  postscript(file=ps.similarityFile,paper="letter")
  
  ## Create heatmap plot of sample similarity
  if ("classCompareCols" %in% names(runSAM.obj)) {
    normData <- as.matrix(runSAM.obj$data[,runSAM.obj$classCompareCols])
  }
  else {
    normData <- as.matrix(runSAM.obj$data[,runSAM.obj$dataCol])
  }
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
  garbageCollect <- capture.output(
    gplots::heatmap.2(distMatrix,symm=TRUE,trace="none",breaks=256,col=colours,
                      margins=c(12,9),cexRow=1,cexCol=1,srtCol=45,
                      ColSideColors=sideColours[gr.col],
                      RowSideColors=sideColours[gr.row])
    )

  ## Create PCA plot of sample similarity
  
  pca <- prcomp(t(rld))
  percentVar = pca$sdev^2/sum(pca$sdev^2)
  d <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], group=groups)
  
  p <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
          geom_point() + geom_text(show.legend=FALSE, check_overlap=TRUE, 
                                   hjust=0, vjust=0, aes(label = rownames(distMatrix))) + 
          xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
          ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
          coord_fixed()
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid::grid.newpage()
  grid::grid.draw(gt)
  
  dev.off()

  pdf(file=pdf.similarityFile,paper="letter")
  garbageCollect <- capture.output(
    gplots::heatmap.2(distMatrix,symm=TRUE,trace="none",breaks=256,col=colours,
                      margins=c(12,9),cexRow=1,cexCol=1,srtCol=45,
                      ColSideColors=sideColours[gr.col],
                      RowSideColors=sideColours[gr.row])
    )
  grid::grid.newpage()
  grid::grid.draw(gt)

  dev.off()
  
  cat("Plots of sample similarity created at ",ps.similarityFile,sep="")
 
}
