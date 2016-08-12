#' @import ggplot2

SampleSimilarity <- function(sig.genes.sam.obj) {

  groups <- sapply(sig.genes.sam.obj$response, function(sample) 
    ifelse(sample == 1, "control", "experiment"))
  
  if ("class.compare.name" %in% names(sig.genes.sam.obj)) {
    ps.similarity <- paste0("./", sig.genes.sam.obj$pipeline.name, "_pipeline/",
                            sig.genes.sam.obj$pipeline.name, 
                            "_sample_similarity_", 
                            sig.genes.sam.obj$class.compare.name, ".ps")  
    pdf.similarity <- paste0("./", sig.genes.sam.obj$pipeline.name, 
                             "_pipeline/",
                             sig.genes.sam.obj$pipeline.name, 
                             "_sample_similarity_",
                             sig.genes.sam.obj$class.compare.name, ".pdf")  
  }
  else {
    ps.similarity <- paste0("./", sig.genes.sam.obj$pipeline.name, "_pipeline/",
                            sig.genes.sam.obj$pipeline.name, 
                            "_sample_similarity.ps")  
    pdf.similarity <- paste0("./", sig.genes.sam.obj$pipeline.name, 
                             "_pipeline/", sig.genes.sam.obj$pipeline.name, 
                             "_sample_similarity.pdf")  
  }

  postscript(file=ps.similarity, paper="letter")
  
  ## Create heatmap plot of sample similarity
  if ("class.compare.cols" %in% names(sig.genes.sam.obj)) {
    norm.data <- 
      as.matrix(sig.genes.sam.obj$data[, sig.genes.sam.obj$class.compare.cols])
  }
  else {
    norm.data <- as.matrix(sig.genes.sam.obj$data[, sig.genes.sam.obj$data.col])
  }
  storage.mode(norm.data) <- "integer"
  rld <- DESeq2::rlog(norm.data, fitType="local")
  dist.matrix <- as.matrix(dist(t(rld)))
  colnames(dist.matrix) <- colnames(rld)
  rownames(dist.matrix) <- colnames(dist.matrix)
  
  hclustfunc <- function(x) hclust(x, method="complete")
  distfunc <- function(x) dist(x, method="euclidean")
  cl.row <- hclustfunc(distfunc(dist.matrix))
  cl.col <- hclustfunc(distfunc(t(dist.matrix)))
  gr.row <- cutree(cl.row, 2)
  gr.col <- cutree(cl.col, 2)
  side.colors <- c("red", "blue")
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  garbageCollect <- capture.output(
    gplots::heatmap.2(dist.matrix, symm=TRUE, trace="none", breaks=256, 
                      col=colors, margins=c(12, 9),cexRow=1, cexCol=1, 
                      srtCol=45, ColSideColors=side.colors[gr.col],
                      RowSideColors=side.colors[gr.row])
    )

  ## Create PCA plot of sample similarity
  
  pca <- prcomp(t(rld))
  percent.var <- pca$sdev^2 / sum(pca$sdev^2)
  d <- data.frame(PC1=pca$x[, 1], PC2=pca$x[, 2], group=groups)
  
  p <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
    geom_point() + geom_text(show.legend=F, check_overlap=T, hjust=0, vjust=0,
                             aes(label=rownames(dist.matrix))) + 
    xlab(paste0("PC1: ", round(percent.var[1] * 100), "% variance")) + 
    ylab(paste0("PC2: ", round(percent.var[2] * 100), "% variance")) + 
    coord_fixed()
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  grid::grid.newpage()
  grid::grid.draw(gt)
  
  dev.off()

  pdf(file=pdf.similarity, paper="letter")
  garbageCollect <- capture.output(
    gplots::heatmap.2(dist.matrix, symm=T, trace="none", breaks=256, col=colors,
                      margins=c(12, 9), cexRow=1, cexCol=1, srtCol=45,
                      ColSideColors=side.colors[gr.col],
                      RowSideColors=side.colors[gr.row])
    )
  grid::grid.newpage()
  grid::grid.draw(gt)

  dev.off()
  
  cat("Plots of sample similarity created at ", ps.similarity,sep="")
 
}