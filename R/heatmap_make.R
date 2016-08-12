#' @import gplots

HeatmapMake <- function(heatmap.data, graph.title, cluster, filenames, 
                        subsets.dir) {  
  
  heatmap.colors <- c("#000033", "#333366", "#666699", "#9999CC", "#CCCCFF",
                      "#EEEEFF", "#FFFFFF", "#FFEEEE", "#FFCCCC", "#FF9999",
                      "#FF6666", "#FF3333", "#CC0000")
  
  if (grepl("gene", cluster)) {
    dendro.status <- "row"
    col.v.status <- FALSE
    row.v.status <- TRUE
  } else if (grepl("sample", cluster)) {
    dendro.status <- "column"
    col.v.status <- TRUE
    row.v.status <- FALSE
  } else if (grepl("both", cluster)) {
    dendro.status <- "both"
    col.v.status <- TRUE
    row.v.status <- TRUE
  }

  # Plot the heatmaps to EPS, PDF, TIFF
  setEPS()
  postscript(filenames$eps)
  garbage.collect <- capture.output(
    heatmap.2(heatmap.data,
              hclustfun=function(x) hclust(x, method="ward.D2"),
              distfun=function(x) dist(x, method="euclidean"),
              main=graph.title,
              trace="none",
              margins=c(12, 9),
              dendrogram=dendro.status,
              density.info="none",
              Rowv=row.v.status,
              Colv=col.v.status,
              breaks=c(-7:-1/7,  1:7/7),
              cexRow=1, cexCol=1, srtCol=45,
              col=heatmap.colors))
  garbage <- dev.off()               # close the EPS device
  

  lhei <- c(0.5, 0.75, 4)
  lwid <- c(1.5, 4)
  lmat <- rbind(c(4, 0), c(0, 3), c(2, 1))

  pdf(filenames$pdf, height=17, width=8.5)
  garbage.collect <- capture.output(
    heatmap.2(heatmap.data,
              hclustfun=function(x) hclust(x, method="ward.D2"),
              distfun=function(x) dist(x, method="euclidean"),
              main = graph.title,
              trace="none",
              margins=c(12, 9),
              dendrogram=dendro.status,
              density.info="none",
              Rowv=row.v.status,
              Colv=col.v.status,
              keysize=1,
              breaks=c(-7:-1/7,  1:7/7),
              cexRow=1, cexCol=1, srtCol=45,
	            lmat=lmat, lwid=lwid, lhei=lhei,
              col=heatmap.colors))
  garbage <- dev.off()               # close the PDF device
  
  suppressWarnings(xfig(filenames$xfig))
  garbage.collect <- capture.output(
    heatmap.2(heatmap.data,
              hclustfun=function(x) hclust(x,method="ward.D2"),
              distfun=function(x) dist(x,method="euclidean"),
              main=graph.title,
              trace="none",         
              margins=c(12, 9),
              dendrogram=dendro.status,
              density.info="none",
              Rowv=row.v.status,
              Colv=col.v.status,
              breaks=c(-7:-1/7, 1:7/7),
              cexRow=1,cexCol=1,srtCol=45,
              col=heatmap.colors))
  garbage <- dev.off()               # close the FIG device
  
  if(capabilities("tiff")) {
    tiff(filenames$tiff)
    garbage.collect <- capture.output(
      heatmap.2(heatmap.data,
                hclustfun=function(x) hclust(x, method="ward.D2"),
                distfun=function(x) dist(x, method="euclidean"),
                main=graph.title, 
                trace="none",         
                margins=c(12, 9),
                dendrogram=dendro.status,
                density.info="none",
                Rowv=row.v.status,
                Colv=col.v.status,
                breaks=c(-7:-1/7, 1:7/7),
                cexRow=1, cexCol=1, srtCol=45,
                col=heatmap.colors))
    garbage <- dev.off()  # close the TIFF device
  }

  cat(paste0("All plots have been created at ./heatmap_output/",subsets.dir,
             "\n"),sep="")
}
