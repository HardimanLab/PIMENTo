heatmapMake <- function(heatmapData,graphTitle,cluster,filenames,pathwaysDir) {
  
  library(gplots)
  
  heatmapColors <- c("#000033","#333366","#666699","#9999CC","#CCCCFF",
                      "#EEEEFF","#FFFFFF","#FFEEEE","#FFCCCC","#FF9999",
                      "#FF6666","#FF3333","#CC0000")
  
  if (grepl("gene",cluster)) {
    dendroStatus="row"
    colV=FALSE
    rowV=TRUE
  } else if (grepl("sample",cluster)) {
    dendroStatus="column"
    colV=TRUE
    rowV=FALSE
  } else if (grepl("both",cluster)) {
    dendroStatus="both"
    colV=TRUE
    rowV=TRUE
  }

  # Plot the heatmaps to EPS, PDF, TIFF
  setEPS()
  postscript(filenames$eps)
  garbageCollect <- capture.output(heatmap.2(heatmapData,
            hclustfun=function(x) hclust(x,method="ward.D2"),
            distfun=function(x) dist(x,method="euclidean"),
            main = graphTitle,
            trace="none",         # turns off trace lines inside the heat map
            margins=c(12,9),     # widens margins around plot
            dendrogram=dendroStatus,
            density.info="none",
            Rowv=rowV,
            Colv=colV,
            breaks=c(-7:-1/7,1:7/7),
            cexRow=1,cexCol=1,srtCol=45,
            col=heatmapColors))
  garbage <- dev.off()               # close the EPS device
  

  lhei = c(0.5,0.75,4)
  lwid = c(1.5,4)
  lmat = rbind(c(4,0),
	       c(0,3),
               c(2,1))
  pdf(filenames$pdf,height=17,width=8.5)
  garbageCollect <- capture.output(heatmap.2(heatmapData,
            hclustfun=function(x) hclust(x,method="ward.D2"),
            distfun=function(x) dist(x,method="euclidean"),
            main = graphTitle,
            trace="none",         # turns off trace lines inside the heat map
            margins=c(12,9),     # widens margins around plot
            dendrogram=dendroStatus,
            density.info="none",
            Rowv=rowV,
            Colv=colV,
            keysize=1,
            breaks=c(-7:-1/7,1:7/7),
            cexRow=1,cexCol=1,srtCol=45,
	    lmat = lmat, lwid = lwid, lhei = lhei,
            col=heatmapColors))
  garbage <- dev.off()               # close the PDF device
  
  if(capabilities("tiff")) {
    tiff(filenames$tiff)
    garbageCollect <- capture.output(heatmap.2(heatmapData,
              hclustfun=function(x) hclust(x,method="ward.D2"),
              distfun=function(x) dist(x,method="euclidean"),
              main = graphTitle,
              trace="none",         # turns off trace lines inside the heat map
              margins=c(12,9),     # widens margins around plot
              dendrogram=dendroStatus,
              density.info="none",
              Rowv=rowV,
              Colv=colV,
              breaks=c(-7:-1/7,1:7/7),
              cexRow=1,cexCol=1,srtCol=45,
              col=heatmapColors))
    garbage <- dev.off()               # close the TIFF device
  }
  
  suppressWarnings(xfig(filenames$xfig))
  garbageCollect <- capture.output(heatmap.2(heatmapData,
            hclustfun=function(x) hclust(x,method="ward.D2"),
            distfun=function(x) dist(x,method="euclidean"),
            main = graphTitle,
            trace="none",         # turns off trace lines inside the heat map
            margins=c(12,9),     # widens margins around plot
            dendrogram=dendroStatus,
            density.info="none",
            Rowv=rowV,
            Colv=colV,
            breaks=c(-7:-1/7,1:7/7),
            cexRow=1,cexCol=1,srtCol=45,
            col=heatmapColors))
  garbage <- dev.off()               # close the FIG device

  cat(paste0("All plots have been created at ./heatmap_output/",pathwaysDir,
             "\n"),sep="")
}
