#' @title Identify significant genes through SAM
#' @description Implement SAM and compute significant genes given delta. Output
#' will consist of all significant genes ordered by increasing q-value and 
#' decreasing d-score.
#' @usage SigGenesSAM(background.subtraction.obj, class.compare.cols, 
#' class.compare.name, fdr.cutoff=0.1, response)
#' @param background.subtraction.obj Object returned from call to 
#' BackgroundSubtraction
#' @param class.compare.cols Vector of column indices indicating which subset 
#' of arrays are to be compared for this comparison
#' @param class.compare.name String title given to the name of the comparison
#' @param fdr.cutoff Max FDR for SAM, will use delta value which results in max 
#' FDR below this cutoff
#' @param response For two class unpaired: vector of 1, 2 values that indicate 
#' group membership. For two class paired: vector of -1, 1, -2, 2, etc. 
#' values that indicate pairings.
#' @return A list with components
#' \item{siggenes.table}{Combined data frame of genes having significant 
#' positive and negative correlation}
#' \item{data.col}{Vector of column indices containing array data}
#' \item{ntext}{Number of leading text columns}
#' \item{response}{Vector of array group membership, 1=control, 2=experimental}
#' \item{pipeline.name}{Name of pipeline generated from input file name sans 
#' extension}
#' \item{data}{Data frame of chosen normalization method data}
#' \item{class.compare.cols}{Value entered through class.compare.cols parameter}
#' \item{class.compare.name}{Value entered through class.compare.name parameter}
#' \item{symbol.index}{Column index that contains gene symbol}
#' @export

SigGenesSAM <- function(background.subtraction.obj, class.compare.cols, 
                        class.compare.name, fdr.cutoff=0.1, response, delta) {
  
  if ((missing(class.compare.cols) & !missing(class.compare.name)) | 
       (missing(class.compare.name) & !missing(class.compare.cols))) {
    stop("Cannot have class.compare.cols set without class.compare.name
           and vice-versa")
  }

  if (missing(class.compare.cols)) {
    data.SAM <- 
      background.subtraction.obj$normalized[, background.subtraction.obj$data.col]
  }
  else {
      data.SAM <- background.subtraction.obj$normalized[, class.compare.cols]
  }
  log.data.SAM <- log2(data.SAM)
  
  genenames <- as.data.frame(background.subtraction.obj$symbol)
  geneid <- as.data.frame(background.subtraction.obj$id)
  
  symbol.index <- background.subtraction.obj$symbol.index
  geneid.index <- background.subtraction.obj$id.index
  pipeline.name <- background.subtraction.obj$pipeline.name
  
  if(length(response) != ncol(data.SAM)) {
    stop("Number of responses does not match number of samples.")
  }
  
  list.SAM = list(x=log.data.SAM, y=response, genenames=genenames, 
                  geneid=geneid, logged2=T)
  
  list.SAM$x <- as.matrix(list.SAM$x)
  
  cat("Beginning SAM processing\n")
  
  if (sum(response < 0) == 0) {
    capture.output(samr.obj <- samr::samr(list.SAM, 
                                          resp.type="Two class unpaired",
                                          s0.perc=50, testStatistic="standard",
                                          nperms=200))
  } else {
    capture.output(samr.obj <- samr::samr(list.SAM, 
                                          resp.type="Two class paired",
                                          s0.perc=50, testStatistic="standard",
                                          nperms=200))
  }  

  
  if(missing(delta)) {
    cat("Calculating delta table\n")
    capture.output(delta.table <- samr::samr.compute.delta.table(samr.obj, 
                                                                 nvals=1000))
    
    delta <- delta.table[which(delta.table[, 5] <= fdr.cutoff)[1], 1]
    while (is.na(delta)) {
      fdr.cutoff <- fdr.cutoff + 0.05
      if (fdr.cutoff == 1.00) {
        stop("Have reached cutoff of 1.00 and no delta found.")
      }
      cat("Cutoff is too stringent, no delta available. Increasing FDR cutoff to ",
          fdr.cutoff, "\n")
      delta <- delta.table[which(delta.table[, 5] <= fdr.cutoff)[1], 1]
    }
  }
  desc.data.SAM <- data.frame(background.subtraction.obj$desc.stats, data.SAM)
  
  # samr.compute.siggenes.table flips genename and geneid
  symbol.id.columns <- c(symbol.index, geneid.index)
  colnames(desc.data.SAM)[symbol.id.columns] <- c("geneid", "genenames")
  siggenes.table <- samr::samr.compute.siggenes.table(samr.obj, delta, 
                                                      desc.data.SAM,
                                                      delta.table, 
                                                      compute.localfdr=T)
  
  cat("\nSAM Results:\nOptimal delta:", delta, 
      "\nNo. of significantly down regulated genes:",
      siggenes.table$ngenes.lo, "\nNo. of significantly up regulated genes:",
      siggenes.table$ngenes.up, "\n\n")
  
  if (sum(nrow(siggenes.table$genes.up), nrow(siggenes.table$genes.lo)) < 1) {
    stop("No significant genes at provided delta")
  }
    
  if (!missing(class.compare.name)) {
    siggenes.file=paste0(pipeline.name, "_pipeline/",
                         pipeline.name, "_siggenes_",
                         class.compare.name, ".csv")
  }
  else {
    siggenes.file=paste0(pipeline.name, "_pipeline/",
                         pipeline.name, "_siggenes.csv")
  }
    
  all.siggenes <- as.matrix(rbind(siggenes.table$genes.up, siggenes.table$genes.lo))
  ordered.all.siggenes <- all.siggenes[order(-(as.numeric(all.siggenes[, 8])),
                                             abs(as.numeric(all.siggenes[, 4])),
                                             decreasing=T), ]
  write.siggenes <- as.data.frame(ordered.all.siggenes)
  write.siggenes$optimalDelta <- as.character(
    c(delta, rep(" ", nrow(ordered.all.siggenes)-1)))
  write.siggenes$fdr.cutoff <- as.character(
    c(fdr.cutoff, rep(" ", nrow(ordered.all.siggenes)-1)))
  write.csv(write.siggenes[,-1], file=siggenes.file, row.names=F)
  cat("Significant gene list available at ./", siggenes.file, "\n", sep="")
  
  colnames(ordered.all.siggenes) <- make.names(colnames(all.siggenes))
  
  if (missing(class.compare.cols)) {
    sam.return.list <- list(siggenes.table=ordered.all.siggenes, 
                            normalized=desc.data.SAM,
                            ntext=background.subtraction.obj$ntext,
                            pipeline.name=pipeline.name,
                            response=response,
                            data.col=background.subtraction.obj$data.col,
                            id.index=geneid.index,
                            symbol.index=symbol.index,
                            fdr.cutoff=fdr.cutoff)
    save(background.subtraction.obj, sam.return.list, 
         file=paste0(pipeline.name, "_pipeline/", "PIMENTo-", pipeline.name, "_", 
                     format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".RData"))
  }
  else {
    ntext <- background.subtraction.obj$ntext
    subsetclass.compare.cols <- c((ntext+1):(ntext+length(class.compare.cols)))
    sam.return.list <- list(siggenes.table=ordered.all.siggenes, 
                            normalized=desc.data.SAM,
                            ntext=background.subtraction.obj$ntext, 
                            pipeline.name=pipeline.name,
                            response=response, 
                            data.col=background.subtraction.obj$data.col, 
                            class.compare.cols=subsetclass.compare.cols, 
                            class.compare.name=class.compare.name,
                            id.index=geneid.index,
                            symbol.index=symbol.index,
                            fdr.cutoff=fdr.cutoff)
    
    save(background.subtraction.obj, sam.return.list, 
         file=paste0(pipeline.name, "_pipeline/", "PIMENTo-", pipeline.name, 
                     "_", class.compare.name, "_", 
                     format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".RData"))
  }
  SampleSimilarity(sam.return.list)

  return(sam.return.list)
}
