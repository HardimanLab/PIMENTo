#' @title Identify significant genes through SAM
#' @description Implement SAM and compute significant genes given delta. Output
#' will consist of all significant genes ordered by increasing q-value and 
#' decreasing d-score.
#' @usage runSAM(backgroundSub.obj, classCompareCols, classCompareName, 
#' fdr.cutoff=0.1, response)
#' @param backgroundSub.obj Object returned from call to 
#' backgroundSub
#' @param classCompareCols Vector of column indices indicating which subset 
#' of arrays are to be compared for this comparison
#' @param classCompareName String title given to the name of the comparison
#' @param fdr.cutoff Max FDR for SAM, will use delta value which results in max 
#' FDR below this cutoff
#' @param response Vector of 1, 2 values that indicate group membership
#' @return A list with components
#' \item{siggenesTable}{Combined data frame of genes having significant 
#' positive and negative correlation}
#' \item{dataCol}{Vector of column indices containing array data}
#' \item{ntext}{Number of leading text columns}
#' \item{response}{Vector of array group membership, 1=control, 2=experimental}
#' \item{pipelineName}{Name of pipeline generated from input file name sans 
#' extension}
#' \item{data}{Data frame of chosen normalization method data}
#' \item{classCompareCols}{Value entered through classCompareCols parameter}
#' \item{classCompareName}{Value entered through classCompareName parameter}
#' \item{symbolIndex}{Column index that contains gene symbol}
#' @export

runSAM <- function(backgroundSub.obj,classCompareCols,classCompareName,
                   fdr.cutoff=0.1,response) {
  
  if ((missing(classCompareCols) & !missing(classCompareName)) | 
       (missing(classCompareName) & !missing(classCompareCols))) {
    stop("Cannot have classCompareCols set without classCompareName
           and vice-versa")
  }

  if (missing(classCompareCols)) {
    dataSAM <- backgroundSub.obj$data[,backgroundSub.obj$dataCol]
  }
  else {
      dataSAM <- backgroundSub.obj$data[,classCompareCols]
  }
  log.dataSAM <- log2(dataSAM)
  genenames <- as.data.frame(backgroundSub.obj$symbol)
  geneid <- as.data.frame(backgroundSub.obj$id)
  
  if(length(response) != ncol(dataSAM)) {
    stop("Number of responses does not match number of samples.")
  }
  
  listSAM = list(x=log.dataSAM,y=response,genenames=genenames,geneid=geneid,
                 logged2=T)
  cat("Beginning SAM processing\n")
  
  if(length(unique(response)) == 2) {
    capture.output(samr.obj <- samr::samr(listSAM,resp.type="Two class unpaired",
                                          s0.perc=50,testStatistic="standard",
                                          nperms=200))
  } else {
    stop("Arrays can only belong to control (1) or experimental (2).")
  }
  
  cat("Calculating delta table\n")
  capture.output(delta.table <- samr::samr.compute.delta.table(samr.obj,nvals=1000))
  
  delta <- delta.table[which(delta.table[,5] <= fdr.cutoff)[1],1]
  while (is.na(delta)) {
    fdr.cutoff <- fdr.cutoff + 0.05
    if (fdr.cutoff == 1.00) {
      stop("Have reached cutoff of 1.00 and no delta found.")
    }
    cat("Cutoff is too stringent, no delta available. Increasing FDR cutoff to ",
        fdr.cutoff, "\n")
    delta <- delta.table[which(delta.table[,5] <= fdr.cutoff)[1],1]
  }
  desc.dataSAM <- data.frame(backgroundSub.obj$descStats,dataSAM)
  
  # samr.compute.siggenes.table flips genename and geneid
  colnames(desc.dataSAM)[c(backgroundSub.obj$symbolIndex,
                           backgroundSub.obj$idInd)] <- c("geneid","genenames")
  siggenes.table <- samr::samr.compute.siggenes.table(samr.obj,delta,desc.dataSAM,
                                                delta.table,compute.localfdr=TRUE)
  
  cat("\nSAM Results:\nOptimal delta:",delta,"\nNo. of significantly down regulated genes:",
      siggenes.table$ngenes.lo,"\nNo. of significantly up regulated genes:",
      siggenes.table$ngenes.up,"\n\n")
  
  if (sum(nrow(siggenes.table$genes.up),nrow(siggenes.table$genes.lo)) < 1) {
    stop("No significant genes at provided delta")
  }
    
  if (!missing(classCompareName)) {
    siggenesFile=paste0(backgroundSub.obj$pipelineName,"_pipeline/",
                        backgroundSub.obj$pipelineName,"_sigGenes-",
                        classCompareName,".csv")
  }
  else {
    siggenesFile=paste0(backgroundSub.obj$pipelineName,"_pipeline/",
                        backgroundSub.obj$pipelineName,"_sigGenes.csv")
  }
    
  allSiggenes <- as.matrix(rbind(siggenes.table$genes.up,siggenes.table$genes.lo))
  ordered.allSiggenes <- allSiggenes[order(-(as.numeric(allSiggenes[,8])),
                                             abs(as.numeric(allSiggenes[,4])),
                                             decreasing=TRUE),]
  write.siggenes <- as.data.frame(ordered.allSiggenes)
  write.siggenes$optimalDelta <- as.character(
    c(delta, rep(" ", nrow(ordered.allSiggenes)-1)))
  write.siggenes$fdr.cutoff <- as.character(
    c(fdr.cutoff, rep(" ", nrow(ordered.allSiggenes)-1)))
  write.csv(write.siggenes,file=siggenesFile,row.names=FALSE)
  cat("Significant gene list available at ./",siggenesFile,"\n",sep="")
  
  colnames(ordered.allSiggenes) <- make.names(colnames(allSiggenes))

  if (missing(classCompareCols)) {
    sam.return.list <- list(siggenesTable=ordered.allSiggenes,data=desc.dataSAM,
              ntext=backgroundSub.obj$ntext,
              pipelineName=backgroundSub.obj$pipelineName,
              response=response,dataCol=backgroundSub.obj$dataCol)
  }
  else {
    ntext <- backgroundSub.obj$ntext
    subsetClassCompareCols <- c((ntext+1):(ntext+length(classCompareCols)))
    sam.return.list <- list(siggenesTable=ordered.allSiggenes,data=desc.dataSAM,
              ntext=backgroundSub.obj$ntext,
              pipelineName=backgroundSub.obj$pipelineName,
              response=response,dataCol=backgroundSub.obj$dataCol,
              classCompareCols=subsetClassCompareCols,idIndex=backgroundSub.obj$idIndex,
              symbolIndex=backgroundSub.obj$symbolIndex,
              classCompareName=classCompareName)
  }
  sampleSimilarity(sam.return.list)

  return(sam.return.list)
}
