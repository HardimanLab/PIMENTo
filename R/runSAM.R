#' @title Identify significant genes through SAM
#' @description Implement SAM and compute significant genes given delta. Output
#' will consist of all significant genes ordered by increasing q-value and 
#' decreasing d-score.
#' @usage runSAM(backgroundSub.obj, response, delta)
#' @param backgroundSub.obj Object returned from call to 
#' backgroundSub
#' @param response Vector of array group membership, 1=control, 2=experimental
#' @param delta Tuning parameter to obtain significant number of genes
#' @return A list with components
#' \item{siggenesTable}{Combined data frame of genes having significant 
#' positive and negative correlation}
#' \item{dataCol}{Vector of column indices containing array data}
#' \item{ntext}{Number of leading text columns}
#' \item{response}{Vector of array group membership, 1=control, 2=experimental}
#' \item{pipelineName}{Name of pipeline generated from input file name sans 
#' extension}
#' \item{data}{Data frame of chosen normalization method data}
#' @export

runSAM <- function(backgroundSub.obj,response,delta) {
   
  dataSAM <- backgroundSub.obj$data[,backgroundSub.obj$dataCol]
  log.dataSAM <- log2(dataSAM)
  genenames <- as.data.frame(backgroundSub.obj$symbol)
  geneid <- as.data.frame(backgroundSub.obj$id)
  
  if(length(response) != ncol(dataSAM))
    stop("Number of responses does not match number of samples.")
  if(missing(delta))
    stop("Provide a delta tuning parameter")
  
  listSAM = list(x=log.dataSAM,y=response,genenames=genenames,geneid=geneid,
                 logged2=T)
  cat("Beginning SAM processing\n")
  
  if(length(unique(response)) == 2) {
    capture.output(samr.obj <- samr::samr(listSAM,resp.type="Two class unpaired",
                                          s0.perc=50,testStatistic="standard",
                                          nperms=200))
  } else if (length(unique(response)) > 2) {
    capture.output(samr.obj <- samr::samr(listSAM,resp.type="Multiclass",
                                          s0.perc=50,testStatistic="standard",
                                          nperms=200))
  } else {
    stop("Arrays can only belong to control (1) or experimental (2).")
  }
  cat("Calculating delta table\n")
  capture.output(delta.table <- samr::samr.compute.delta.table(samr.obj,nvals=200))
  
  desc.dataSAM <- backgroundSub.obj$data
  # samr.compute.siggenes.table flips genename and geneid
  colnames(desc.dataSAM)[c(backgroundSub.obj$symbolIndex,
                           backgroundSub.obj$idInd)] <- c("geneid","genenames")
  print(head(desc.dataSAM))
  siggenes.table <- samr::samr.compute.siggenes.table(samr.obj,delta,desc.dataSAM,
                                                delta.table,compute.localfdr=TRUE)
  
  cat("\nSAM Results:\nDelta:",delta,"\nNo. of significantly down regulated genes:",
      siggenes.table$ngenes.lo,"\nNo. of significantly up regulated genes:",
      siggenes.table$ngenes.up,"\n\n")
  if (sum(nrow(siggenes.table$genes.up),nrow(siggenes.table$genes.lo)) < 1)
    stop("No significant genes at provided delta")
  
  siggenesFile=paste0(backgroundSub.obj$pipelineName,"_pipeline/",
                      backgroundSub.obj$pipelineName,"_sigGenes.csv")
  allSiggenes <- as.matrix(rbind(siggenes.table$genes.up,siggenes.table$genes.lo))
  ordered.allSiggenes <- allSiggenes[order(-(as.numeric(allSiggenes[,8])),
                                             abs(as.numeric(allSiggenes[,4])),
                                             decreasing=TRUE),]
  write.csv(ordered.allSiggenes,file=siggenesFile,row.names=FALSE)
  cat("Significant gene list available at ./",siggenesFile,sep="")
  
  colnames(ordered.allSiggenes) <- make.names(colnames(allSiggenes))
    
  return(list(siggenesTable=ordered.allSiggenes,data=desc.dataSAM,
              ntext=backgroundSub.obj$ntext,
              pipelineName=backgroundSub.obj$pipelineName,
              response=response,dataCol=backgroundSub.obj$dataCol))
}