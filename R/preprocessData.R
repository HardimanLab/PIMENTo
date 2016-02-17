#' @title Generate MA plots of raw and normalized data
#' @description Normalize input raw data using quantile and mloess methods. Plots
#' of the normalized data along with a dendrogram librarclustering all samples will 
#' be stored in newly created pipeline directory.
#' @usage preprocessPlots(inputFile, fileSheet=1, ntext=2, dataCol, symbolIndex=1,
#'  idIndex=2)
#' @param inputFile Path to the microarray expression file, be it .xlsx or .csv
#' @param fileSheet Sheet number in the spreadsheet with data
#' @param ntext Number of leading text columns
#' @param dataCol Range of columns which contain data (indexing begins with
#' first column of file)
#' @param symbolIndex Column index which contains gene symbols
#' @param idIndex Column index which contains gene ID's
#' @return A list with components
#' \item{ntext}{Number of leading text columns}
#' \item{dataCol}{Vector of column indices containing array data}
#' \item{id}{Vector containing gene ID's}
#' \item{idIndex}{Column index containing gene ID information}
#' \item{symbol}{Vector containing gene symbols}
#' \item{symbolIndex}{Column index containing gene symbol information}
#' \item{descStats}{Vector of column indices containing descriptive statistics}
#' \item{pipelineName}{Name of pipeline generated from input file name sans 
#' extension}
#' \item{mloess}{Data rame of quantile normalized data}
#' \item{quantile}{Data rame of quantile normalized data}
#' @export

preprocessData <- function(inputFile,fileSheet=1,ntext=2,dataCol,symbolIndex=1,
                           idIndex=2) {

  ext <- tools::file_ext(inputFile)
  if (ext=="xlsx")
    X <- openxlsx::read.xlsx(inputFile,sheet=fileSheet)
  else if (ext=="csv")
    X <- read.csv(inputFile,header=TRUE,sep=",",stringsAsFactors=F)
  else if (ext=="txt")
    X <- read.table(inputFile,header=TRUE,sep="\t",stringsAsFactors=F)
  else
    stop("Input file must be a .xlsx spreadsheet, comma-separated .csv, or 
         tab-seperated .txt")
  
  pipelineName <- tools::file_path_sans_ext(inputFile)
  panels <- length(dataCol)
  format <- c(2,3)

  labels <- colnames(X)[dataCol]
  textCol <- 1:ntext  
  data <- X[,dataCol]
  
  if(ntext==1) {
    symbolIndex = 1
    idIndex = 2
    dataCol <- dataCol+1
    descStats <- X[,rep(textCol,2)]
  } else {
    descStats <- X[,textCol]
  }
  
  dataList = list()
  dataList$ntext=ntext
  dataList$dataCol=dataCol
  dataList$id=X[,idIndex]
  dataList$idIndex=idIndex
  dataList$symbol=X[,symbolIndex]
  dataList$symbolIndex=symbolIndex
  dataList$descStats=descStats
  dataList$pipelineName=pipelineName
  
  if (!dir.exists(paste0(getwd(),"/",pipelineName,"_pipeline"))) {
    cat("Creating output directory at ./",pipelineName,"_pipeline\n",sep="")
    dir.create(file.path(getwd(),paste0(pipelineName,"_pipeline")))
  }
  
  for (normType in c("raw","loess","quantile")) {
    ps.outputFile=paste0("./",pipelineName,"_pipeline/",pipelineName,"_",normType,".ps")
    pdf.outputFile=paste0("./",pipelineName,"_pipeline/",pipelineName,"_",normType,".pdf")

    if (normType == "raw") {
      dataNorm <- data
      ps.clusterFile=paste0("./",pipelineName,"_pipeline/",pipelineName,"_cluster.ps")
      pdf.clusterFile=paste0("./",pipelineName,"_pipeline/",pipelineName,"_cluster.pdf")
      
      dist <- dist(t(dataNorm))
      
      postscript(file=ps.clusterFile,paper="letter")
      plot(hclust(dist,method="ward.D2"))
      invisible(dev.off())
      
      pdf(file=pdf.clusterFile,paper="letter")
      plot(hclust(dist,method="ward.D2"))
      invisible(dev.off())
      
      cat("Dendrogram of raw data plotted at",ps.clusterFile,"\n")
      cat("Raw data plots created at",ps.outputFile,"\n")
    } else if (normType == "loess") {
      capture.output(dataNorm <- affy::normalize.loess(as.matrix(data)))
      dataList$mloess <- cbind(descStats,dataNorm)
      cat("MLOESS normalized data plots created at",ps.outputFile,"\n")
    } else if (normType == "quantile") {
      dataNorm <- limma::normalizeQuantiles(as.matrix(data))
      dataList$quantile <- cbind(descStats,dataNorm)
      cat("Quantile normalized data plots created at",ps.outputFile,"\n")
    }

    # Postscript output
    postscript(file=ps.outputFile,paper="letter")
    par(mfrow=format,pty="s")
        
    for (i in 1:panels) {     
      x<-as.numeric(dataNorm[,((i-1) %% panels)+1])
      y<-as.numeric(dataNorm[,(i %% panels)+1])
      A<-(log2(x)+log2(y))/2
      M<-log2(y)-log2(x)
      plot(A,M,xlim=c(4,14),ylim=c(-3,3),cex=1,pch=".",xlab=labels[((i-1) %% 
          panels)+1],main=labels[(i %% panels)+1],font.main=1,cex.main=1)
      lines(c(0,18),c(0,0),col="cyan")
    }
    invisible(dev.off()) 
    
    # PDF output
    pdf(file=pdf.outputFile,paper="letter")
    par(mfrow=format,pty="s")
    
    for (i in 1:panels) {     
      x<-as.numeric(dataNorm[,((i-1) %% panels)+1])
      y<-as.numeric(dataNorm[,(i %% panels)+1])
      A<-(log2(x)+log2(y))/2
      M<-log2(y)-log2(x)
      plot(A,M,xlim=c(4,14),ylim=c(-3,3),cex=1,pch=".",xlab=labels[((i-1) %% 
           panels)+1],main=labels[(i %% panels)+1],font.main=1,cex.main=1)
      lines(c(0,18),c(0,0),col="cyan")
    }
    invisible(dev.off()) 
    
  }
  return(dataList)
  
}
