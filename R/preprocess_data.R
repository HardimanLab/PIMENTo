#' @title Generate MA plots of raw and normalized data
#' @description Normalize input raw data using quantile and mloess methods. 
#' Plots of the normalized data along with a dendrogram clustering all samples 
#' will be stored in newly created pipeline directory.
#' @usage preprocessPlots(input.file, file.sheet=1, ntext=2, data.col, 
#' symbol.index=1, id.index=2)
#' @param input.file Path to the microarray expression file, be it .xlsx or .csv
#' @param file.sheet Sheet number in the spreadsheet with data
#' @param ntext Number of leading text columns
#' @param data.col Range of columns which contain data (indexing begins with
#' first column of file)
#' @param symbol.index Column index which contains gene symbols
#' @param id.index Column index which contains gene ID's
#' @param batch.vector Character vector indicating to which batch each sample
#' belongs
#' @return A list with components
#' \item{ntext}{Number of leading text columns}
#' \item{data.col}{Vector of column indices containing array data}
#' \item{id}{Vector containing gene ID's}
#' \item{id.index}{Column index containing gene ID information}
#' \item{symbol}{Vector containing gene symbols}
#' \item{symbol.index}{Column index containing gene symbol information}
#' \item{desc.stats}{Vector of column indices containing descriptive statistics}
#' \item{pipeline.name}{Name of pipeline generated from input file name sans 
#' extension}
#' \item{mloess}{Data rame of quantile normalized data}
#' \item{quantile}{Data rame of quantile normalized data}
#' @export

PreprocessData <- function(input.file, file.sheet=1, ntext=2, data.col, 
                           symbol.index=1, id.index=2, batch.vector=NA) {

  ext <- tools::file_ext(input.file)
  if (ext == "xlsx")
    X <- openxlsx::read.xlsx(input.file, sheet=file.sheet)
  else if (ext == "csv")
    X <- read.csv(input.file, header=TRUE, sep=",", stringsAsFactors=F)
  else if (ext == "txt")
    X <- read.table(input.file, header=TRUE, sep="\t", stringsAsFactors=F)
  else
    stop("Input file must be a .xlsx spreadsheet, comma-separated .csv, or 
         tab-seperated .txt")
  
  pipeline.name <- tools::file_path_sans_ext(input.file)
  panels <- length(data.col)
  format <- c(2, 3)
  remove.batch <- F

  if (!is.na(batch.vector)) {
    if (!length(batch.vector == length(data.col))) {
      stop("Length of batch vector must equal number of data columns.")
    }
    remove.batch <- T
    X <- limma::removeBatchEffect(X, batch=batch.vector)
  }
  
  labels <- colnames(X)[data.col]
  text.col <- 1:ntext  
  data <- X[, data.col]
  
  if(ntext == 1) {
    symbol.index = 1
    id.index = 2
    data.col <- data.col + 1
    desc.stats <- X[, rep(text.col, 2)]
  } else {
    desc.stats <- X[, text.col]
  }
  
  data.list = list()
  data.list$ntext=ntext
  data.list$data.col=data.col
  data.list$id=X[, id.index]
  data.list$id.index=id.index
  data.list$symbol=X[, symbol.index]
  data.list$symbol.index=symbol.index
  data.list$desc.stats=desc.stats
  data.list$pipeline.name=pipeline.name
  
  if (!dir.exists(paste0(getwd(), "/", pipeline.name, "_pipeline"))) {
    cat("Creating output directory at ./", pipeline.name, "_pipeline\n", sep="")
    dir.create(file.path(getwd(), paste0(pipeline.name, "_pipeline")))
  }
  
  for (normType in c("raw", "mloess", "quantile")) {
    if (!remove.batch) {
      ps.plots <- paste0("./", pipeline.name, "_pipeline/", pipeline.name, "_", 
                         normType, "_plots.ps")
      pdf.plots <- paste0("./", pipeline.name, "_pipeline/", pipeline.name, "_", 
                          normType, "_plots.pdf")
    } else {
      ps.plots <- paste0("./", pipeline.name, "_pipeline/", pipeline.name, "_", 
                         normType, "_remove_batch_plots.ps")
      pdf.plots <- paste0("./", pipeline.name, "_pipeline/", pipeline.name, "_", 
                          normType, "_remove_batch_plots.pdf")
    }

    if (normType == "raw") {
      data.norm <- data      
      dist <- dist(t(data.norm))
      
      ps.cluster <- paste0("./", pipeline.name, "_pipeline/", pipeline.name,
                           "_cluster.ps")
      pdf.cluster <- paste0("./", pipeline.name, "_pipeline/", pipeline.name, 
                            "_cluster.pdf")

      postscript(file=ps.cluster, paper="letter")
      plot(hclust(dist, method="ward.D2"))
      garbage <- dev.off()
      
      pdf(file=pdf.cluster, paper="letter")
      plot(hclust(dist, method="ward.D2"))
      garbage <- dev.off()
      
      cat("Dendrogram of raw data plotted at", ps.cluster, "\n")
      cat("Raw data plots created at", ps.plots, "\n")
      
    } else if (normType == "mloess") {
      mat <- as.matrix(data)
      mat[mat == 0] <- 1

      capture.output(data.norm <- affy::normalize.loess(mat))

      data.list$mloess <- cbind(desc.stats, data.norm)
      write.csv(data.list$mloess, paste0("./", pipeline.name, 
                                        "_pipeline/", pipeline.name, 
                                        "_mloess_data.csv"), row.names=F)
      cat("MLOESS normalized data plots created at", ps.plots, "\n")
      
    } else if (normType == "quantile") {
      data.norm <- limma::normalizeQuantiles(as.matrix(data))

      data.list$quantile <- cbind(desc.stats, data.norm)
      write.csv(data.list$quantile, paste0("./", pipeline.name, 
                                          "_pipeline/", pipeline.name, 
                                          "_quantile_data.csv"), row.names=F)
      cat("Quantile normalized data plots created at", ps.plots, "\n")
    }

    # Postscript output
    postscript(file=ps.plots, paper="letter")
    par(mfrow=format, pty="s")
        
    for (i in 1:panels) {     
      x <- as.numeric(data.norm[, ((i-1) %% panels) + 1])
      y <- as.numeric(data.norm[, (i %% panels) + 1])
      A <- (log2(x) + log2(y)) / 2
      M <- log2(y) - log2(x)
      plot(A, M, xlim=c(4, 14), ylim=c(-3, 3), cex=1, pch=".", 
           xlab=labels[((i - 1) %% panels) + 1], 
           main=labels[(i %% panels) + 1], font.main=1, cex.main=1)
      lines(c(0, 18), c(0, 0), col="cyan")
    }
    garbage <- dev.off() 
    
    # PDF output
    pdf(file=pdf.plots, paper="letter")
    par(mfrow=format, pty="s")
    
    for (i in 1:panels) {     
      x <- as.numeric(data.norm[, ((i - 1) %% panels) + 1])
      y <- as.numeric(data.norm[, (i %% panels) + 1])
      A <- (log2(x) + log2(y)) / 2
      M <- log2(y) - log2(x)
      plot(A, M, xlim=c(4, 14), ylim=c(-3, 3), cex=1, pch=".", 
           xlab=labels[((i - 1) %% panels) + 1], 
           main=labels[(i %% panels) + 1], font.main=1, cex.main=1)
      lines(c(0, 18), c(0, 0), col="cyan")
    }
    garbage <- dev.off() 
  }
  return(data.list)
}