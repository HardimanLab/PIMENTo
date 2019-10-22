# PIMENTo
PIMENTo is an R package to normalize, analyze, and visualize microarray data

# Installation
You can install `PIMENTo` using the `devtools` package as such:
```coffee
install.packages("devtools")
source("https://bioconductor.org/biocLite.R")
biocLite(c("affy", "limma", "DESeq2","samr","gplots"))
library(devtools)
devtools::install_github("HardimanLab/PIMENTo")
library(PIMENTo)
```
