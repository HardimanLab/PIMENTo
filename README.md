# PIMENTo
PIMENTo is an R package to normalize, analyze, and visualize microarray data

# Installation
You can install `PIMENTo` using the `devtools` package as such:
```coffee
# Install Dependencies from CRAN
install.packages(c('devtools', 'BiocManager'), repos='https://cloud.r-project.org/')

# Install Dependencies from Bioconductor
BiocManager::install(c('impute', 'DESeq2', 'affy', 'limma'))

# Install PIMENTo
devtools::install_github("HardimanLab/PIMENTo")

# Load package
library(PIMENTo)
```
