# PIMENTo
PIMENTo is an R package to normalize, analyze, and visualize microarray data

# Native R Installation
You can install `PIMENTo` using `devtools` and `BiocManager`:
```R
# Install Dependencies from CRAN
install.packages(c('devtools', 'BiocManager'), repos='https://cloud.r-project.org/')

# Install Dependencies from Bioconductor
BiocManager::install(c('impute', 'DESeq2', 'affy', 'limma'))

# Install PIMENTo
devtools::install_github("HardimanLab/PIMENTo")

# Load package
library(PIMENTo)
```
# Run PIMENTo with Docker
To start an interactive `PIMENTo` session and mount the current <br/>
working directory, `$PWD`, within the `Docker` container to /home/pimento/ <br/>
Run the following:
```dockerfile
docker run -v "$PWD:/home/pimento/" --rm -it hardimanlab/pimento
```
<sub>*Note that the --rm flag means that the container will autoremove <br/>
after you close the PIMENTo session.*</sub> <br/>

# Links
[PIMENTo - GitHub](https://github.com/HardimanLab/PIMENTo)<br/>
[PIMENTo - DockerHub](https://hub.docker.com/r/hardimanlab/pimento)
