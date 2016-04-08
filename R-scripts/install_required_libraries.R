################################################################
## Install libraries required for this project

source("https://bioconductor.org/biocLite.R")

if (!require(edgeR)) {
  biocLite("edgeR") ## Install package only if it cannot be loaded (because not installed yet)
  library(edgeR)
}

if (!require(DESeq2)) {
  biocLite("DESeq2") ## Install package only if it cannot be loaded (because not installed yet)
  library(DESeq2)
}



