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

if (!require(DESeq)){
  biocLite("DESeq") ## install packege only if it cannot be loaded (Because not installed yet)
  library(DESeq)
  
}

if (!require(TCC)){
  biocLite("TCC") ## install packege only if it cannot be loaded (Because not installed yet)
  library(TCC)
  
}


if (!require (ggplot2)){
  install.packages('ggplot2')
  library(ggplot2)
  
}
