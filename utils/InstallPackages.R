## install packages - only need once
source("http://bioconductor.org/biocLite.R")
biocLite("snpStats")
biocLite("SNPRelate")
biocLite("rtracklayer")
biocLite("biomaRt")
biocLite("GeneGeneInteR")

install.packages(c("plyr", "GenABEL", "LDheatmap", "doParallel", "ggplot2", "coin", "igraph", "devtools"))

library("devtools")
install_url("http://cran.r-project.org/src/contrib/Archive/postgwas/postgwas_1.11.tar.gz")

# MSigDB - pathway databases
## https://github.com/oganm/MSigDB
devtools::install_github('oganm/MSigDB')

## R Interface to Python
install.packages("reticulate")

