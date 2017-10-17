## https://github.com/oganm/MSigDB
# devtools::install_github('oganm/MSigDB')
library(MSigDB)
names(MSigDB)

head(names(MSigDB$C2_CURATED))

head(MSigDB$HALLMARK$HALLMARK_COMPLEMENT)

getMSigInfo('HALLMARK_COMPLEMENT')

## Gene Synonym: https://github.com/oganm/geneSynonym
library(devtools)
install_github('oganm/geneSynonym')

library(geneSynonym)
humanSyno('MOG')

###
##
## KEGG pathway database API connection
## KEGGREST
## https://www.bioconductor.org/packages/devel/bioc/vignettes/KEGGREST/inst/doc/KEGGREST-vignette.html
##
## The KEGG REST API is built on some simple operations: info, list, find, get, conv, and link. 
## The corresponding R functions in KEGGREST are: 
##  keggInfo(), keggList(), keggFind(), keggGet(), keggConv, and keggLink().
##
library(KEGGREST)
listDatabases()
