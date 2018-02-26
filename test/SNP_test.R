#install.packages("rsnps")
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("org.Hs.eg.db")
# install.packages("HardyWeinberg")
# install.packages("powerGWASinteraction")

library(rsnps)

ncbi_snp_query('rs4845140')
genotypes('rs4845140', userid='1,6,8', df=TRUE)


library(rentrez)

gene_search = entrez_search(db   = "gene",
                            term = "IL10[Gene Name]) AND (Homo sapiens[Organism])", retmax = 1)
geneId <- gene_search$ids 

snp_links <- entrez_link(dbfrom='gene', id=geneId, db='snp')

# access results with $links
length(snp_links$links$gene_snp)

head(snp_links$links$gene_snp)

snp_links$links$gene_snp # $gene_clinvar
## with multi gene ids:
multi_snp_links <- entrez_link(dbfrom='gene', id=c("5728", "374654"), db='snp', by_id=TRUE) # by_id = True makes the IDs grouped by ids.
lapply(multi_snp_links, function(x) head(x$links$gene_snp))


## 
snp_search = entrez_search(db = 'snp',
                           term = "rs1800896[Reference SNP ID] OR rs1800871[Reference SNP ID] OR rs1800872[Reference SNP ID]")
snp_search$count
snp_ids <- snp_search$ids
gene_links <- entrez_link(dbfrom = 'snp', id =snp_ids, db = 'gene')
snpGene_ID <- gene_links$links$snp_gene
geneName <- entrez_search(db = 'gene', term = paste0(snpGene_ID, "[uid]"), 
                          rettype = 'gene_table')

geneName$ids
geneName$count
geneName$file
geneName$QueryTranslation

### get gene symbol from Entrez Gene ID
library(org.Hs.eg.db)
snpGene_name = unlist(mget(x=snpGene_ID, envir=org.Hs.egGENENAME))

### Power calculation
library("powerGWASinteraction")

n <- 1000
power <- 0.95
model <- list()
model$prev <- 0.001
model$pGene1 <- 0.05
model$pGene2 <- 0.04
model$nSNP <- 1000
model$beta.LOR <- c(0.4, 0.3, 0.3)
alpha1 <- 0.05
powerGG(n = n, 
        #power = power,
        model = model, alpha1 = alpha1)
###