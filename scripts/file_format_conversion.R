## install packages as needed - only need once
## utils/InstallPacakges.R
#####
source("utils/utility.R")
###################
# Step 1. Load raw data and preprocessing
###################
input.dir <- "../Data/"
output.dir <- "../R_Output/"

# input
gwas.fn <- lapply(c(bed = 'bed', bim='bim', fam='fam', gds='gds', map='map'), function(n) sprintf("%s/22bed.%s", input.dir, n))
clinical.fn <- sprintf("%scldt.csv", input.dir)
onethou.fn <- lapply(c(ped='ped', nosex='nosex'), function(x) sprintf("%s22bed.%s", input.dir, x))

# output
gwaa.fname <- sprintf("%s22bedOut.txt", output.dir)





## 
library(snpStats)
# read PLINK files to create a list
genotypes <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = c("-9"))
# genotype data
genotype <- genotypes$genotypes
print(genotype) # SnpMatrix with  1983 rows and  123536 columns

# map data
map <- genotypes$map
print(map)

snpsum <- col.summary(genotype)
summary(snpsum)

# by SNPs
par(mfrow = c(1, 2))
hist(snpsum$MAF)
hist(snpsum$z.HWE)

# by samples
sample.qc <- row.summary(genotype)
summary(sample.qc)

par(mfrow = c(1, 1))
plot(sample.qc)

## Analysis
as(genotype[fam$member[1:5], 101207], "character")

which(grepl("rs12157537*", colnames(genotype)))

# original input
original_gwas.fn <- lapply(c(bed = 'bed', bim='bim', fam='fam', gds='gds', map='map'), function(n) sprintf("%s/22_original.%s", input.dir, n))
#clinical.fn <- sprintf("%scldt.csv", input.dir)
onethou.fn <- lapply(c(ped='ped', nosex='nosex'), function(x) sprintf("%s22bed.%s", input.dir, x))

# read PLINK files to create a list
original_genotypes <- read.plink(original_gwas.fn$bed, original_gwas.fn$bim, original_gwas.fn$fam, na.strings = c("-9"))
# genotype data
original_genotype <- original_genotypes$genotypes
print(original_genotype) # Sn

as(original_genotype[fam$member[1:5], "rs12157537"], "character")



library(rentrez)
# available DBs
entrez_dbs()

# Brief description of what the database is
entrez_db_summary("snp")

# Set of search terms that can used with this database
entrez_db_searchable("snp")

## Search db: at least need db and term
snpInfo <- entrez_search(db   = "snp",
                         term = "rs148791235[RS]")
snpInfo$
#

###
# R interface to Python
library("reticulate")
# Sys.which("python3")
use_python("/usr/local/bin/python3")

# import SciPy (it will be automatically discovered in "r-reticulate")
scipy <- import("scipy")
#pandas <- import("pandas")


