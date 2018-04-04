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
genotypes <- read.plink(gwas.fn$bed, gwas.fn$bim, gwas.fn$fam, na.strings = ("-9"))
# genotype data
genotype <- genotypes$genotypes
print(genotype) # SnpMatrix with  1983 rows and  123536 columns

# SNP info
genoBim <- genotypes$map 
colnames(genoBim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(genoBim))

rm(genotypes) # remove original reads, to free up space

### outcome
clinical_avail.fn <- sprintf("%savailable_cases.csv", output.dir)
hla.fn <- sprintf("%savailable_cases_HLA.csv", output.dir)

Outcome_avail <- read.csv(clinical_avail.fn, stringsAsFactors = F)
HLA_tab <- read.csv(hla.fn, stringsAsFactors = F)

#####
# Step 2. SNP-level filtering
#####

snpsum.col <- col.summary(genotype) # Summarize rows or columns of a snp matrix
print(head(snpsum.col))
#                           Calls Call.rate Certain.calls       RAF        MAF         P.AA       P.AB      P.BB
# rs148791235:16059973:C:A  1100 0.5547151             1 0.9886364 0.01136364 0.0000000000 0.02272727 0.9772727
# rs201278103:16081207:G:A   994 0.5012607             1 0.9818913 0.01810865 0.0000000000 0.03621730 0.9637827
# 22:16149693:C:T            994 0.5012607             1 0.9884306 0.01156942 0.0000000000 0.02313883 0.9768612
# rs369746900:16220993:C:G  1003 0.5057993             1 0.9860419 0.01395813 0.0009970090 0.02592223 0.9730808
# 22:16284380:T:TA          1079 0.5441251             1 0.9712697 0.02873031 0.0009267841 0.05560704 0.9434662
# rs199952431:16287339:C:T  1036 0.5224407             1 0.9864865 0.01351351 0.0000000000 0.02702703 0.9729730
#                           z.HWE
# rs148791235:16059973:C:A  0.3812212
# rs201278103:16081207:G:A  0.5814547
# 22:16149693:C:T           0.3690273
# rs369746900:16220993:C:G -1.8458639
# 22:16284380:T:TA         -0.1193081
# rs199952431:16287339:C:T  0.4409172

# thresholds
callrate <- 0.95 # call rate
minor <- 0.01    # MAF

use_snp <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= callrate)
use_snp[is.na(use_snp)] <- FALSE
cat(ncol(genotype)-sum(use_snp), "SNPs will be removed due to low MAF or call rate.\n") # 32542 removed

#
genotype <- genotype[, use_snp]
snpsum.col <- snpsum.col[use_snp, ]
print(genotype) # 90994 SNPs left

################
# Step 3. Sample-level filtering
################
library(SNPRelate)
library(plyr)

snpsum.row <- row.summary(genotype) # sample statistics (call.rate, heterozygosity)
# Add the F stat (inbreeding coefficient) to the snpsum.row
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity * (ncol(genotype)) * Call.rate)
snpsum.row$hetF <- 1 - (hetObs/hetExp)

head(snpsum.row)

## thresholds
sampCall <- 0.95  # sample call rate cut-off
hetcutoff <- 0.1  # Inbreeding coefficient cut-off

use_sample <- with(snpsum.row, !is.na(Call.rate) & Call.rate > sampCall & abs(hetF) <= hetcutoff)
use_sample[is.na(use_sample)] <- FALSE
cat(nrow(genotype)-sum(use_sample), "subjects will be removed due to low sample call or inbreeding coefficient.\n")
# 169 subjects are removed

# passed subs
genotype <- genotype[use_sample, ]
Outcome_avail <- filter_outcome(clinical_avail.fn, genotype)

## Checking for relatedness
ld.threshold <- 0.2  # LD cut-off
kin.threshold <- 0.01 # kinship cut-off 

# Create gdf file, required for SNPRelate function
snpgdsBED2GDS(gwas.fn$bed, gwas.fn$fam, gwas.fn$bim, gwas.fn$gds)
genofile <- openfn.gds(gwas.fn$gds, readonly = F)

# Automatically added "-1" sample suffixes are removed
# gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
# gds.ids <- sub("-1", "", gds.ids)
# add.gdsn(genofile, "sample.id", gds.ids, replace = T)

# Prune SNPs for IBD analysis
set.seed(1000)
geno.sample.ids <- rownames(genotype)
snp.ids <- colnames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.threshold, 
                          sample.id = geno.sample.ids,
                          snp.id = snp.ids)
snpset.ibd <- unlist(snpSUB, use.names = F)
cat(length(snpset.ibd), "will be used in IBD analysis.\n") # 2566

### IBD
# Find IBD coefficients using Method of Moments procedure
ibd <- snpgdsIBDMoM(genofile, kinship = T,
                    sample.id = geno.sample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
ibdcoeff <- snpgdsIBDSelection(ibd)
head(ibdcoeff)

## candicates of relatedness
ibdcoeff <- ibdcoeff[ibdcoeff$kinship >= kin.threshold, ]

# iteratively remove samples with high kinship starting with the sample with the most pairing
related.samples <- NULL
while(nrow(ibdcoeff)>0){
  
  # count the number of occurrences of each and take the top one
  sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, 'x']
  cat("Removing samples", as.character(rm.sample), "too closely related to", sample.counts[1, 'freq'], 'other samples.\n')
  
  # remove from ibdcoeff adn add to the list
  ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample, ]
  related.samples <- c(as.character(rm.sample), related.samples)
}
genotype <- genotype[!(rownames(genotype) %in% related.samples), ]

# checking for ancestry 
# find PCA matrix
pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids, snp.id = snpset.ibd, num.thread = 1)

# create data frame of the first two principal components
pctab <- data.frame(sample.id = pca$sample.id,
                    PC1 = pca$eigenvect[, 1],
                    PC2 = pca$eigenvect[, 2], 
                    stringsAsFactors = F)
# Plot the first two principal components
plot(pctab$PC2, pctab$PC1, xlab = "Principal Component 2", ylab="Principal Component 1", 
     main = "Ancestry Plot")
