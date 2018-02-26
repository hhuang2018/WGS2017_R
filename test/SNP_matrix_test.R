#source("https://bioconductor.org/biocLite.R")
#biocLite("snpStats")
#biocLite("snpMatrix")
outputFp <- "../../HLI/GenotypeTables/"

### Load Genotype data
GeneName = "IL10"
load(paste0(outputFp, GeneName, "_DRpair_genotypes.RData"))

IL10_list <- GeneGenoType_List

GeneName = "IL10RB"
load(paste0(outputFp, GeneName, "_DRpair_genotypes.RData"))

IL10RB_list <- GeneGenoType_List

###
#IL10_samples <- colnames(IL10_list)[-(1:5)]
#IL10RB_samples <- colnames(IL10RB_list)[-(1:5)]
# 
## Convert to SNP matrix
geneMat <- function(GeneName, GENE_list){
  Samples <- colnames(GENE_list)[-(1:5)]
  numSamp <- length(Samples)
  numSNV <- dim(GENE_list)[1]
  
  Gene_SNP_numeric <- matrix(data = NA, nrow = numSamp, ncol = numSNV)
  GVHD_Label <- vector(mode = "character", length = numSamp)
  SNPName <- vector(mode = "character", length = numSNV)
  
  for(id in 1:numSamp){
    
    GVHD_Label[id] <- strsplit(Samples[id], "_")[[1]][1]
    
    for(jd in 1:numSNV){
      SNPName[jd] <- paste0(GeneName, "-", GENE_list$CHROM[jd], "-", GENE_list$POS[jd])
      DR_GT <- unlist(GENE_list[Samples[id]])[jd]
      D_GT <- strsplit(DR_GT, "-")[[1]][1]
      R_GT <- strsplit(DR_GT, "-")[[1]][2]
      if(is.na(D_GT) | is.na(R_GT)){
        Gene_SNP_numeric[id, jd] <- 4
      }else if(D_GT == R_GT){
        Gene_SNP_numeric[id, jd] <- 1
      }else{
        
        if(is.na(D_GT)){
          D_GT <- c(GENE_list$REF[jd], GENE_list$REF[jd])
        }else{
          D_GT <- strsplit(D_GT, '/')[[1]]
        }
        
        if(is.na(R_GT)){
          R_GT <- c(GENE_list$REF[jd], GENE_list$REF[jd])
        }else{
          R_GT <- strsplit(R_GT, '/')[[1]]
        }
        
        # same genotype --- 0  1
        # 1 mismatch ------ 1  2
        # 2 mismatches ---- 2  3
        Gene_SNP_numeric[id, jd] <- 3 - length(intersect(D_GT, R_GT)) 
        
      }
    }
    
  }
  
  colnames(Gene_SNP_numeric) <- SNPName
  rownames(Gene_SNP_numeric) <- Samples
  return(list(Gene_mat = Gene_SNP_numeric, Label = GVHD_Label))
}
## Donor and Recipient
DR_GenoType <- function(GeneName, GENE_list){
  Samples <- colnames(GENE_list)[-(1:5)]
  numSamp <- length(Samples)
  numSNV <- dim(GENE_list)[1]
  
  Gene_SNP_donor <- matrix(data = NA, nrow = numSamp, ncol = numSNV)
  Gene_SNP_recipient <- matrix(data = NA, nrow = numSamp, ncol = numSNV)
  GVHD_Label <- vector(mode = "character", length = numSamp)
  SNPName <- vector(mode = "character", length = numSNV)
  
  for(id in 1:numSamp){
    
    GVHD_Label[id] <- strsplit(Samples[id], "_")[[1]][1]
    
    for(jd in 1:numSNV){
      SNPName[jd] <- paste0(GeneName, "-", GENE_list$CHROM[jd], "-", GENE_list$POS[jd])
      DR_GT <- unlist(GENE_list[Samples[id]])[jd]
      D_GT <- strsplit(DR_GT, "-")[[1]][1]
      R_GT <- strsplit(DR_GT, "-")[[1]][2]
      if(is.na(D_GT)){
        Gene_SNP_donor[id, jd] <- 0
      }else{
        D_GTs <- strsplit(D_GT, "/")
        Gene_SNP_donor[id, jd] <- D_GT
      }
      if(is.na(R_GT)){
        Gene_SNP_recipient[id, jd] <- paste0(GENE_list$REF[jd], "/", GENE_list$REF[jd])
      }else{
        Gene_SNP_recipient[id, jd] <- R_GT
      }
    }
    
  }
  
  colnames(Gene_SNP_donor) <- SNPName
  rownames(Gene_SNP_donor) <- Samples
  
  colnames(Gene_SNP_recipient) <- SNPName
  rownames(Gene_SNP_recipient) <- Samples
  return(list(Donor = Gene_SNP_donor, Recipient = Gene_SNP_recipient, Label = GVHD_Label))
}
## 
#GeneName <- "IL10"
#GENE_list <- IL10_list
IL10_Gene_List <- geneMat("IL10", IL10_list)

IL10_rm_ids <- vector(mode = "numeric", length=0)
for(id in 1:dim(IL10_Gene_List$Gene_mat)[2]){
  if(length(which(IL10_Gene_List$Gene_mat[,id] == 4)) > 20){
    IL10_rm_ids <- c(IL10_rm_ids, id)
  }
}

filtered_IL10_SNV <- IL10_Gene_List$Gene_mat[, -rm_ids]

IL10_Gene_List_GT <- DR_GenoType("IL10", IL10_list)
IL10_Donor <- IL10_Gene_List_GT$Donor
IL10_Recipient <- IL10_Gene_List_GT$Recipient


##
IL10RB_Gene_List <- geneMat("IL10RB", IL10RB_list)

IL10RB_rm_ids <- vector(mode = "numeric", length=0)
for(id in 1:dim(IL10RB_Gene_List$Gene_mat)[2]){
  if(length(which(IL10RB_Gene_List$Gene_mat[,id] == 4)) > 20){
    IL10RB_rm_ids <- c(IL10RB_rm_ids, id)
  }
}

filtered_IL10RB_SNV <- IL10RB_Gene_List$Gene_mat[, -rm_ids]

IL10RB_Gene_List_GT <- DR_GenoType("IL10RB", IL10RB_list)
IL10RB_Donor <- IL10RB_Gene_List_GT$Donor
IL10RB_Recipient <- IL10RB_Gene_List_GT$Recipient

##
Gene1_Recipient <- new("SnpMatrix", IL10_Recipient)
Gene2_Donor <- new("SnpMatrix", IL10RB_Donor)
##
Gene1_filtered <- new("SnpMatrix", filtered_IL10_SNV)
Gene2_filtered <- new("SnpMatrix", filtered_IL10RB_SNV)
##
Gene1 <- new("SnpMatrix", IL10_Gene_List$Gene_mat)
Gene2 <- new("SnpMatrix", IL10RB_Gene_List$Gene_mat)

Label <- as.factor(IL10_Gene_List$Label)
######
# Multidimensional methods at the gene level
# 1. PCA
PCA.test(Y=Label, G1=Gene1_Recipient, G2=Gene2_Donor,threshold=0.7,
         method="GenFreq")
PCA.test(Y=gene.pair$Y, G1=Gene1_Recipient, G2=Gene2_Donor, threshold=0.7,
         method="Std")

# 2.  Canonical Correlation Analysis (CCA)
set.seed(1234)
CCA.test(Y=gene.pair$Y, G1=gene.pair$G1, G2=gene.pair$G2,n.boot=500)

# 3.  Kernel Canonical Correlation Analysis (KCCA)
set.seed(1234)
KCCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,
          kernel="rbfdot",sigma = 0.05,n.boot=500)

KCCA.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,
          kernel="polydot",degree = 1, scale = 1, offset = 1)

# 4. Partial Least Square Path Modeling (PLSPM)
set.seed(1234)
PLSPM.test(Y=gene.pair$Y, G1=gene.pair$G1,G2=gene.pair$G2,n.perm=1000)

# 5. Composite Linkage Disequilibrium (CLD)
CLD.test(Y=Label, G1=Gene1_Recipient,G2=Gene2_Donor, n.perm=1000)

# 6. Gene-Based Information Gain Method (GBIGM)
GBIGM.test(Y=Label, G1=Gene1,G2=Gene2,n.perm=2000)


##  From SNP-SNP interaction to Gene-Gene interaction testing
# 1. minP
set.seed(1234)
minP.test(Y=Label, G1=Gene1_filtered,G2=Gene2_filtered)

# 2. GATES
set.seed(1234)
gates.test(Y=Label, G1=Gene1_filtered,G2=Gene2_filtered,me.est="ChevNy")

set.seed(1234)
gates.test(Y=Label, G1=Gene1_filtered,G2=Gene2_filtered,alpha=0.05,me.est="Keff")

set.seed(1234)
gates.test(Y=Label, G1=Gene1_filtered,G2=Gene2_filtered,me.est="LiJi")

set.seed(1234)
gates.test(Y=Label, G1=Gene1_filtered,G2=Gene2_filtered,me.est="Galwey")

# 3. tTS and tProd
set.seed(1234)
tTS.test(Y=Label, G1=Gene1_filtered,G2=Gene2_filtered,tau=0.5,n.sim=10000)

set.seed(1234)
tProd.test(Y=Label, G1=Gene1_filtered,G2=Gene2_filtered,tau=0.05,n.sim=1000)


