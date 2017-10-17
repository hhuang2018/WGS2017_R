################# All HLI cohort All MiHA
Individual_MiHA_table_fp <- "../Output/wwangMiHAIP_DtoR/RestrictedMiHA/"
Individual_MiHA_table_files <- list.files(Individual_MiHA_table_fp, pattern = "\\.RData")
num_files <- length(Individual_MiHA_table_files)

# save(Restircted_MiHA_db, file = "../Data/Restricted_MiHA_database_RMduplicated.RData")

load("../Data/Restricted_MiHA_database_RMduplicated.RData") ## Restircted_MiHA_db
load("../Data/ID_table_wCaseID.RData")          ## ID_table
load("../Data/HLI_reformatted_HLA_table_corrected.RData") ## reformated_HLA_typing_list

num_Restricted_MiHA <- dim(Restircted_MiHA_db)[1]

aGVHD_Restricted_MiHA_table <- data.frame(MiHA_SNP = Restircted_MiHA_db$HLA_SNP, stringsAsFactors = F)
aGVHD_Restricted_MiHA_table <- cbind(aGVHD_Restricted_MiHA_table, matrix(0, nrow = num_Restricted_MiHA, ncol = 32))
aGVHD_Restricted_MiHA_table_count <- data.frame(HLA_SNP = Restircted_MiHA_db$HLA_SNP, 
                                                SNP =  Restircted_MiHA_db$rs.number, 
                                                Count = numeric(length(Restircted_MiHA_db$HLA_SNP)),
                                                stringsAsFactors = F)
nGVHD_Restricted_MiHA_table <- aGVHD_Restricted_MiHA_table
nGVHD_Restricted_MiHA_table_count <- aGVHD_Restricted_MiHA_table_count

for(id in 1:num_files){
  
  eval(parse(text=paste0("load (\"", Individual_MiHA_table_fp, Individual_MiHA_table_files[id], "\")"))) # MiHA_table
  
  groupID <- as.numeric(unlist(strsplit(x = Individual_MiHA_table_files[id], split = "_"))[1])
  if(groupID == 0) stop(id)
  HLA_table_id <- which(reformated_HLA_typing_list$GroupID == groupID)
  groupType <- reformated_HLA_typing_list$GroupType[HLA_table_id]
  num_pos <- dim(MiHA_table)[1]
  
  for(jd in 1:num_pos){
    
    MiHA_db_ind <- which(Restircted_MiHA_db$CHROM_POS %in% MiHA_table$CHROM_POS[jd])
    # if(length(MiHA_db_ind) > 1) stop(jd)
    Restricted_HLA <- unique(Restircted_MiHA_db$HLA[MiHA_db_ind])
    # length(which(reformated_HLA_typing_list[HLA_table_id, 8:17] %in% Restricted_HLA))
    matched_HLA_id <- which(Restricted_HLA %in% as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))
    
    if(length(matched_HLA_id) == 1){ # restricted MiHA Presence
      # if(is.element(Restricted_HLA, as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))){ # restricted MiHA Presence
      
      for(kd in 1:length(matched_HLA_id)){
        
        if(groupType =="a"){
          MiHA_db_ind_id <- which(Restircted_MiHA_db$HLA[MiHA_db_ind] %in% Restricted_HLA[matched_HLA_id[kd]])
          
          aGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] <- aGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] + 1
          firstZeroIndex1 <- which(aGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], -1] == 0)[1] + 1
          aGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], firstZeroIndex1] <- groupID
          
        }else{
          MiHA_db_ind_id <- which(Restircted_MiHA_db$HLA[MiHA_db_ind] %in% Restricted_HLA[matched_HLA_id[kd]])
          
          nGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] <- nGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] + 1
          firstZeroIndex2 <- which(nGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], -1] == 0)[1] + 1
          nGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], firstZeroIndex2] <- groupID
          
        }
        
      }
      
    }
    
  }
  
  
}

total_counts <- cbind(aGVHD_Restricted_MiHA_table_count, nGVHD_Restricted_MiHA_table_count$Count)
colnames(total_counts)[3:4] <- c("aGVHD", "nGVHD")
# write.csv(total_counts, file = "../FirstPaper/Table/HLI_Restricted_MiHA_All_counts_corrected.csv", row.names = F)

total_counts_rm_0s <-total_counts[-which(rowSums(total_counts[,3:4]) == 0), ]
# write.csv(total_counts_rm_0s, file = "../FirstPaper/Table/HLI_Restricted_MiHA_All_counts_rm0_CORRECTED.csv", row.names = F)

#### add MiHA-HLA-SNP-gene
num_MiHas <- dim(total_counts_rm_0s)[1]
comprehensive_MiHA_count_table <- data.frame(MiHA = character(num_MiHas),
                                             Restricted_HLA = character(num_MiHas),
                                             SNP = character(num_MiHas),
                                             Gene = character(num_MiHas),
                                             aGVHD = numeric(num_MiHas),
                                             nGVHD = numeric(num_MiHas),
                                             LLR = numeric(num_MiHas),
                                             stringsAsFactors = F)
for(id in 1:num_MiHas){
  
  index <- which(Restircted_MiHA_db$HLA_SNP %in% total_counts_rm_0s$HLA_SNP[id])
  comprehensive_MiHA_count_table$MiHA[id] <- as.character(Restircted_MiHA_db$KnownMiHA[index])
  comprehensive_MiHA_count_table$Restricted_HLA[id] <- as.character(Restircted_MiHA_db$HLA[index])
  comprehensive_MiHA_count_table$SNP[id] <- as.character(Restircted_MiHA_db$rs.number[index])
  comprehensive_MiHA_count_table$Gene[id] <- as.character(Restircted_MiHA_db$Gene[index])
  comprehensive_MiHA_count_table$aGVHD[id] <- total_counts_rm_0s$aGVHD[id]
  comprehensive_MiHA_count_table$nGVHD[id] <- total_counts_rm_0s$nGVHD[id]
  comprehensive_MiHA_count_table$LLR[id] <- log10(total_counts_rm_0s$aGVHD[id]/total_counts_rm_0s$nGVHD[id])
  
}

#### HLA counts
Restricted_HLA_count <- read.csv(file = "../FirstPaper/Table/Restricted_HLA_summary_all_0630.csv")
colnames(Restricted_HLA_count)[1] <- "HLA"

### Hypergeometric distribution
# dhyper(x, m, n, k, log = FALSE) # density function
# phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE) # distribution function
# qhyper(p, m, n, k, lower.tail = TRUE, log.p = FALSE) # quantile function
# rhyper(nn, m, n, k)  # random generation
##
# x - the number of white balls drawn without replacement from an urn;; the number of aGVHD with RestrictedHLA+ MiHA SNP
# m - the number of white balls in the urn. ;; the number of aGVHD in the HLA type
# n - the number of black balls in the urn. ;; the number of nGVHD in the HLA type
# k - the number of balls drawn from the urn. ;; the number of aGVHD+nGVHD with RestrictedHLA + MiHA SNPt

comprehensive_MiHA_count_table$pValue <- 0
for(id in 1:num_MiHas){
  HLA_ind <- which(Restricted_HLA_count$HLA %in% comprehensive_MiHA_count_table$Restricted_HLA[id])
  comprehensive_MiHA_count_table$pValue[id] <- phyper(comprehensive_MiHA_count_table$aGVHD[id],
                                                      Restricted_HLA_count$aGVHD[HLA_ind],
                                                      Restricted_HLA_count$nGVHD[HLA_ind],
                                                      comprehensive_MiHA_count_table$aGVHD[id]+comprehensive_MiHA_count_table$nGVHD[id])
  
}
comprehensive_MiHA_count_table$fdr_p <- 0
fdr_rate <- 0.05
sorted_rank <- order(comprehensive_MiHA_count_table$pValue, decreasing = F)
comprehensive_MiHA_count_table$fdr_p <- sorted_rank/num_MiHas*fdr_rate
#### Log-Odds ratio plot
comprehensive_MiHA_count_table$MiHA_HLA_SNP_Gene <- sapply(1:dim(comprehensive_MiHA_count_table)[1], 
                                                           function(x) paste0(comprehensive_MiHA_count_table[x, 1:4],
                                                                              collapse = "<>"))
rm_0ids <- which(comprehensive_MiHA_count_table$LLR == 0)
MiHA_LLR <- comprehensive_MiHA_count_table[-rm_0ids, c(10,7)]
mc_LLR <- MiHA_LLR[order(MiHA_LLR$LLR, decreasing = TRUE), ]
mc_LLR <- within(mc_LLR, MiHA_HLA_SNP_Gene <- factor(MiHA_HLA_SNP_Gene, levels=factor(mc_LLR$MiHA_HLA_SNP_Gene)))

library(ggplot2)
ggplot(mc_LLR, aes(x = MiHA_HLA_SNP_Gene, ymax = LLR, ymin = 0, color = "#D55E00")) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  ggtitle("LOD of HLA restricted MiHAs [log10(aGVHD/nGVHD)]") +
  theme_bw() +
  xlab("") +
  # scale_colour(values = "#D55E00") +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none") +
  coord_flip() 


#### Cooccurence table
#################
################
# counts table
################
library(corrplot)

table2matrix <- function(MiHA_count_table){
  #MiHA_count_table <-  aGVHD_Restricted_MiHA_table
  
  num_SNPs <- dim(MiHA_count_table)[1]
  
  groupIDs <- unique(c(as.matrix(MiHA_count_table[, -1])))
  groupIDs <- groupIDs[-which(groupIDs == 0)]
  num_case <- length(groupIDs)
  
  MiHA_SNP_Matrix <- matrix(data = 0, nrow = num_case, ncol = num_SNPs)
  rownames(MiHA_SNP_Matrix) <- as.character(groupIDs)
  colnames(MiHA_SNP_Matrix) <- MiHA_count_table$MiHA_SNP
  for(id in 1:num_SNPs){
    
    col_ind <- which(colnames(MiHA_SNP_Matrix) %in% MiHA_count_table$MiHA_SNP[id])
    tb_gpIDs <- MiHA_count_table[id, -1]
    tb_gpIDs <- as.character(tb_gpIDs[-which(tb_gpIDs == 0)])
    row_ids <- which(rownames(MiHA_SNP_Matrix) %in% tb_gpIDs)
    MiHA_SNP_Matrix[row_ids, col_ind] = MiHA_SNP_Matrix[row_ids, col_ind] + 1
  }
  return(MiHA_SNP_Matrix)
}

aGVHD_SNP_Mat <- table2matrix(aGVHD_Restricted_MiHA_table)
nGVHD_SNP_Mat <- table2matrix(nGVHD_Restricted_MiHA_table)

##########
# correlation_mat <- cor(SNP_mat[, -c(1,2)], use = "pairwise.complete.obs", method = "pearson")
# 
comprehensive_MiHA_count_table$HLA_SNP <- sapply(1:num_MiHas, function(x) paste0(comprehensive_MiHA_count_table$Restricted_HLA[x], 
                                                                                 "-", comprehensive_MiHA_count_table$SNP[x]))

correlation_pairs <- function(SNP_mat){
  
  correlation_coeff <- cor(SNP_mat, method = "spearman")
  num_Snps <- dim(SNP_mat)[2]
  p_value_mat <- matrix(data = 1, nrow = num_Snps, ncol = num_Snps)
  rownames(p_value_mat) <- colnames(SNP_mat)
  colnames(p_value_mat) <- colnames(SNP_mat)
  for(id in 1:(num_Snps-1)){
    
    for(jd in (id+1):num_Snps){
      
      if(correlation_coeff[id, jd] > 0 ){
        alt = "greater"
      }else{alt = "less"}
      
      p_value_mat[id, jd] <- cor.test(SNP_mat[, id], SNP_mat[, jd], alternative = alt, method = "spearman")[["p.value"]]
      p_value_mat[jd, id] <- cor.test(SNP_mat[, id], SNP_mat[, jd], alternative = alt, method = "spearman")[["p.value"]]
      
    }
    
  }
  temp_p_mat <- p_value_mat
  temp_p_mat[lower.tri(temp_p_mat)] <- 1
  vec_ind <- which(temp_p_mat < 0.05)
  
  rid <- floor(vec_ind/num_Snps) 
  row_ind <- vec_ind - rid * num_Snps  
  col_ind <- rid + 1
  
  sig_pairs <- sapply(1:length(row_ind), function(x) paste(p_value_mat[row_ind[x], col_ind[x]], 
                                                           rownames(p_value_mat)[row_ind[x]], 
                                                           colnames(p_value_mat)[col_ind[x]], 
                                                           correlation_coeff[row_ind[x], col_ind[x]]))
  
  return(list(cor_coeff_mat = correlation_coeff, p_values = p_value_mat, sigPairs = sig_pairs))
}

checkSigGenePairs <- function(sigPairs, MiHA_convert_table){
  
  num_pairs <- length(sigPairs)
  GenePairs <- data.frame(MiHA_HLA_SNP_Gene1 = character(num_pairs),
                          MiHA_HLA_SNP_Gene2 = character(num_pairs),
                          correlation_coeff = numeric(num_pairs),
                          p_values = numeric(num_pairs),
                          stringsAsFactors = F) 
  for(id in 1:num_pairs){
    
    PairInfo <- unlist(strsplit(sigPairs[id], " "))
    GenePairs$correlation_coeff[id] <- round(as.numeric(PairInfo[4]), digits = 3)
    GenePairs$p_values[id] <- round(as.numeric(PairInfo[1]), digits = 3)
    
    GenePairs$MiHA_HLA_SNP_Gene1[id] <- MiHA_convert_table$MiHA_HLA_SNP_Gene[which(MiHA_convert_table$HLA_SNP %in% PairInfo[2])]
    GenePairs$MiHA_HLA_SNP_Gene2[id] <- MiHA_convert_table$MiHA_HLA_SNP_Gene[which(MiHA_convert_table$HLA_SNP %in% PairInfo[3])]
    
    GenePairs
  }
  
  return(GenePairs)
}

## A*02:01 restricted MiHAs
aSNP_mat_A0201 <- aGVHD_SNP_Mat[, which(grepl("A\\*02:01", colnames(aGVHD_SNP_Mat)))]
aSNP_mat_A0201 <- aSNP_mat_A0201[-which(rowSums(aSNP_mat_A0201) == 0), ]
aSNP_mat_A0201 <- aSNP_mat_A0201[, -which(colSums(aSNP_mat_A0201) == 0)]
corrgram(aSNP_mat_A0201, order = F, lower.panel = panel.shade,
         upper.panel = panel.pie, cor.method = "spearman")
colSums(aSNP_mat_A0201)

nSNP_mat_A0201 <- nGVHD_SNP_Mat[, which(grepl("A\\*02:01", colnames(nGVHD_SNP_Mat)))]
nSNP_mat_A0201 <- nSNP_mat_A0201[-which(rowSums(nSNP_mat_A0201) == 0), ]
nSNP_mat_A0201 <- nSNP_mat_A0201[, -which(colSums(nSNP_mat_A0201) == 0)]
corrgram(nSNP_mat_A0201, order = F, lower.panel = panel.shade,
         upper.panel = panel.pie, cor.method = "spearman")
colSums(nSNP_mat_A0201)

aGVHD_correlation_list <- correlation_pairs(aSNP_mat_A0201)
nGVHD_correlation_list <- correlation_pairs(nSNP_mat_A0201)

checkSigGenePairs(aGVHD_correlation_list[["sigPairs"]], comprehensive_MiHA_count_table[, 10:11])
checkSigGenePairs(nGVHD_correlation_list[["sigPairs"]], comprehensive_MiHA_count_table[, 10:11])

corrplot(aGVHD_correlation_list[["cor_coeff_mat"]], type = "upper", method = "circle")
corrplot(nGVHD_correlation_list[["cor_coeff_mat"]], type = "lower", method = "circle")

corrplot(aGVHD_correlation_list[["cor_coeff_mat"]], type = "upper", method = "ellipse", 
         p.mat = aGVHD_correlation_list[["p_values"]], sig.level=0.05, insig = "pch", pch.col = "red", pch.cex = 0.5)
corrplot(nGVHD_correlation_list[["cor_coeff_mat"]], type = "upper", method = "ellipse", 
         p.mat = nGVHD_correlation_list[["p_values"]], sig.level=0.05, insig = "pch", pch.col = "red", pch.cex = 0.5)

## B*07:02
aSNP_mat_B0702 <- aGVHD_SNP_Mat[, which(grepl("B\\*07:02", colnames(aGVHD_SNP_Mat)))]
aSNP_mat_B0702 <- aSNP_mat_B0702[-which(rowSums(aSNP_mat_B0702) == 0), ]
aSNP_mat_B0702 <- aSNP_mat_B0702[, -which(colSums(aSNP_mat_B0702) == 0)]

nSNP_mat_B0702 <- nGVHD_SNP_Mat[, which(grepl("B\\*07:02", colnames(nGVHD_SNP_Mat)))]
nSNP_mat_B0702 <- nSNP_mat_B0702[-which(rowSums(nSNP_mat_B0702) == 0), ]
nSNP_mat_B0702 <- nSNP_mat_B0702[, -which(colSums(nSNP_mat_B0702) == 0)]

aGVHD_correlation_list <- correlation_pairs(aSNP_mat_B0702)
nGVHD_correlation_list <- correlation_pairs(nSNP_mat_B0702)

checkSigGenePairs(aGVHD_correlation_list[["sigPairs"]], comprehensive_MiHA_count_table[, 10:11])
checkSigGenePairs(nGVHD_correlation_list[["sigPairs"]], comprehensive_MiHA_count_table[, 10:11])


corrplot(aGVHD_correlation_list[["cor_coeff_mat"]], type = "upper", method = "circle")
corrplot(nGVHD_correlation_list[["cor_coeff_mat"]], type = "lower", method = "circle")

corrplot(aGVHD_correlation_list[["cor_coeff_mat"]], type = "lower", method = "ellipse", 
         p.mat = aGVHD_correlation_list[["p_values"]], sig.level=0.05, insig = "pch", pch.col = "red", pch.cex = 0.5)
corrplot(nGVHD_correlation_list[["cor_coeff_mat"]], type = "lower", method = "ellipse", 
         p.mat = nGVHD_correlation_list[["p_values"]], sig.level=0.05, insig = "pch", pch.col = "red", pch.cex = 0.5)



##########
#### RandomForest variable importance
library(randomForest)

# Feat <- rbind(aGVHD_SNP_Mat, nGVHD_SNP_Mat)
Feat <- rbind(aGVHD_SNP_Mat, nGVHD_SNP_Mat)
Labels <- c(rep("a", dim(aGVHD_SNP_Mat)[1]), 
            rep("n", dim(nGVHD_SNP_Mat)[1]))
Labels <- as.factor(Labels)
feat_lab <- cbind(Feat, Labels)

RF_feat <- randomForest(x = Feat, y = Labels, importance=TRUE)
# round(importance(RF_feat), 2)
imp <- importance(RF_feat, type =1, scale=T) 
importances_order <- order(imp, decreasing = T)
ordered_feat_set <- colnames(Feat)[importances_order]
ordered_feat_imp <- imp[importances_order]

feat_imp_table <- data.frame(MiHA_HLA_SNP_Gene = character(length(ordered_feat_set)),
                             HLA_SNP = ordered_feat_set,
                             FeatScore = ordered_feat_imp,
                             stringsAsFactors = F)
feat_imp_table$MiHA_HLA_SNP_Gene <- sapply(1:length(ordered_feat_set), function(x) comprehensive_MiHA_count_table$MiHA_HLA_SNP_Gene[which(comprehensive_MiHA_count_table$HLA_SNP %in% ordered_feat_set[x])])


################################### GWAS ### 
#################
load("../Data/GWAS_HLI_overlapping_SNPs_wHLA.RData") ## overlapping_SNPs
overlapping_SNPs$CHROM_POS <- overlapping_SNPs$chr_POS
load("../Data/ID_table_wCaseID.RData")          ## ID_table

GWAS_avail_cases <- read.csv(file = "../ClinVar/GWASH/GWASH_available_IDs.csv")
############################
## MiHA Combination
AML_cases_table <- GWAS_avail_cases

Num_AML_cases <- dim(AML_cases_table)[1]
mismatch_MiHA_table <- matrix(data = 0, nrow = Num_AML_cases, ncol = 2*num_SNPs+2)
colnames(mismatch_MiHA_table) <- c("CaseNumber", "aGVHD", 
                                   as.vector(sapply(1:num_SNPs, function(x) paste0(MiHA_SNPs$SNP_ID[x], c("_PresAbs", "_Number"))))) 
## presense-absence of the MihA mimatch + number of mismatched alleles


mismatch_MiHA_table <- as.data.frame(mismatch_MiHA_table)
counter <- 0
for(id in 1:Num_AML_cases){
  
  index <- which(MiHA_genotypes$case_id %in% AML_cases_table$bmt_case[id])
  
  if(length(index) == 2){ # matching cases
    counter <- counter + 1
    mismatch_MiHA_table$CaseNumber[counter] <- AML_cases_table$bmt_case[id]
    mismatch_MiHA_table$aGVHD[counter] <- AML_cases_table$agvhi24[id]
    
    for(jd in 1:num_SNPs){
      
      mismatch_MiHA_table[counter, 4+2*(jd-1)] <- identical(MiHA_genotypes[index[1], SNP_index[1+2*(jd-1)]], MiHA_genotypes[index[2], SNP_index[1+2*(jd-1)]]) + 
        identical(MiHA_genotypes[index[1], SNP_index[2*jd]], MiHA_genotypes[index[2], SNP_index[2*jd]])
      
      mismatch_MiHA_table[counter, 3+2*(jd-1)] <- as.integer(identical(MiHA_genotypes[index[1], SNP_index[1+2*(jd-1)]], MiHA_genotypes[index[2], SNP_index[1+2*(jd-1)]]) | 
                                                               identical(MiHA_genotypes[index[1], SNP_index[2*jd]], MiHA_genotypes[index[2], SNP_index[2*jd]]))
      
    }
    
    
  }
  
  
}

mismatch_MiHA_table <- mismatch_MiHA_table[1:counter, ]

all_possible_pairs_num <- num_SNPs * (num_SNPs -1) / 2
# all_possible_pairs_names <- character(all_possible_pairs_num)
a_pair_wise_combination <- as.data.frame(matrix(data = 0, nrow = all_possible_pairs_num, ncol = 5), stringsAsFactors = F) # col: abs - abs 0-0;  pres-abs 1-0;  abs-pres 0-1; pres - pres 1-1
colnames(a_pair_wise_combination) <- c("SNPComb", "abs-abs", "pres-abs", "abs-pres", "pres-pres")
n_pair_wise_combination <- as.data.frame(matrix(data = 0, nrow = all_possible_pairs_num, ncol = 5), stringsAsFactors = F)
colnames(n_pair_wise_combination) <- c("SNPComb", "abs-abs", "pres-abs", "abs-pres", "pres-pres")

aGVHD_index <- which(mismatch_MiHA_table$aGVHD == 1) # 146
nGVHD_index <- which(mismatch_MiHA_table$aGVHD == 0) # 186
counter <- 0 
for(id in 1:(num_SNPs-1)){
  
  for(jd in (id+1):num_SNPs){
    # if(id != jd){
    
    counter <- counter + 1
    
    sum_tab <- as.data.frame(table(mismatch_MiHA_table[aGVHD_index, c(3 + 2*(id-1), 3 + 2*(jd-1))]), stringsAsFactors = F)
    
    SNP_1 <- gsub("_PresAbs", "", colnames(sum_tab)[1])
    SNP_2 <- gsub("_PresAbs", "", colnames(sum_tab)[2])
    
    a_pair_wise_combination[counter, 1] <- paste0(SNP_1, "-", SNP_2)
    a_pair_wise_combination[counter, 2:5] <- as.numeric(sum_tab$Freq)
    
    if(!identical(as.numeric(sum_tab[2, c(1,2)]),  c(1, 0))) cat("aGVHD: id = ", id, "; jd = ", jd, "\n")
    
    nsum_tab <- as.data.frame(table(mismatch_MiHA_table[nGVHD_index, c(3 + 2*(id-1), 3 + 2*(jd-1))]), stringsAsFactors = F)
    
    nSNP_1 <- gsub("_PresAbs", "", colnames(nsum_tab)[1])
    nSNP_2 <- gsub("_PresAbs", "", colnames(nsum_tab)[2])
    
    n_pair_wise_combination[counter, 1] <- paste0(nSNP_1, "-", nSNP_2)
    n_pair_wise_combination[counter, 2:5] <- as.numeric(nsum_tab$Freq)
    
    if(!identical(as.numeric(nsum_tab[2, c(1,2)]),  c(1, 0))) cat("nGVHD: id = ", id, "; jd = ", jd, "\n")
    
    
  }
  
  # }
  
}

save(mismatch_MiHA_table, file = "../Data/GWASH_mismatch_MiHA_table_AllSamples.RData")
##### 45 MiHA SNPs --- Outcome vs rs-rs association test
a_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_summary_AP.csv", stringsAsFactors = F)
n_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_summary_AP.csv", stringsAsFactors = F)

a_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_ordered_AP.csv", stringsAsFactors = F)
n_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_ordered_AP.csv", stringsAsFactors = F)

a_overlapped_combination <- intersect(a_joint_SNPs$X, a_pair_wise_combination$SNPComb) # 102 
n_overlapped_combination <- intersect(n_joint_SNPs$X, n_pair_wise_combination$SNPComb) # 102

a_ordered_name <- sapply(1:length(a_joint_SNPs$X), function(x) paste0(sort(unlist(strsplit(a_joint_SNPs$X[x], "-"))), collapse = "-"))
n_ordered_name <- sapply(1:length(n_joint_SNPs$X), function(x) paste0(sort(unlist(strsplit(n_joint_SNPs$X[x], "-"))), collapse = "-"))

aaa_ordered_name <- sapply(1:length(a_pair_wise_combination$SNPComb), function(x) paste0(sort(unlist(strsplit(a_pair_wise_combination$SNPComb[x], "-"))), collapse = "-"))
nnn_ordered_name <- sapply(1:length(n_pair_wise_combination$SNPComb), function(x) paste0(sort(unlist(strsplit(n_pair_wise_combination$SNPComb[x], "-"))), collapse = "-"))

a_overlapped_combination_ordered <- intersect(a_ordered_name, aaa_ordered_name) # 190
n_overlapped_combination_ordered <- intersect(n_ordered_name, nnn_ordered_name) # 190

## verification of overlapping MiHA combination between the HLI cohort and GWAS cohort
aind <- which(!(aaa_ordered_name %in% a_overlapped_combination_ordered))
aind2 <- which(!(a_ordered_name %in% a_overlapped_combination_ordered))

a_num_pairs <- sapply(1:length(a_joint_SNPs$X), function(x) length(which((unlist(strsplit(a_joint_SNPs$X[x], "-"))) %in% MiHA_SNPs$SNP_ID))) # length(which(a_num_pairs == 2)): 190
n_num_pairs <- sapply(1:length(n_joint_SNPs$X), function(x) length(which((unlist(strsplit(n_joint_SNPs$X[x], "-"))) %in% MiHA_SNPs$SNP_ID))) # length(which(n_num_pairs == 2)): 190

########## verification ends here ####################

####### combine two tables

## aGVHD
a_GWAS_ind <- which((aaa_ordered_name %in% a_overlapped_combination_ordered))
a_HLI_ind <- which((a_ordered_name %in% a_overlapped_combination_ordered))

a_combined_table <- a_pair_wise_combination[a_GWAS_ind, ]
# abs-abs; pres-abs; abs-pres; pres-pres
num_overlapped <- length(a_GWAS_ind)
for(id in 1:num_overlapped){
  
  index <- which(a_ordered_name %in% aaa_ordered_name[a_GWAS_ind[id]])
  
  a_combined_table$`abs-abs`[id] <- a_combined_table$`abs-abs`[id] + a_joint_SNPs$absence.absence[index]
  a_combined_table$`pres-abs`[id] <- a_combined_table$`pres-abs`[id] + a_joint_snps_ordered$presence.HET..absence[index] + a_joint_snps_ordered$presence.HOM..absence[index]
  a_combined_table$`abs-pres`[id] <- a_combined_table$`abs-pres`[id] + a_joint_snps_ordered$absence.presence.HET.[index] + a_joint_snps_ordered$absence.presence.HOM.[index]
  a_combined_table$`pres-pres`[id] <- a_combined_table$`pres-pres`[id] + a_joint_snps_ordered$presence.HET..presence.HOM.[index] + a_joint_snps_ordered$presence.HOM..presence.HET.[index]
  
}
a_combined_table <- rbind(a_combined_table, a_pair_wise_combination[-a_GWAS_ind, ])

# nGVHD
n_GWAS_ind <- which((nnn_ordered_name %in% n_overlapped_combination_ordered))
n_HLI_ind <- which((n_ordered_name %in% n_overlapped_combination_ordered))

n_combined_table <- n_pair_wise_combination[n_GWAS_ind, ]
# abs-abs; pres-abs; abs-pres; pres-pres
num_overlapped <- length(n_GWAS_ind)
for(id in 1:num_overlapped){
  
  index <- which(n_ordered_name %in% nnn_ordered_name[n_GWAS_ind[id]])
  
  n_combined_table$`abs-abs`[id] <- n_combined_table$`abs-abs`[id] + n_joint_SNPs$absence.absence[index]
  n_combined_table$`pres-abs`[id] <- n_combined_table$`pres-abs`[id] + n_joint_snps_ordered$presence.HET..absence[index] + n_joint_snps_ordered$presence.HOM..absence[index]
  n_combined_table$`abs-pres`[id] <- n_combined_table$`abs-pres`[id] + n_joint_snps_ordered$absence.presence.HET.[index] + n_joint_snps_ordered$absence.presence.HOM.[index]
  n_combined_table$`pres-pres`[id] <- n_combined_table$`pres-pres`[id] + n_joint_snps_ordered$presence.HET..presence.HOM.[index] + n_joint_snps_ordered$presence.HOM..presence.HET.[index]
  
}
n_combined_table <- rbind(n_combined_table, n_pair_wise_combination[-n_GWAS_ind, ])
