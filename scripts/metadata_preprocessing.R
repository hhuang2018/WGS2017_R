## install packages as needed - only need once
## utils/InstallPacakges.R
#####

input.dir <- "../Data/"
output.dir <- "../R_Output/"

# input
gwas.fn <- lapply(c(bed = 'bed', bim='bim', fam='fam', map='map'), function(n) sprintf("%s/22bed.%s", input.dir, n))
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
Outcome <- read.csv(clinical.fn)

individuals <- rownames(genotype)

num_individual <- length(individuals)
keep_ridx <- vector(mode = "numeric", length = 0)
keep_didx <- keep_ridx
for(id in 1:num_individual){
  
  if(nchar(individuals[id]) == 9){ # rid
    
    idx <- which(Outcome$rid %in% as.numeric(gsub("-", "", individuals[id])))
    
    if(length(idx) == 1){
      
      keep_ridx <- c(keep_ridx, idx)
      
    }else if(length(idx) > 1){
      
      cat(paste0("Duplicated Recipient ", individuals[id], "\n"))
      
      keep_ridx <- c(keep_ridx, idx) 
    }
    
  }else if (nchar(individuals[id]) == 11){ # did
    
    idx <- which(Outcome$did %in% as.numeric(gsub("-", "", individuals[id])))
    
    if(length(idx) == 1){
      
      keep_didx <- c(keep_didx, idx)
      
    }else if(length(idx) > 1){
      
      cat(paste0("Duplicated Donor ", individuals[id], "\n"))
      keep_didx <- c(keep_didx, idx)
      
    }
    
  }else {
    
    cat(paste0(individuals[id], " is an invalid ID. \n"))
    
  }
  
}

# length(keep_ridx) # 988 recipients avaliable
# length(keep_didx) # 996 donors available ### Duplicated Donor 0255-7626-5
# length(intersect(keep_ridx, keep_didx)) # 988 paired cases available

case_idx <- intersect(keep_ridx, keep_didx)

Outcome_avail <- Outcome[case_idx, ]
write.csv(Outcome_avail, file = "../R_Output/available_cases.csv", row.names = F)

##### 
table(Outcome_avail$agvhi24)
# 0   1 
# 531 456 
table(Outcome_avail$agvhi34)
# 0   1 
# 754 233
table(Outcome_avail$cgvhi)
# 0   1 
# 478 509
table(Outcome_avail$disease)
# 10  20  40  50 
# 332 123 343 190 
table(Outcome_avail$cgvhi[Outcome_avail$disease == 40])

table(Outcome_avail$sexmatch)
# 1   2   3   4 
# 416 280 129 163 
# 1: Male -> Male
# 2: Male -> Female
# 3: Female -> Male
# 4: Female -> Female

###### 
table(Outcome_avail$agvhi34[Outcome_avail$sexmatch == 2]) 

table(Outcome_avail$disease[Outcome_avail$sexmatch == 2])


table(Outcome_avail$cgvhi[Outcome_avail$sexmatch == 2 & Outcome_avail$disease == 40])

###############################
hla_tab <- read.csv("../Data/GWAS_HLA.csv")

HLA_tab <- hla_tab[hla_tab$bmt_case_num %in% intersect(hla_tab$bmt_case_num, Outcome_avail$bmt_case), ]
write.csv(HLA_tab, file = "../R_Output/available_cases_HLA.csv", row.names = F)
