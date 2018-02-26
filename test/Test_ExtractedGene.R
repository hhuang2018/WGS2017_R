#require(snpStats)

fileFp <- "../../HLI/GeneData/"
outputFp <- "../../HLI/GenotypeTables/"

GeneName <- 'IL10'
files <- list.files(fileFp, pattern = paste0("\\",GeneName, ".vcf$"))

load("../../HLI/Metadata/HLI_available_pairs_dis_table.RData")

pairNum <- dim(Available_paired_table)[1]

GeneGenoType_List <- data.frame(CHROM = character(0),
                                POS = numeric(0),
                                REF = character(0),
                                ALT = character(0),
                                ANN = character(0),
                                CHRPOS = character(0),
                                GT = character(0))
for(id in 1:pairNum){
  
  idx <- grep(paste0("_", Available_paired_table$GroupID[id],"_"), files)
  paired_GenoType <- data.frame(CHROM = character(0),
                                POS = numeric(0),
                                REF = character(0),
                                ALT = character(0),
                                GT = character(0),
                                ANN = character(0),
                                CHRPOS = character(0))
  for(jd in 1:length(idx)){
    Record <- read.delim(paste0(fileFp, files[idx[jd]]), sep = "\t", header = F, stringsAsFactors = F)
    colnames(Record) <- c("CHROM", "POS", "RS", "REF", "ALT", "QC", "FT", "ANN", "FORMAT", "GT")
    Record$ANN2 <- sapply(1:dim(Record)[1], function(x) strsplit(Record$ANN[x], ";")[[1]][3])
    Record$GT2 <- sapply(1:dim(Record)[1], function(x) paste0(sapply(strsplit(strsplit(Record$GT[x], ":")[[1]][1],"/")[[1]], function(y) if(y=="0") Record$REF[x] else Record$ALT[x]), collapse = "/"))
    type <- strsplit(files[idx[jd]], "_")[[1]][3]
    
    Record$CHRPOS <- sapply(1:dim(Record)[1], function(x) paste0(Record$CHROM[x], "-",Record$POS[x]))
    if(length(paired_GenoType$CHROM) >0){
      paired_GenoType$CHRPOS <- sapply(1:dim(paired_GenoType)[1], function(x) paste0(paired_GenoType$CHROM[x], "-",paired_GenoType$POS[x]))
      for(kd in 1:dim(Record)[1]){
        kdx <- which(paired_GenoType$CHRPOS %in% Record$CHRPOS[kd])
        if(length(kdx) == 1){ 
          ## Donor_GenoType - Recipient_GenoType
          ifelse(type == "D", paired_GenoType$GT[kdx] <- paste0(Record$GT2[kd], paired_GenoType$GT[kdx]), 
                 paired_GenoType$GT[kdx] <- paste0(paired_GenoType$GT[kdx], Record$GT2[kd]))
        }else{
          temp_GenoType <- data.frame(CHROM = character(1),
                                        POS = numeric(1),
                                        REF = character(1),
                                        ALT = character(1),
                                        GT = character(1),
                                        ANN = character(1),
                                      CHRPOS = character(1))
          temp_GenoType$CHROM <- Record$CHROM[kd]
          temp_GenoType$POS <- Record$POS[kd]
          temp_GenoType$REF <- Record$REF[kd]
          temp_GenoType$ALT <- Record$ALT[kd]
          temp_GenoType$GT <- ifelse(type == "D", paste0(Record$GT2[kd], "-"), paste0("-", Record$GT2[kd]))
          temp_GenoType$ANN <- Record$ANN2[kd]
          temp_GenoType$CHRPOS <- Record$CHRPOS[kd]
          paired_GenoType <- rbind(paired_GenoType, temp_GenoType)
          
        }
      }
    }else{
      
      paired_GenoType <-  Record[, c("CHROM", "POS", "REF", "ALT", "GT2", "ANN2")]
      colnames(paired_GenoType) <- c("CHROM", "POS", "REF", "ALT", "GT", "ANN")
      if(type == 'D'){
        paired_GenoType$GT <- paste0(paired_GenoType$GT, '-')
      }else{
        paired_GenoType$GT <- paste0('-',paired_GenoType$GT)
      }
      
    }
  }
  
  colName <- paste0(Available_paired_table$GroupType[id], "_", Available_paired_table$GroupID[id])
  if(length(GeneGenoType_List$CHROM) >0){
    List_CHRPOS <- sapply(1:dim(GeneGenoType_List)[1], function(x) paste0(GeneGenoType_List$CHROM[x], "-",GeneGenoType_List$POS[x]))
    
    eval(parse(text=paste0("GeneGenoType_List$", colName, " <- ''")))
    
    for(kd in 1:dim(paired_GenoType)[1]){
      kdx <- which(List_CHRPOS %in% paired_GenoType$CHRPOS[kd])
      if(length(kdx) == 1){ 
        
        eval(parse(text=paste0("GeneGenoType_List$", colName, "[", as.character(kdx),"] <- paired_GenoType$GT[kd]")))
        
      }else{
        temp_GenoType <- data.frame(CHROM = character(1),
                                    POS = numeric(1),
                                    REF = character(1),
                                    ALT = character(1),
                                    ANN = character(1),
                                    CHRPOS = character(1))
        temp_GenoType$CHROM <- paired_GenoType$CHROM[kd]
        temp_GenoType$POS <- paired_GenoType$POS[kd]
        temp_GenoType$REF <- paired_GenoType$REF[kd]
        temp_GenoType$ALT <- paired_GenoType$ALT[kd]
        temp_GenoType$GT <- paired_GenoType$GT[kd]
        temp_GenoType$ANN <- paired_GenoType$ANN[kd]
        temp_GenoType$CHRPOS <- paired_GenoType$CHRPOS[kd]
        #rbind(temp_GenoType, GeneGenoType_List, stringsAsFactors=FALSE)
        GeneGenoType_List[nrow(GeneGenoType_List)+1, c("CHROM", "POS", "REF", "ALT", "ANN")] <- temp_GenoType[c("CHROM", "POS", "REF", "ALT", "ANN")]
        
        eval(parse(text=paste0("GeneGenoType_List$", colName, "[", as.character(nrow(GeneGenoType_List)),"] <- temp_GenoType$GT")))
        
      }
    }
    
    
  }else{
    
    GeneGenoType_List <- paired_GenoType[, c("CHROM", "POS", "REF", "ALT", "ANN", "GT")]
    colnames(GeneGenoType_List)[dim(GeneGenoType_List)[2]] <- colName
  }
  
}
save(GeneGenoType_List, file = paste0(outputFp, GeneName, "_DRpair_genotypes.RData"))

write.csv(GeneGenoType_List, file = paste0(outputFp, GeneName, "_DRpair_genotypes.csv"))
write.table(GeneGenoType_List, file = paste0(outputFp, GeneName, "_DRpair_genotypes.txt"), sep = "\t", row.names = F)
