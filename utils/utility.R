
filter_outcome <- function(clinical.fn, genotype){
  
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
  
  cat(nrow(Outcome)-length(case_idx), "cases has been filtered out.\n")
  
  return(Outcome_avail)
}