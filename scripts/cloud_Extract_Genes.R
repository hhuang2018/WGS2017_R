### Extract Specific Gene reads
# list .vcf.gz files

original_file_dir = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq/"
destination_dir = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_Extracted_genes/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
filenames <- gsub(".vcf.gz", "", filenames)

GeneList <- data.frame(GeneNames = c('IL10', 'IL10RB'),
                       GenePOS = c("chr1:206,767,602-206,772,494", "chr21:33,266,358-33,297,234"))

load("../Data/HLI_available_pairs_dis_table.RData")

for(id in 1:dim(GeneList)[1]){
  extracted_filenames <- list.files(destination_dir, pattern = "\\.vcf.gz$")
  extracted_filenames <- gsub(paste0(GeneList$GeneNames[id], ".vcf.gz"), "", extracted_filenames) 
  
  Not_extracted_index <- which(!(filenames %in% extracted_filenames))
  
  num_files <- length(Not_extracted_index)
  
  for(jd in 1:num_files){
    system(paste0("tabix ", original_file_dir, filenames[Not_extracted_index[jd]],".vcf.gz ", GeneList$GenePOS[id], 
                  " > ", destination_dir, filenames[Not_extracted_index[jd]], "_", GeneList$GeneNames[id],
                  ".vcf"),
           intern = TRUE)
  print(paste0("Extracting File: ", filenames[Not_extracted_index[jd]]))
  }
  
}



