
paired_vcf_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"

all_vcf_files <- list.files(paired_vcf_dir, pattern = "\\.vcf.gz$")

output_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/vcf_chr6/"

temp_output_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/temp_a_vcf_chr6/"

existed_files <- list.files(temp_output_dir, pattern = "\\.vcf.gz$")
groupIDs <- gsub("_chr6.vcf.gz", "", existed_files)

all_groupIDs <- gsub(".vcf.gz", "", all_vcf_files)
remain_ids <- which(!(all_groupIDs %in% groupIDs))
all_groupIDs <- all_groupIDs[remain_ids]
all_vcf_files <- all_vcf_files[remain_ids]

num_files <- length(all_vcf_files)

count <- 0

system(paste0("mkdir ", temp_output_dir))
for(id in 1:num_files){
  
  if(grepl("a_", all_vcf_files[id])){
    count <- count + 1
    
    # if(count <= 30){
    
    ptm <- proc.time()
    
    cat("Paired file #", count, "\n")
    out.filename <- gsub(".vcf.gz", "_chr6.vcf.gz", all_vcf_files[id])
    system(paste0("tabix -h ", paired_vcf_dir, all_vcf_files[id], " chr6 | bgzip > ", temp_output_dir, out.filename))
    system(paste0("cd ", temp_output_dir, "; tabix -p vcf ", out.filename))
    cat(paste0("tabix -h ", paired_vcf_dir, all_vcf_files[id], " chr6 | bgzip > ", temp_output_dir, out.filename), "\n")
    cat(paste0("cd ", temp_output_dir, "; tabix -p vcf ", out.filename), "\n")
    print(proc.time()-ptm)
  } else break
  # }
}
cat("Extraction done! \n")
all_vcf_files <- list.files(temp_output_dir, pattern = "\\.vcf.gz$")
cat("All files #", length(all_vcf_files), "\n")
file_list <- paste(all_vcf_files, collapse = " ")

ptm <- proc.time()
system(paste0("cd ", temp_output_dir, "; vcf-merge -R 0/0 ", file_list, " | bgzip -c > ", output_dir, "aGVHD_all_pairs_chr6.vcf.gz"))
system(paste0("cd ", output_dir, "; tabix -p vcf aGVHD_all_pairs_chr6.vcf.gz"))
#system(paste0("rm -r ", temp_output_dir))
cat(paste0("cd ", temp_output_dir, "; vcf-merge -R 0/0 ", file_list, " | bgzip -c > ", output_dir, "aGVHD_all_pairs_chr6.vcf.gz"), "\n")
cat(paste0("cd ", output_dir, "; tabix -p vcf aGVHD_all_pairs_chr6.vcf.gz"), "\n")
#cat(paste0("rm -r ", temp_output_dir),"\n")

print(proc.time() - ptm)
