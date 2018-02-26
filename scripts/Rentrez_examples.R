## Rentrez tutorial:
## https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
##
## Entrez Programming Utilities Help:
## https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/
##
##

library(rentrez)
# available DBs
entrez_dbs()

# Brief description of what the database is
entrez_db_summary("snp")

# Set of search terms that can used with this database
entrez_db_searchable("snp")

## Search db: at least need db and term
entrez_search(db   = "pubmed",
              term = "(vivax malaria[MeSH]) AND (folic acid antagonists[MeSH])")

### Adavance counting: count publication
search_year <- function(year, term){
  query <- paste(term, "AND (", year, "[PDAT])")
  entrez_search(db="pubmed", term=query, retmax=0)$count
}

year <- 2008:2017
papers <- sapply(year, search_year, term="CAR-T", USE.NAMES=FALSE)

plot(year, papers, type='b', main="The Rise of the CAR-T")


### Finding cross-references : entrez_link()
all_the_links <- entrez_link(dbfrom='gene', id=351, db='all')
all_the_links
all_the_links$links  # format: [source_database]_[linked_database]

all_the_links$links$gene_pmc[1:10]

all_the_links$links$gene_clinvar
all_cv <- entrez_summary(db="clinvar", id=all_the_links$links$gene_clinvar)
all_cv
extract_from_esummary(all_cv, "gene_sort")

### 
# clinvar example
res <- entrez_search(db = "clinvar", term = "leukemia", retmax=10)
cv <- entrez_summary(db="clinvar", id=res$ids)
cv
extract_from_esummary(cv, "title", simplify=FALSE)
extract_from_esummary(cv, "trait_set")[1:2]
extract_from_esummary(cv, "gene_sort")

## Fetching full records: entrez_fetch()
# Fetch DNA sequences in fasta format
gene_ids <- c(351, 11647)
linked_seq_ids <- entrez_link(dbfrom="gene", id=gene_ids, db="nuccore")
linked_transripts <- linked_seq_ids$links$gene_nuccore_refseqrna
head(linked_transripts)
all_recs <- entrez_fetch(db="nuccore", id=linked_transripts, rettype="fasta")
class(all_recs)
nchar(all_recs)
cat(strwrap(substr(all_recs, 1, 500)), sep="\n")

write(all_recs, file="my_transcripts.fasta")
# Alternatively writing:
temp <- tempfile()
write(all_recs, temp)
parsed_recs <- ape::read.dna(all_recs, temp)

## Fetch a parsed XML document
Tt <- entrez_search(db="taxonomy", term="(Tetrahymena thermophila[ORGN]) AND Species[RANK]")
tax_rec <- entrez_fetch(db="taxonomy", id=Tt$ids, rettype="xml", parsed=TRUE)
class(tax_rec)
tax_list <- XML::xmlToList(tax_rec)
tax_list$Taxon$GeneticCode

tt_lineage <- tax_rec["//LineageEx/Taxon/ScientificName"]
tt_lineage[1:4]

XML::xpathSApply(tax_rec, "//LineageEx/Taxon/ScientificName", XML::xmlValue)


#############
# IL10 and IL10RB snps
for(ind in seq(from=1, to=10000, by=50 )){
  SnpList <- entrez_search(db   = "snp",
                           term = "(IL10[GENE]) AND (human[ORGN])", retstart = ind, retmax = 50)
  #length(SnpList$ids)
  ids <- SnpList$ids
  pop_summ <- entrez_summary(db="snp", id=ids)
  
  MAF = extract_from_esummary(pop_summ, "global_maf")
  if(length(which(MAF != ""))){
    idx = which(MAF != "")
    print(MAF[idx])
  }
}

# 
# esummary result with 46 items:
# [1] uid                   snp_id                organism              allele_origin        
# [5] global_maf            global_population     global_samplesize     suspected            
# [9] clinical_significance genes                 acc                   chr                  
# [13] weight                handle                fxn_class             validated            
# [17] gtype                 nonref                docsum                het                  
# [21] srate                 tax_id                chrrpt                orig_build           
# [25] upd_build             createdate            updatedate            pop_class            
# [29] method_class          snp3d                 linkout               ss                   
# [33] locsnpid              allele                snp_class             chrpos               
# [37] contigpos             text                  lookup                sort_priority        
# [41] snp_id_sort           clinical_sort         human_sort            cited_sort           
# [45] weight_sort           chrpos_sort   

# > pop_summ$`3021093`$chr
# [1] "1"
# > pop_summ$`3021093`$snp3d
# [1] ""
# > pop_summ$`3021093`$snp_class
# [1] "snp"
# > pop_summ$`3021093`$snp_id_sort
# [1] "0003021093"
# > pop_summ$`3021093`$snp_id
# [1] 3021093
# > pop_summ$`3021093`$genes
# name gene_id
# 1 IL10    3586
# > pop_summ$`1800896`$global_maf
# [1] "C=0.2722/1363"
# > pop_summ$`1800896`$suspected
# [1] ""
# > pop_summ$`1800896`$ss
# [1] 276190046
# > pop_summ$`1800896`$lookup
# [1] "2216332910"
# > pop_summ$`1800896`$validated
# [1] "by-1000G,by-2hit-2allele,by-cluster,by-frequency,by-hapmap"
# > pop_summ$`1800896`$chrpos
# [1] "1:206773552"
# > pop_summ$`1800896`$docsum
# [1] "HGVS=NC_000001.10:g.206946897T&gt;C,NC_000001.11:g.206773552T&gt;C,NG_012088.1:g.3943A&gt;G,NM_000572.2:c.-1117A&gt;G|SEQ=CAACACTACTAAGGCTTCTTTGGGA[A/G]GGGGAAGTAGGGATAGGTAAGAGGA|GENE=IL10:3586"