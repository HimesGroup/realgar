library(DBI, quietly = T)
library(RSQLite, quietly=T)

#SQLite for omics database ---
ngs_db <- dbConnect(SQLite(), dbname="realgar_data/realgar-omics.sqlite")

#SQLite for SNP-gene annotation ---
annot_db <- dbConnect(SQLite(), dbname="realgar_data/realgar-gwas-snp-annot.sqlite")
all_snps <- dbGetQuery(annot_db, "SELECT snp from gene_annot")[,"SNP"]

#SNP with eQTL annotation ---
eqtl_db <- dbConnect(SQLite(), dbname="realgar_data/realgar-eqtl-summary.sqlite")

#SQLite for gwas data from database ---
gwas_db <- dbConnect(SQLite(), dbname="realgar_data/realgar-gwas-hg38-normal.sqlite")
gwas_nominal_db <- dbConnect(SQLite(), dbname="realgar_data/realgar-gwas-hg38-nominal.sqlite")
gwas_genomewide_db <- dbConnect(SQLite(), dbname="realgar_data/realgar-gwas-hg38-genomewide.sqlite")

#SQLite for RNA-Seq counts ---
count_db <- dbConnect(SQLite(), dbname="transcriptomics/realgar-transcriptomics.sqlite")

## FUNCTIONS ------
## Get databases from gwas files and filter selected gene
get_query_db <- function(name, curr_gene){
  gene <- paste0("'",curr_gene,"'")
  query <- paste0("SELECT * FROM ",name," WHERE symbol=", gene)
  res <- dbSendQuery(gwas_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

get_query_nominal_db <- function(name, curr_gene){
  gene <- paste0("'",curr_gene,"'")
  query <- paste0("SELECT * FROM ",name," WHERE symbol=", gene)
  res <- dbSendQuery(gwas_nominal_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

get_query_genomewide_db <- function(name, curr_gene){
  gene <- paste0("'",curr_gene,"'")
  query <- paste0("SELECT * FROM ",name," WHERE symbol=", gene)
  res <- dbSendQuery(gwas_genomewide_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

## Get eQTL annotation for GWAS SNP
get_query_eqtl_db <- function(snp) {
  snp <- paste0("'", snp, "'")
  query <- paste0("SELECT * FROM eqtl WHERE SNP=", snp)
  res <- dbSendQuery(eqtl_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

## Annotate SNP data frame with eQTL annotation
get_eqtl_annot_to_GWAS_tb <- function(dat) {
  eqtl_tb <- do.call(rbind, lapply(dat$SNP, get_query_eqtl_db))
  dplyr::left_join(dat, eqtl_tb)
}


## Get snp data file. No longer use this function
get_snp <- function(name){
  query <- paste0("SELECT snp FROM ",name)
  res <- dbSendQuery(gwas_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(as.vector(data$snp))
}

## Get match for specific entered snp id
get_matches <- function(snp){
  rsid <- paste0("'",snp,"'")
  query <- paste0("SELECT snp,end,symbol FROM gene_annot WHERE snp = ",rsid)
  res <- dbSendQuery(annot_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

## Match snp id selected to gene location database and select gene nearest by distance to selected snp
# join_gene_snp <- function(all_matches){
#   genes <- unique(as.vector(all_matches$symbol))
#   data_table <- data.frame()
#   for (i in genes){
#     gene <- paste0("'",i,"'")
#     query <- paste0("SELECT DISTINCT symbol,start FROM gene_locations WHERE symbol = ", gene)
#     res <- dbSendQuery(gwas_db, query)
#     data <- dbFetch(res) 
#     dbClearResult(res)
#     data_table <- rbind(data_table,data)
#   }
#   data_table <- data_table[which(!duplicated(data_table$symbol)),]
#   all_matches <- merge(all_matches, data_table, by = "symbol")
#   all_matches$dist <- abs(all_matches$start - all_matches$end) # here, "end" is snp position, "start" is gene start
#   return(unique(all_matches$symbol[which(all_matches$dist==min(all_matches$dist))])) # choose the gene symbol whose start is the smallest absolute distance away
# }

join_gene_snp <- function(all_matches){
  # gene data is located from all_genes (realgar_data/gene_symbol_coords_hg38.RDS)
  genes <- unique(as.vector(all_matches$symbol))
  snp_pos <- unique(all_matches$end)
  data_table <- all_genes %>%
    dplyr::filter(symbol%in%genes) %>%
    dplyr::mutate(symbol=as.character(symbol)) %>%
    dplyr::mutate(gene_start=ifelse(strand==-1, End, Start), gene_end=ifelse(strand==-1, Start, End)) %>% # change start and end position for genes in minus strand
    dplyr::mutate(start_snp=abs(gene_start-snp_pos), end_snp=abs(gene_end-snp_pos)) # distance from snp to gene's boundaries
  start_snp_min <- min(data_table$start_snp)
  end_snp_min <- min(data_table$end_snp)
  if (start_snp_min>20000) { # if gene start is far away from snp, choose genes that are nearest to snp
    min_dist <- min(start_snp_min, end_snp_min)
    gene_sel_table <- data_table %>%
      dplyr::filter(start_snp==min_dist | end_snp==min_dist)
  } else { # if snp is within 20kb to gene start, choose this gene
    min_dist = start_snp_min
    gene_sel_table <- data_table %>%
      dplyr::filter(start_snp==min_dist)
  }
  gene_sel <- gene_sel_table$symbol
  return(gene_sel)
}

## Get counts for specific entered gene symbol
get_query_counts <- function(name, curr_gene){
  gene <- paste0("'",curr_gene,"'")
  query <- paste0("SELECT * FROM ",name," WHERE gene_symbol = ",gene)
  res <- dbSendQuery(count_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

## Get DE results for specific gene Ensembl id
get_query_DE <- function(name, curr_gene){
  gene <- paste0("'",curr_gene,"'")
  query <- paste0("SELECT * FROM ",name," WHERE Gene = ",gene)
  res <- dbSendQuery(count_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

## Get pheno info for specific gene Ensembl id
get_query_pheno <- function(name) {
  query <- paste0("SELECT * FROM ",name)
  res <- dbSendQuery(count_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}



