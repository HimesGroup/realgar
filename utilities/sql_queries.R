library(DBI, quietly = T)
library(RSQLite, quietly=T)

#SQLite for omics database ---
ngs_db <- dbConnect(SQLite(), dbname="realgar_data/realgar-omics.sqlite")

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
  query <- paste0("SELECT * FROM ",name," WHERE symbol = ",gene)
  res <- dbSendQuery(gwas_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

get_query_nominal_db <- function(name, curr_gene){
  gene <- paste0("'",curr_gene,"'")
  query <- paste0("SELECT * FROM ",name," WHERE symbol = ",gene)
  res <- dbSendQuery(gwas_nominal_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

get_query_genomewide_db <- function(name, curr_gene){
  gene <- paste0("'",curr_gene,"'")
  query <- paste0("SELECT * FROM ",name," WHERE symbol = ",gene)
  res <- dbSendQuery(gwas_genomewide_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}



## Get snp data file
get_snp <- function(name){
  query <- paste0("SELECT snp FROM ",name)
  res <- dbSendQuery(gwas_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(as.vector(data$snp))
}

## Get match for specific entered snp id
get_matches <- function(snp, name){
  rsid <- paste0("'",snp,"'")
  query <- paste0("SELECT snp,end,symbol FROM ",name," WHERE snp = ",rsid)
  res <- dbSendQuery(gwas_db, query)
  data <- dbFetch(res) 
  dbClearResult(res)
  return(data)
}

## Match snp id selected to gene location database and select gene nearest by distance to selected snp
join_gene_snp <- function(all_matches){
  genes <- unique(as.vector(all_matches$symbol))
  data_table <- data.frame()
  for (i in genes){
    gene <- paste0("'",i,"'")
    query <- paste0("SELECT DISTINCT symbol,start FROM gene_locations WHERE symbol = ", gene)
    res <- dbSendQuery(gwas_db, query)
    data <- dbFetch(res) 
    dbClearResult(res)
    data_table <- rbind(data_table,data)
  }
  data_table <- data_table[which(!duplicated(data_table$symbol)),]
  all_matches <- merge(all_matches, data_table, by = "symbol")
  all_matches$dist <- abs(all_matches$start - all_matches$end) # here, "end" is snp position, "start" is gene start
  return(unique(all_matches$symbol[which(all_matches$dist==min(all_matches$dist))])) # choose the gene symbol whose start is the smallest absolute distance away
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



