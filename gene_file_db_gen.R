# Getting gene information from cummeRbund is slow
# This file extracts all the gene information from
# cummeRbund databases and stores each file into it's
# own R object on disk.
rm(list=ls())
library("cummeRbund")
library("magrittr")
library("plyr")
library("parallel")
library("doParallel")

#Database to work on:
db_id <- "SRP043162"
db_dir <- paste0("databases/", db_id)

#Read Cufflinks DB and extract gene list
cuff.db = readCufflinks(dir=db_dir, rebuild=FALSE)
gene_list <- cuff.db %>% genes %>% featureNames %>% strsplit("[[:blank:]\"]+") %>% unlist
gene_list %>% length %>% cat #Output length to user

saveGene <- function(gene_name) {
    out.file <- paste0(db_dir,"/",gene_name,".rds")
    cuff.db <- readCufflinks(dir=db_dir, rebuild=FALSE)
    cuff.db %>% getGene(geneId=gene_name) %>% saveRDS(file=out.file)
}

#Setup parallel processing cluster on local machine
core_count <- detectCores()
cl <- makeCluster(core_count)
registerDoParallel(cl)

#Prep Cluster Processes
clusterEvalQ(cl, {
    library("cummeRbund")
    library("magrittr")
    library("plyr")
})
clusterExport(cl, varlist = c("db_id","db_dir"))
#clusterEvalQ(cl, {
#    cuff.db = readCufflinks(dir=db_dir, rebuild=FALSE)
#})
#Extract genes in parallel
adply(gene_list, 1, saveGene, .progress="text", .parallel = TRUE)

#Destroy cluster
stopCluster(cl)


