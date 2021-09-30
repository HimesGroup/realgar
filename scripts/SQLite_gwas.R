#Install package
require("RSQLite")
library(DBI)
library(dplyr)
library(dbplyr)
library(feather)
library(viridis)

#Read in files
#load info for gene tracks: gene locations, SNPs, etc.with gene mappings

#SNP data from GRASP - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp <- read.delim("/mnt/volume_nyc3_01/realgar_files/bedtools_GWAS/results/grasp_output_for_app.hg38.genes.bed", 
                  header=F, col.names = c("chromosome","start","end","snp","p","pmid","symbol")) 

#SNP data from EVE - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_eve <- read.delim("/mnt/volume_nyc3_01/realgar_files/bedtools_GWAS/results/eve_data_realgar.hg38.genes.bed",
                      header=F, col.names = c("chromosome","start","end","snp","meta_P","meta_P_EA","meta_P_AA","meta_P_LAT","symbol")) 

#SNP data from GABRIEL - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_gabriel <- read.delim("/mnt/volume_nyc3_01/realgar_files/bedtools_GWAS/results/gabriel_data_realgar.hg38.genes.bed",
                          header=F, col.names = c("chromosome","start","end","snp","P_ran","symbol")) 

#SNP data from Ferreira - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_fer <- read.delim("/mnt/volume_nyc3_01/realgar_files/bedtools_GWAS/results/allerg_GWAS_data_realgar.hg38.genes.bed",
                      header=F, col.names = c("chromosome","start","end","snp","PVALUE","symbol")) 

#SNP data from TAGC - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_TAGC <- read.delim("/mnt/volume_nyc3_01/realgar_files/bedtools_GWAS/results/TAGC_data_realgar.hg38.genes.bed",
                       header=F, col.names = c("chromosome","start","end","snp","p_ran_multi","p_ran_euro","symbol")) 

#SNP data from UKBB - asthma, copd, aco results combined - matched to gene symbols using bedtools
snp_UKBB <- read.delim("/mnt/volume_nyc3_01/realgar_files/bedtools_GWAS/results/ukbb.all.sig.hg38.genes.bed",
                       header=F, col.names = c("chromosome","start","end","snp","ref","alt","P","OR","SE","Phenotype","symbol")) 

#Gene locations
gene_locations <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/gene_locations_hg38.RDS")
gene_locations <- gene_locations %>% dplyr::mutate(chromosome = paste0("chr",chromosome)) 

#save for app
genes <- gene_locations %>% 
  dplyr::select(chromosome,start,end,strand,symbol) %>% 
  dplyr::rename("Chromosome"="chromosome","Start"="start","End"="end") %>% 
  unique()
saveRDS(genes, "/home/avantika/gene_symbol_coords_hg38.RDS") # mv to /mnt/volume_nyc3_01/realgar_data/

breaks <- c(seq(0,8,by=0.001), Inf) # this sets max universally at 8 (else highest one OF THE SUBSET would be the max)

#SNP
snp <- dplyr::mutate(snp, neg_log_p = -log10(p))
snp$color <- inferno(8002)[as.numeric(cut(snp$neg_log_p, breaks = breaks))]


#SNP_EVE
snp_eve <- dplyr::mutate(snp_eve, neg_log_meta_p = -log10(meta_P),
                         neg_log_meta_p_aa = -log10(meta_P_AA),
                         neg_log_meta_p_ea = -log10(meta_P_EA),
                         neg_log_meta_p_lat = -log10(meta_P_LAT))
snp_eve$color_meta_P <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p,breaks = breaks))]
snp_eve$color_meta_P_AA <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p_aa,breaks = breaks))]
snp_eve$color_meta_P_EA <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p_ea,breaks = breaks))]
snp_eve$color_meta_P_LAT <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p_lat,breaks = breaks))]

#SNP gabriel
snp_gabriel <- dplyr::mutate(snp_gabriel, neg_log_p = -log10(P_ran))
snp_gabriel$color <- inferno(8002)[as.numeric(cut(snp_gabriel$neg_log_p, breaks = breaks))]

#SNP fer
snp_fer <- dplyr::mutate(snp_fer, neg_log_p = -log10(PVALUE))
snp_fer$color <- inferno(8002)[as.numeric(cut(snp_fer$neg_log_p, breaks = breaks))]

#SNP TAGC
snp_TAGC <- dplyr::mutate(snp_TAGC, neg_log_p_multi = -log10(p_ran_multi),
                          neg_log_p_euro = -log10(p_ran_euro))
snp_TAGC$color_p_ran_multi <- inferno(8002)[as.numeric(cut(snp_TAGC$p_ran_multi,breaks = breaks))]
snp_TAGC$color_p_ran_euro <- inferno(8002)[as.numeric(cut(snp_TAGC$p_ran_euro,breaks = breaks))]

#UKBB 
snp_UKBB <- dplyr::mutate(snp_UKBB, neg_log_p = -log10(P))
snp_UKBB$color <- inferno(8002)[as.numeric(cut(snp_UKBB$neg_log_p, breaks = breaks))]

#color tfbs based on binding score - used in tracks
#create color scheme based on encode binding score & snp p-values
#tfbs$color <- inferno(50)[as.numeric(cut(tfbs$score,breaks = 50))]

#Put it in database
#db = dbConnect(SQLite(), dbname="/mnt/volume_nyc3_01/realgar_data/realgar-gwas.sqlite")
db = dbConnect(SQLite(), dbname="/home/avantika/realgar-gwas-hg38.sqlite") #sudo mv to /mnt/volume_nyc3_01/realgar_data/

#Check table
dbListTables(db)

#Add in database
dbWriteTable(conn=db, name="snp", snp, row.names=F)
dbWriteTable(conn=db, name="snp_eve", snp_eve, row.names=F)
dbWriteTable(conn=db, name="snp_gabriel", snp_gabriel, row.names=F)
dbWriteTable(conn=db, name="snp_fer", snp_fer, row.names=F)
dbWriteTable(conn=db, name="snp_TAGC", snp_TAGC, row.names=F)
dbWriteTable(conn=db, name="snp_UKBB", snp_UKBB, row.names=F)
#dbWriteTable(conn=db, name="tfbs", tfbs, row.names=F)
#dbWriteTable(conn=db, name="gene_locations", gene_locations, row.names=F)

#Check table
dbListTables(db)

#Disconnect
dbDisconnect(db)

