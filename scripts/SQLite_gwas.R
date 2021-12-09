#Install package
require("RSQLite")
library(DBI)
library(dplyr)
library(dbplyr)
library(feather)
library(viridis)

#Read in files
#load info for gene tracks: gene locations, SNPs, etc.with gene mappings

#Gene locations (Run the following codes once)-----------
# gene_locations <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/gene_locations_hg38.RDS")
# gene_locations <- gene_locations %>% dplyr::mutate(chromosome = paste0("chr",chromosome)) 
# 
# #save for app
# genes <- gene_locations %>% 
#   dplyr::select(chromosome,start,end,strand,symbol) %>% 
#   dplyr::rename("Chromosome"="chromosome","Start"="start","End"="end") %>% 
#   unique()
# saveRDS(genes, "/home/avantika/gene_symbol_coords_hg38.RDS") # mv to /mnt/volume_nyc3_01/realgar_data/
#-----------

#SNP data from GRASP - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/grasp_output_for_app.hg38.RDS") %>%
  dplyr::rename(meta_P=p)

#SNP data from EVE - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_eve <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/eve_data_realgar.hg38.RDS")
snp_eve_all <- snp_eve %>% dplyr::select("chromosome","start","end","snp","meta_P", "symbol", "eqtl_lung", "eqtl_blood", "eqtl_muscle") %>% dplyr::filter(!is.na(meta_P))
snp_eve_ea <- snp_eve %>% dplyr::select("chromosome","start","end","snp","meta_P_EA", "symbol", "eqtl_lung", "eqtl_blood", "eqtl_muscle") %>% dplyr::rename(meta_P=meta_P_EA) %>% dplyr::filter(!is.na(meta_P))
snp_eve_aa <- snp_eve %>% dplyr::select("chromosome","start","end","snp","meta_P_AA", "symbol", "eqtl_lung", "eqtl_blood", "eqtl_muscle") %>% dplyr::rename(meta_P=meta_P_AA) %>% dplyr::filter(!is.na(meta_P))
snp_eve_la <- snp_eve %>% dplyr::select("chromosome","start","end","snp","meta_P_LAT", "symbol", "eqtl_lung", "eqtl_blood", "eqtl_muscle") %>% dplyr::rename(meta_P=meta_P_LAT) %>% dplyr::filter(!is.na(meta_P))


#SNP data from GABRIEL - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_gabriel <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/gabriel_data_realgar.hg38.RDS") %>%
  dplyr::rename(meta_P=P_ran)

#SNP data from Ferreira - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_fer <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/allerg_GWAS_data_realgar.hg38.RDS") %>%
  dplyr::rename(meta_P=PVALUE)

#SNP data from TAGC - lifted over from hg19 to hg38 - matched to gene symbols using bedtools 
snp_TAGC <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/TAGC_data_realgar.hg38.RDS")
snp_TAGC_multi <- snp_TAGC %>% dplyr::select("chromosome","start","end","snp","p_ran_multi", "symbol", "eqtl_lung", "eqtl_blood", "eqtl_muscle") %>%
  dplyr::rename(meta_P=p_ran_multi) %>% dplyr::filter(!is.na(meta_P))

snp_TAGC_euro <- snp_TAGC %>% dplyr::select("chromosome","start","end","snp","p_ran_euro", "symbol", "eqtl_lung", "eqtl_blood", "eqtl_muscle") %>%
  dplyr::rename(meta_P=p_ran_euro) %>% dplyr::filter(!is.na(meta_P))

#SNP data from UKBB - asthma, copd, aco results combined - matched to gene symbols using bedtools
snp_UKBB_asthma <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/ukbb.asthma.sig.hg38.RDS") %>%
  dplyr::rename(meta_P=P) %>% dplyr::filter(!is.na(meta_P))
snp_UKBB_copd <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/ukbb.copd.sig.hg38.RDS") %>%  dplyr::rename(meta_P=P) %>% dplyr::filter(!is.na(meta_P))
snp_UKBB_aco <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/ukbb.aco.sig.hg38.RDS") %>%
  dplyr::rename(meta_P=P) %>% dplyr::filter(!is.na(meta_P))

# add color columns
breaks <- c(seq(0,8,by=0.001), Inf) # this sets max universally at 8 (else highest one OF THE SUBSET would be the max)

addcolor_func <- function(dat) {
  dat <- dplyr::mutate(dat, neg_log_p = -log10(meta_P))
  dat$color <- inferno(8002)[as.numeric(cut(dat$neg_log_p, breaks = breaks))]
  return(dat)
}

#GRASP SNP
snp <- addcolor_func(snp)

#SNP_EVE
snp_eve_all <- addcolor_func(snp_eve_all)
snp_eve_aa <- addcolor_func(snp_eve_aa)
snp_eve_ea <- addcolor_func(snp_eve_ea)
snp_eve_la <- addcolor_func(snp_eve_la)

#SNP gabriel
snp_gabriel <- addcolor_func(snp_gabriel)

#SNP fer
snp_fer <- addcolor_func(snp_fer)

#SNP TAGC
snp_TAGC_multi <- addcolor_func(snp_TAGC_multi)
snp_TAGC_euro <- addcolor_func(snp_TAGC_euro)

#UKBB 
snp_UKBB_asthma <- addcolor_func(snp_UKBB_asthma)
snp_UKBB_copd <- addcolor_func(snp_UKBB_copd)
snp_UKBB_aco <- addcolor_func(snp_UKBB_aco)

# apply p-value thresholds
datsel_func <- function(dat, pval_thr=0.05) {
  dat %>% dplyr::filter(meta_P<pval_thr)
}

pval_thr1 = 10^(-5)
#GRASP SNP - nominal significance
snp_nominal <- datsel_func(snp, pval_thr = pval_thr1)

#SNP_EVE - nominal significance
snp_eve_all_nominal <- datsel_func(snp_eve_all, pval_thr = pval_thr1)
snp_eve_aa_nominal <- datsel_func(snp_eve_aa, pval_thr = pval_thr1)
snp_eve_ea_nominal <- datsel_func(snp_eve_ea, pval_thr = pval_thr1)
snp_eve_la_nominal <- datsel_func(snp_eve_la, pval_thr = pval_thr1)

#SNP gabriel - nominal significance
snp_gabriel_nominal <- datsel_func(snp_gabriel, pval_thr = pval_thr1)

#SNP fer - nominal significance
snp_fer_nominal <- datsel_func(snp_fer, pval_thr = pval_thr1)

#SNP TAGC - nominal significance
snp_TAGC_multi_nominal <- datsel_func(snp_TAGC_multi, pval_thr = pval_thr1)
snp_TAGC_euro_nominal <- datsel_func(snp_TAGC_euro, pval_thr = pval_thr1)

#UKBB  - nominal significance
snp_UKBB_asthma_nominal <- datsel_func(snp_UKBB_asthma, pval_thr = pval_thr1)
snp_UKBB_copd_nominal <- datsel_func(snp_UKBB_copd, pval_thr = pval_thr1)
snp_UKBB_aco_nominal <- datsel_func(snp_UKBB_aco, pval_thr = pval_thr1)


pval_thr2 = 5*10^(-8)
#GRASP SNP - gwas significance
snp_gwas <- datsel_func(snp, pval_thr = pval_thr2)

#SNP_EVE - gwas significance
snp_eve_all_gwas <- datsel_func(snp_eve_all, pval_thr = pval_thr2)
snp_eve_aa_gwas <- datsel_func(snp_eve_aa, pval_thr = pval_thr2)
snp_eve_ea_gwas <- datsel_func(snp_eve_ea, pval_thr = pval_thr2)
snp_eve_la_gwas <- datsel_func(snp_eve_la, pval_thr = pval_thr2)

#SNP gabriel - gwas significance
snp_gabriel_gwas <- datsel_func(snp_gabriel, pval_thr = pval_thr2)

#SNP fer - gwas significance
snp_fer_gwas <- datsel_func(snp_fer, pval_thr = pval_thr2)

#SNP TAGC - gwas significance
snp_TAGC_multi_gwas <- datsel_func(snp_TAGC_multi, pval_thr = pval_thr2)
snp_TAGC_euro_gwas <- datsel_func(snp_TAGC_euro, pval_thr = pval_thr2)

#UKBB  - gwas significance
snp_UKBB_asthma_gwas <- datsel_func(snp_UKBB_asthma, pval_thr = pval_thr2)
snp_UKBB_copd_gwas <- datsel_func(snp_UKBB_copd, pval_thr = pval_thr2)
snp_UKBB_aco_gwas <- datsel_func(snp_UKBB_aco, pval_thr = pval_thr2)

#color tfbs based on binding score - used in tracks
#create color scheme based on encode binding score & snp p-values
#tfbs$color <- inferno(50)[as.numeric(cut(tfbs$score,breaks = 50))]

#Put it in database
#db = dbConnect(SQLite(), dbname="/mnt/volume_nyc3_01/realgar_data/realgar-gwas.sqlite")
db = dbConnect(SQLite(), dbname="realgar-gwas-hg38-normal.sqlite") #sudo mv to /mnt/volume_nyc3_01/realgar_data/

#Check table
dbListTables(db)

#Add in database
dbWriteTable(conn=db, name="snp", snp, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_all", snp_eve_all, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_aa", snp_eve_aa, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_ea", snp_eve_ea, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_la", snp_eve_la, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_gabriel", snp_gabriel, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_fer", snp_fer, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_TAGC_multi", snp_TAGC_multi, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_TAGC_euro", snp_TAGC_euro, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_UKBB_asthma", snp_UKBB_asthma, row.names=F, overwrite=T)
#dbWriteTable(conn=db, name="tfbs", tfbs, row.names=F)
#dbWriteTable(conn=db, name="gene_locations", gene_locations, row.names=F)

#Check table
dbListTables(db)

#Disconnect
dbDisconnect(db)


#Put it in database - nominal significance
#db = dbConnect(SQLite(), dbname="/mnt/volume_nyc3_01/realgar_data/realgar-gwas.sqlite")
db = dbConnect(SQLite(), dbname="realgar-gwas-hg38-nominal.sqlite") #sudo mv to /mnt/volume_nyc3_01/realgar_data/

#Check table
dbListTables(db)

#Add in database
dbWriteTable(conn=db, name="snp", snp_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_all", snp_eve_all_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_aa", snp_eve_aa_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_ea", snp_eve_ea_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_la", snp_eve_la_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_gabriel", snp_gabriel_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_fer", snp_fer_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_TAGC_multi", snp_TAGC_multi_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_TAGC_euro", snp_TAGC_euro_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_UKBB_asthma", snp_UKBB_asthma_nominal, row.names=F, overwrite=T)
#dbWriteTable(conn=db, name="tfbs", tfbs, row.names=F)
#dbWriteTable(conn=db, name="gene_locations", gene_locations, row.names=F)

#Check table
dbListTables(db)

#Disconnect
dbDisconnect(db)

#Put it in database - nominal significance
#db = dbConnect(SQLite(), dbname="/mnt/volume_nyc3_01/realgar_data/realgar-gwas.sqlite")
db = dbConnect(SQLite(), dbname="realgar-gwas-hg38-genomewide.sqlite") #sudo mv to /mnt/volume_nyc3_01/realgar_data/

#Check table
dbListTables(db)

#Add in database
dbWriteTable(conn=db, name="snp", snp_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_all", snp_eve_all_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_aa", snp_eve_aa_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_ea", snp_eve_ea_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_eve_la", snp_eve_la_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_gabriel", snp_gabriel_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_fer", snp_fer_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_TAGC_multi", snp_TAGC_multi_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_TAGC_euro", snp_TAGC_euro_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_UKBB_asthma", snp_UKBB_asthma_gwas, row.names=F, overwrite=T)
#dbWriteTable(conn=db, name="tfbs", tfbs, row.names=F)
#dbWriteTable(conn=db, name="gene_locations", gene_locations, row.names=F)

#Check table
dbListTables(db)

#Disconnect
dbDisconnect(db)

