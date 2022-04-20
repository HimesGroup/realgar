#Install package
require("RSQLite")
library(DBI)
library(dplyr)
library(dbplyr)
library(feather)
library(viridis)

# Read in GWAS SNP files
# Original files are save under the lab server: /projects/AsthmaApp/REALGAR/GWAS/hg38_annot/[study].hg38.realgar.final.bed.
# The files are transferred to the R server: /mnt/volume_nyc3_01/realgar_files/hg38_annot

readbed_func <- function(bed_fn) {
  dat <- read.table(bed_fn, header=T, sep="\t", comment.char = "")
  names(dat)[1] <- "chromosome"
  return(dat)
}

#SNP data from GRASP
snp_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/GRASP.hg38.realgar.final.bed"
snp <- readbed_func(snp_fn) %>%
  dplyr::rename(meta_P=pval) %>%
  dplyr::select(chromosome, end, SNP, symbol, meta_P, pmid) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)

#SNP data from EVE
snp_eve_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/EVE.hg38.realgar.final.bed"
snp_eve <- readbed_func(snp_eve_fn)

snp_eve_all <- snp_eve %>%
  dplyr::select(chromosome, end, SNP, symbol, meta.p) %>%
  dplyr::rename(meta_P=meta.p) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


snp_eve_ea <- snp_eve %>%
  dplyr::select(chromosome, end, SNP, symbol, EA.p) %>%
  dplyr::rename(meta_P=EA.p) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


snp_eve_aa <- snp_eve %>%
  dplyr::select(chromosome, end, SNP, symbol, AA.p) %>%
  dplyr::rename(meta_P=AA.p) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


snp_eve_la <- snp_eve %>%
  dplyr::select(chromosome, end, SNP, symbol, LA.p) %>%
  dplyr::rename(meta_P=LA.p) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


#SNP data from GABRIEL
snp_gabriel_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/GABRIEL.hg38.realgar.final.bed"
snp_gabriel <- readbed_func(snp_gabriel_fn) %>%
  dplyr::select(chromosome, end, SNP, symbol, meta.p) %>%
  dplyr::rename(meta_P=meta.p) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


#SNP data from Ferreira
snp_fer_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/Ferreira.hg38.realgar.final.bed"
snp_fer <- readbed_func(snp_fer_fn) %>%
  dplyr::select(chromosome, end, SNP, symbol, pval) %>%
  dplyr::rename(meta_P=pval) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


#SNP data from TAGC
snp_TAGC_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/TAGC.hg38.realgar.final.bed"
snp_TAGC <- readbed_func(snp_TAGC_fn)

snp_TAGC_multi <- snp_TAGC %>%
  dplyr::select(chromosome, end, SNP, symbol, pval_multi) %>%
  dplyr::rename(meta_P=pval_multi) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


snp_TAGC_euro <- snp_TAGC %>%
  dplyr::select(chromosome, end, SNP, symbol, pval_euro) %>%
  dplyr::rename(meta_P=pval_euro) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


#SNP data from UKBB
snp_UKBB_asthma_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/ukbb.asthma.hg38.realgar.final.bed"
snp_UKBB_asthma <- readbed_func(snp_UKBB_asthma_fn) %>%
  dplyr::select(chromosome, end, SNP, symbol, P) %>%
  dplyr::rename(meta_P=P) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


snp_UKBB_copd_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/ukbb.copd.hg38.realgar.final.bed"
snp_UKBB_copd <- readbed_func(snp_UKBB_copd_fn) %>%
  dplyr::select(chromosome, end, SNP, symbol, P) %>%
  dplyr::rename(meta_P=P) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


snp_UKBB_aco_fn <- "/mnt/volume_nyc3_01/realgar_files/hg38_annot/ukbb.aco.hg38.realgar.final.bed"
snp_UKBB_aco <- readbed_func(snp_UKBB_aco_fn) %>%
  dplyr::select(chromosome, end, SNP, symbol, P) %>%
  dplyr::rename(meta_P=P) %>%
  dplyr::filter(!is.na(meta_P)&meta_P!="") %>%
  unique() %>%
  dplyr::arrange(chromosome, end)


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

#Put it in database
db = dbConnect(SQLite(), dbname="/mnt/volume_nyc3_01/realgar_files/hg38_annot/sqilte_results/realgar-gwas-hg38-normal.sqlite") #sudo mv to /mnt/volume_nyc3_01/realgar_data/

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
dbWriteTable(conn=db, name="snp_UKBB_copd", snp_UKBB_copd, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_UKBB_aco", snp_UKBB_aco, row.names=F, overwrite=T)

#Check table
dbListTables(db)

#Disconnect
dbDisconnect(db)


#Put it in database - nominal significance
db = dbConnect(SQLite(), dbname="/mnt/volume_nyc3_01/realgar_files/hg38_annot/sqilte_results/realgar-gwas-hg38-nominal.sqlite") #sudo mv to /mnt/volume_nyc3_01/realgar_data/

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
dbWriteTable(conn=db, name="snp_UKBB_copd", snp_UKBB_copd_nominal, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_UKBB_aco", snp_UKBB_aco_nominal, row.names=F, overwrite=T)

#Check table
dbListTables(db)

#Disconnect
dbDisconnect(db)

#Put it in database - nominal significance
db = dbConnect(SQLite(), dbname="/mnt/volume_nyc3_01/realgar_files/hg38_annot/sqilte_results/realgar-gwas-hg38-genomewide.sqlite") #sudo mv to /mnt/volume_nyc3_01/realgar_data/

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
dbWriteTable(conn=db, name="snp_UKBB_copd", snp_UKBB_copd_gwas, row.names=F, overwrite=T)
dbWriteTable(conn=db, name="snp_UKBB_aco", snp_UKBB_aco_gwas, row.names=F, overwrite=T)

#Check table
dbListTables(db)

#Disconnect
dbDisconnect(db)
