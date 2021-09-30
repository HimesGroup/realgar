library(dplyr)

#Read in files
#load info for gene tracks: gene locations, TFBS, SNPs, etc.
#tfbs <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/tfbs_for_app.RDS") #TFBS data from ENCODE - matched to gene ids using bedtools
snp <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/grasp_output_for_app.hg38.RDS") #SNP data from GRASP - matched to gene ids using bedtools
snp_eve <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/eve_data_realgar.hg38.RDS") #SNP data from EVE - was already in hg19 - matched to gene ids using bedtools 
snp_gabriel <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/gabriel_data_realgar.hg38.RDS") #SNP data from GABRIEL - lifted over from hg17 to hg19 - matched to gene ids using bedtools 
snp_fer <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/allerg_GWAS_data_realgar.hg38.RDS") #SNP data from Ferreira - already in hg19 - matched to gene ids using bedtools
snp_TAGC <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/TAGC_data_realgar.hg38.RDS") #SNP data from TAGC - already in hg19 - matched to gene ids using bedtools

#UKBiobank
ukbb_asthma <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/ukbb.asthma.sig.hg38.RDS") %>% dplyr::mutate(Phenotype = "Asthma")
ukbb_copd <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/ukbb.copd.sig.hg38.RDS") %>% dplyr::mutate(Phenotype = "COPD")
ukbb_aco <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/ukbb.aco.sig.hg38.RDS") %>% dplyr::mutate(Phenotype = "ACO")
snp_UKBB <- rbind(ukbb_asthma, ukbb_copd, ukbb_aco)

#Write bed files for all
#snp
snp <- snp %>% dplyr::select(-symbol)
write.table(snp, "/home/avantika/bedtools_GWAS/grasp_output_for_app.hg38.bed", row.names = F, quote=F, col.names = F, sep="\t")

#eve
snp_eve$start <- as.integer(snp_eve$end - 1)
snp_eve <- snp_eve %>% dplyr::select(-symbol)
write.table(snp_eve, "/home/avantika/bedtools_GWAS/eve_data_realgar.hg38.bed", row.names = F, quote=F, col.names = F, sep="\t")

#gabriel
snp_gabriel <- snp_gabriel %>% dplyr::select(-symbol)
write.table(snp_gabriel, "/home/avantika/bedtools_GWAS/gabriel_data_realgar.hg38.bed", row.names = F, quote=F, col.names = F, sep="\t")

#fer
snp_fer$start <- as.integer(snp_fer$end - 1)
snp_fer <- snp_fer %>% dplyr::select(-symbol)
write.table(snp_fer, "/home/avantika/bedtools_GWAS/allerg_GWAS_data_realgar.hg38.bed", row.names = F, quote=F, col.names = F, sep="\t")

#TAGC
snp_TAGC$start <- as.integer(snp_TAGC$end - 1)
snp_TAGC <- snp_TAGC %>% dplyr::select(-symbol)
write.table(snp_TAGC, "/home/avantika/bedtools_GWAS/TAGC_data_realgar.hg38.bed", row.names = F, quote=F, col.names = F, sep="\t")

#ukbb
snp_UKBB$start <- as.integer(snp_UKBB$end - 1)
write.table(snp_UKBB, "/home/avantika/bedtools_GWAS/ukbb_all_sig_hg38.bed", row.names = F, quote=F, col.names = F, sep="\t")

#Gene locations
gene_locations <- readRDS("/mnt/volume_nyc3_01/realgar_files/RDS_hg38/gene_locations_hg38.RDS")
gene_locations_df <- gene_locations %>% dplyr::mutate(start = ifelse(strand == -1, start, TSS-20000),
                                                      end = ifelse(strand == -1, TSS+20000, end), 
                                                      chromosome = paste0("chr",chromosome))
write.table(gene_locations_df, "/home/avantika/bedtools_GWAS/gene_locations_hg38_span.bed", row.names = F, quote=F, col.names = F,sep="\t")

#Map genes to snp using bedtools in the terminal
# 1. bedtools intersect -wa -wb -a ukbb_all_sig_hg38.bed -b gene_locations_hg38_span.bed | cut -f 1,2,3,4,5,6,7,8,9,10,18 | uniq > results/ukbb.all.sig.hg38.genes.bed
# 2. bedtools intersect -wa -wb -a grasp_output_for_app.hg38.bed -b gene_locations_hg38_span.bed | cut -f 1,2,3,4,5,6,14 | uniq > results/grasp_output_for_app.hg38.genes.bed
# 3. bedtools intersect -wa -wb -a eve_data_realgar.hg38.bed -b gene_locations_hg38_span.bed | cut -f 1,2,3,4,5,6,7,8,16 | uniq > results/eve_data_realgar.hg38.genes.bed
# 4. bedtools intersect -wa -wb -a gabriel_data_realgar.hg38.bed -b gene_locations_hg38_span.bed | cut -f 1,2,3,4,5,13 | uniq > results/gabriel_data_realgar.hg38.genes.bed
# 5. bedtools intersect -wa -wb -a allerg_GWAS_data_realgar.hg38.bed -b gene_locations_hg38_span.bed | cut -f 1,2,3,4,5,13 | uniq > results/allerg_GWAS_data_realgar.hg38.genes.bed
# 6. bedtools intersect -wa -wb -a TAGC_data_realgar.hg38.bed -b gene_locations_hg38_span.bed | cut -f 1,2,3,4,5,6,14 | uniq > results/TAGC_data_realgar.hg38.genes.bed
