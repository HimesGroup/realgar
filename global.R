library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(feather, quietly = T)
library(dplyr, quietly = T)
library(tidyr, quietly = T)
library(data.table, quietly = T)
library(tidyverse, quietly = T)
library(ggplot2, quietly = T)
library(cowplot, quietly = T)
library(lattice, quietly = T)
library(stringr, quietly = T)
library(viridis, quietly = T) 
library(scales, quietly = T)
library(stringi, quietly = T) # use function stri_dup to add white space
library(DT, quietly = T) 
library(magrittr, quietly = T)
library(karyoploteR, quietly = T) #Need R > 3.5
library(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = T)
source("utilities/sql_queries.R")
source("utilities/meta.R")
source("utilities/comb_pval.R")
#source("utilities/name_convert.R")
source("utilities/forestplots.R")
source("utilities/karyoplot.R")
#source("utilities/gviz_gene_track.R")


##############################################
## Choices for tissue, asthma and treatment ##
##############################################

#Tissue types
tissue_choices <-c("Airway smooth muscle"="ASM", "Bronchial epithelium"="BE", "Large airway epithelium"="LAE", "Lens epithelium" = "LEC",
                   "BEAS-2B" = "BEAS-2B","Nasal epithelium"="NE","Small airway epithelium"="SAE", "Trachea"="trachea", "Buccal mucosa"="Buccal", 
                   "Whole lung"="Lung","Skeletal muscle myotubes"="myotubes","CD4"="CD4", "CD8"="CD8", "MCF10A-Myc" = "MCF10A-Myc","CD3" = "CD3", 
                   "A549" = "A549","Lymphoblastoid cell" = "LCL","Macrophage" = "MACRO", "Alveolar macrophages" = "AM", 
                   "Peripheral blood mononuclear cell"="PBMC","White blood cell"="WBC","Whole blood"="Blood",
                   "Lymphoblastic leukemia cell" = "chALL","Osteosarcoma U2OS cell" = "U2OS", "WI38 fibroblast"="WI38")

tissue_selected <- c("Airway smooth muscle"="ASM", "Bronchial epithelium"="BE")

#Disease types
asthma_choices <- c("Allergic asthma"="allergic_asthma", "Asthma"="asthma",
                    "Fatal asthma"="fatal_asthma", "Mild to moderate asthma"="mild_to_moderate_asthma","Severe asthma"="severe_asthma",
                    "Mild asthma with rhinitis"="rhinitis_mild_asthma","Severe asthma with rhinitis"="rhinitis_severe_asthma",
                    "Non-allergic asthma"="non_allergic_asthma")

asthma_selected <- c("Asthma"="asthma", "Severe asthma"="severe_asthma")

#Treatment choices
treatment_choices <- c("β2-agonist"="BA", "Glucocorticoid" = "GC",
                       "Phosphodiesterase inhibitor"="PDE","Vitamin D"="vitD")
                        #"Smoking"="smoking", "E-cigarette" = "ecig")
#smoking choices
smoking_choices <- c("Cigarette"="cig", "E-cigarette" = "ecig")

# pollutant choices
pollutant_choices <- c("Polycyclic aromatic hydrocarbons"="PAH", "Particulate matter"="pollutant")

#exposure choices
treatment_choices <- c(treatment_choices, smoking_choices, pollutant_choices)
#treatment_choices <- c(treatment_choices, smoking_choices)
treatment_selected <- c("Glucocorticoid" = "GC", "Cigarette"="cig")

#experiment choices
experiment_choices=c("Cell-based assay"="invitro","Human response study"="invivo")
experiment_selected=c("Cell-based assay"="invitro")                 

#GWAS options
gwas_choices <- c("Ferreira"="snp_fer_subs","GABRIEL"="snp_gabriel_subs","GRASP"="snp_subs",
                  "EVE all subjects"="snp_eve_all_subs", "EVE African Americans"="snp_eve_aa_subs", "EVE European Americans"="snp_eve_ea_subs", "EVE Latinos"="snp_eve_la_subs",
                  "TAGC Multiancestry"="snp_TAGC_multi_subs", "TAGC European ancestry"="snp_TAGC_euro_subs",
                  "UKBiobank Asthma"="snp_UKBB_asthma_subs", "UKBiobank COPD"="snp_UKBB_copd_subs", "UKBiobank ACO"="snp_UKBB_aco_subs")
#pval_choices = c("0.05"="normal", "1x10<sup>-5</sup>"="nominal","5x10<sup>-8</sup>"="genomewide")
pval_select <- c("0.05"="normal")
names(pval_select) = paste0("0.05", stri_dup(intToUtf8(160), 18))

pval_for_select <-  tibble(value=c("normal","nominal","genomewide"),
                               label=c("0.05","1x10-5","5x10-8"),
                               html=c("0.05", "1x10<sup>-5</sup>" , "5x10<sup>-8</sup>"))
pval_for_select[[2]][1] = pval_for_select[[3]][1] <- paste0("0.05", stri_dup(intToUtf8(160), 18))
pval_for_select[[2]][2] = paste0("1x10-5", stri_dup(intToUtf8(160), 18))
pval_for_select[[2]][3] = paste0("5x10-8", stri_dup(intToUtf8(160), 18))

#Gene list
# all_genes <- read_feather("realgar_data/gene_list.feather")
# gene_list <- as.vector(all_genes$V1)
# rm(all_genes)

#all_genes <- readRDS("/mnt/volume_nyc3_01/realgar_data/gene_symbol_POS.RDS")
#all_genes <- readRDS("/mnt/volume_nyc3_01/realgar_files/gene_symbol_coords_hg19.RDS")
all_genes <- readRDS("realgar_data/gene_symbol_coords_hg38.RDS")
gene_list <- as.vector(all_genes$symbol)

#Gene choices (Not in use)
#genec <- read_feather("realgar_data/Sig_gene_list.feather")
#gene_choices <- as.vector(genec$V1)
#rm(genec)

####################
## READ IN FILES ##
####################

# load descriptions of all gene expression and GWAS datasets
Alldata_Info <- read_feather("realgar_data/Microarray_data_infosheet_latest_R.feather")
#Alldata_Info <- read.csv("realgar_data/Microarray_data_infosheet_latest_R.csv")

#then split off into gene expression and GWAS dataset info - else forest plot text columns get messed up
GWAS_Dataset_Info <- Alldata_Info[which(Alldata_Info$App == "GWAS"),]
Dataset_Info <- Alldata_Info[which(!(Alldata_Info$App == "GWAS")),]

#Remove big data ---
rm(Alldata_Info)

Dataset_Info$PMID <- as.character(Dataset_Info$PMID) #else next line does not work
Dataset_Info[is.na(Dataset_Info$PMID), "PMID"] <- ""
Dataset_Info$Report <- as.character(c("QC"))

##BA_PDE dataset ---
BA_PDE_Info <- Dataset_Info %>% dplyr::filter(Asthma == "BA_PDE")

##ChIP-Seq dataset
chipseq_dataset <- read_feather("realgar_data/realgar_ChIPSeq_datasets.feather")

####################
## GWAS SNP data ##
####################

#load info for gene tracks: gene locations, TFBS, SNPs, etc.

#from feather files ---
# chrom_bands <- read_feather("realgar_data/chrom_bands.feather") #chromosome band info for ideogram - makes ideogram load 25 seconds faster

####################
## GWAS SNP data ##
####################

# load GR-binding sites
GRbinding <- read_feather("realgar_data/GR_binding_sites_sig.feather")


###########################
## Transcriptomic data ##
###########################

#Load data files - gene names and dataset info
# "lcte" appended to beginning of filename stands for "lung cell transcriptome explorer"
sras <- read_feather("transcriptomics/asthmagenes_deseq2/lcte_dataset_info_asm.feather") %>% tibble::as_tibble()
all_genes_te <- read_feather("transcriptomics/lcte_gene_names.feather") %>% tibble::as_tibble()
unfiltered_genes <- read_feather("transcriptomics/lcte_sleuth_unfiltered_genes.feather") %>% tibble::as_tibble()

rnaseq_choices <- c("SRP033351 (airway smooth muscle)" = "SRP033351", "SRP043162 (airway smooth muscle)" = "SRP043162","SRP098649 (airway smooth muscle)" = "SRP098649",
                    "SRP157114 (bronchial epithelium)" = "SRP157114", "SRP237772 (bronchial epithelium)" = "SRP237772", "SRP216947 (bronchial epithelium)" = "SRP216947",
                    "SRP005411 (small airway epithelium)"="SRP005411", "SRP277255 (bronchial epithelium)"="SRP277255")

#Deseq2 results : log2FC, padj and conditions- for datatable 
# de <- list()
# #Deseq2 count results - by gene for plots
# tpms <- list()
# for (study in unname(rnaseq_choices)) {
#   de[[study]] <- read_feather(paste0("transcriptomics/asthmagenes_deseq2/", study, "/", study, "_de.feather")) %>% tibble::as_tibble()
#   tpms[[study]]  <- read_feather(paste0("transcriptomics/asthmagenes_deseq2/", study, "/", study, "_pheno+counts_updated.feather")) %>% tibble::as_tibble()
# }

# make a list of gene symbols in all datasets for checking whether gene symbol entered is valid - used later on
#deseq2_filtered_genes <- unlist(lapply(unname(rnaseq_choices), function(study)de[[study]][,"gene_symbol"])) %>% unique()
deseq2_filtered_genes_tb <- read_feather("transcriptomics/lcte_sleuth_filtered_genes.feather")
deseq2_filtered_genes <- deseq2_filtered_genes_tb$gene_symbol
