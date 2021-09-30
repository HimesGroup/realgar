
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(feather, quietly = T)
library(dplyr, quietly = T)
library(data.table, quietly = T)
library(ggplot2, quietly = T)
library(cowplot, quietly = T)
library(lattice, quietly = T)
library(stringr, quietly = T)
library(viridis, quietly = T) 
library(scales, quietly = T)
library(DT, quietly = T) 
library(magrittr, quietly = T)
source("utilities/sql_queries.R")
source("utilities/meta.R")
source("utilities/comb_pval.R")
source("utilities/name_convert.R")
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
                   "Lymphoblastic leukemia cell" = "chALL","Osteosarcoma U2OS cell" = "U2OS")

#Disease types
asthma_choices <- c("Allergic asthma"="allergic_asthma", "Asthma"="asthma",
                    "Fatal asthma"="fatal_asthma", "Mild to moderate asthma"="mild_to_moderate_asthma","Severe asthma"="severe_asthma",
                    "Mild asthma with rhinitis"="rhinitis_mild_asthma","Severe asthma with rhinitis"="rhinitis_severe_asthma",
                    "Non-allergic asthma"="non_allergic_asthma")

#Treatment choices
treatment_choices <- c("β2-agonist"="BA", "Glucocorticoid" = "GC",
                       "Phosphodiesterase inhibitor"="PDE","Vitamin D"="vitD")
                        #"Smoking"="smoking", "E-cigarette" = "ecig")
#smoking choices
smoking_choices <- c("Cigarette"="cig", "E-cigarette" = "ecig")
                      
#GWAS options
gwas_choices <- c("EVE"="snp_eve_subs","Ferreira"="snp_fer_subs","GABRIEL"="snp_gabriel_subs","GRASP"="snp_subs","TAGC"="snp_TAGC_subs","UKBiobank"="snp_UKBB_subs")

#Gene list
# all_genes <- read_feather("realgar_data/gene_list.feather")
# gene_list <- as.vector(all_genes$V1)
# rm(all_genes)

#all_genes <- readRDS("/mnt/volume_nyc3_01/realgar_data/gene_symbol_POS.RDS")
#all_genes <- readRDS("/mnt/volume_nyc3_01/realgar_files/gene_symbol_coords_hg19.RDS")
all_genes <- readRDS("realgar_data/gene_symbol_coords_hg38.RDS")
gene_list <- as.vector(all_genes$symbol)

#Gene choices
genec <- read_feather("realgar_data/Sig_gene_list.feather")
gene_choices <- as.vector(genec$V1)
rm(genec)

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

####################
## GWAS SNP data ##
####################

#load info for gene tracks: gene locations, TFBS, SNPs, etc.

#from feather files ---
chrom_bands <- read_feather("realgar_data/chrom_bands.feather") #chromosome band info for ideogram - makes ideogram load 25 seconds faster

###########################
## Transcriptomic data ##
###########################

#Load data files - gene names and dataset info
# "lcte" appended to beginning of filename stands for "lung cell transcriptome explorer"
sras <- read_feather("transcriptomics/asthmagenes_deseq2/lcte_dataset_info_asm.feather") %>% tibble::as_tibble()
all_genes_te <- read_feather("transcriptomics/lcte_gene_names.feather") %>% tibble::as_tibble()
unfiltered_genes <- read_feather("transcriptomics/lcte_sleuth_unfiltered_genes.feather") %>% tibble::as_tibble()

#Deseq2 results : log2FC, padj and conditions- for datatable 
de <- list()
de[["SRP033351"]] <- read_feather("transcriptomics/asthmagenes_deseq2/SRP033351/SRP033351_de.feather") %>% tibble::as_tibble()
de[["SRP043162"]] <- read_feather("transcriptomics/asthmagenes_deseq2/SRP043162/SRP043162_de.feather") %>% tibble::as_tibble()
de[["SRP098649"]] <- read_feather("transcriptomics/asthmagenes_deseq2/SRP098649/SRP098649_de.feather") %>% tibble::as_tibble()
de[["SRP005411"]] <- read_feather("transcriptomics/asthmagenes_deseq2/SRP005411/SRP005411_de.feather") %>% tibble::as_tibble()

#Deseq2 count results - by gene for plots
tpms <- list()

tpms[["SRP033351"]]  <- read_feather("transcriptomics/asthmagenes_deseq2/SRP033351/SRP033351_pheno+counts_updated.feather") %>% tibble::as_tibble()

tpms[["SRP043162"]] <- read_feather("transcriptomics/asthmagenes_deseq2/SRP043162/SRP043162_pheno+counts_updated.feather") %>% tibble::as_tibble()

tpms[["SRP098649"]] <- read_feather("transcriptomics/asthmagenes_deseq2/SRP098649/SRP098649_pheno+counts_updated.feather") %>% tibble::as_tibble()

tpms[["SRP005411"]] <- read_feather("transcriptomics/asthmagenes_deseq2/SRP005411/SRP005411_pheno+counts.feather") %>% tibble::as_tibble()
 
# make a list of gene symbols in all datasets for checking whether gene symbol entered is valid - used later on
deseq2_filtered_genes <- unique(c(de$SRP005411$gene_symbol, de$SRP043162$gene_symbol, de$SRP033351$gene_symbol, de$SRP005411$gene_symbol))
