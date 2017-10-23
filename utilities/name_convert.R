###
# Usage
###

# This script is used for converting short names to full names. Make sure names in each vector are in correct order!
tissue_short_name <- c("ASM", "BE", "BAL", "CD4", "CD8", "LEC", "chALL", "LCL", "MACRO", "MCF10A-Myc", "NE", "U2O", "PBMC", "SAE","WBC", "Lung")
tissue_full_name <- c("Airway smooth muscle", "Bronchial epithelium", "Bronchoalveolar lavage", "CD4", "CD8", "Lens epithelium", "Lymphoblastic leukemia cell", "Lymphoblastoid cell", "Macrophage", "MCF10A-Myc", "Nasal epithelium", "Osteosarcoma U2OS cell", "Peripheral blood mononuclear cell","Small airway epithelium", "White blood cell","Whole lung")
asthma_short_name <- c("allergic_asthma", "asthma", "asthma_and_rhinitis","fatal_asthma", "mild_asthma", "non_allergic_asthma", "severe_asthma")
asthma_full_name <- c("Allergic asthma", "Asthma", "Asthma and rhinitis", "Fatal asthma", "Mild asthma", "Non-allergic asthma", "Severe asthma")
treatment_short_name <- c("BA", "GC", "smoking", "vitD")
treatment_full_name <- c("Beta-agonist treatment", "Glucocorticoid treatment", "Smoking", "Vitamin D treatment")

short_name <- c(tissue_short_name, asthma_short_name, treatment_short_name)
full_name <- c(tissue_full_name, asthma_full_name, treatment_full_name)
names(full_name) <- short_name

nameconvert <- function(name) {
	return(unname(full_name[name]))
}
