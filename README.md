## Reducing Associations by Linking Genes And transcriptomic Results (REALGAR)

REALGAR is a an app that integrates disease-specific, tissue-specific omics results into a tool for designing functional validation studies.  REALGAR uses asthma as a disease model to demonstrate how integrating omics results with GWAS data leads to improved understanding of gene associations.

Authors: Avantika R Diwadkar, Mengyuan Kan, Maya Shumyatcher, Blanca Himes.

### Rationale

Genetic association studies (genome-wide association studies, whole-genome sequencing studies, etc.), have uncovered many novel disease-associated variants, but few of these disease-associated regions have led to clinically actionable insights.  The gap between genetic associations and practical applications is due in part to the difficulty of designing functional validation studies that capture the complexity of the mechanisms in question.   

REALGAR is an integrated resource of disease-specific and tissue-specific results from omics studies which simplifies functional validation experiment design. Using asthma as a disease model, this app brings together omics data, including 1) genome-wide association study (GWAS) data from GRASP, GABRIEL, EVE, Ferreira et al, TAGC, and UK Biobank, 2) transcriptomic data (microarray and RNA-Seq) from the Gene Expression Omnibus (GEO), 3) epigenomic data (ChIP-Seq) of glucocorticoid receptor (GR)-binding sites and putative glucocorticoid response element (GRE) motifs derived from GEO and Sequence Read Archive (SRA), and 4) eQTL data of lung, whole blood and skeletal muscle tissues from Genotype-Tissue Expression (GTEx), allowing researchers to access a breadth of information with a click. REALGAR's disease-specific and tissue-specific information helps guide validation experiments for gene associations, with the aim of discovering clinically actionable insights.

### Dependencies

REALGAR was created using RStudio's Shiny.  To run REALGAR locally, R should be available.  The following packages should be installed: shiny, shinythemes, data.table, dplyr, DT, feather, ggplot2, karyoploteR, lattice, stringi, stringr, tidyverse, tidyr, TxDb.Hsapiens.UCSC.hg38.knownGene, and viridis.

### Citations

If you use REALGAR in your research, please cite the following papers:

Kan M, Diwadkar AR, Saxena S, Shuai H, Joo J, Himes BE. **REALGAR: a web app of integrated respiratory omics data**. *"Bioinformatics*. 2022 Sep 15;38(18):4442-4445. PMID: [35863045](https://www.ncbi.nlm.nih.gov/pubmed/35863045)

Shumyatcher M, Hong R, Levin J, Himes BE. **Disease-Specific Integration of Omics Data to Guide Functional Validation of Genetic Associations**. *AMIA Annu Symp Proc*. 2018;2017:1589–1596. 
PMID: [29854229](https://www.ncbi.nlm.nih.gov/pubmed/29854229)

Analysis of gene expression microarray and RNA-Seq data was performed as described in this paper:

Kan M, Shumyatcher M, Diwadkar A, Soliman G, Himes BE. **Integration of Transcriptomic Data Identifies Global and Cell-Specific Asthma-Related Gene Expression Signatures**. *AMIA Annu Symp Proc*. 2018;2018:1338–1347. 
PMID: [30815178](https://www.ncbi.nlm.nih.gov/pubmed/30815178
 
Analysis of ChIP-Seq data was performed as described in this paper:

Diwadkar AR, Kan M, Himes BE. **Facilitating Analysis of Publicly Available ChIP-Seq Data for Integrative Studies**. *AMIA Annu Symp Proc*. 2019;2019:371-379. PMID:[32308830](https://www.ncbi.nlm.nih.gov/pubmed/32308830)
