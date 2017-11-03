# detach("package:Gviz", unload=TRUE) # this is to keep RStudio happy - run if loading app more than once in same session - keep commented out otherwise
# if load Gviz 2x in same session (i.e. close & re-run app), get "object of type 'closure' is not subsettable" error
# should not be an issue when running app from the website
# cat(file=stderr(), as.character(Sys.time()),"packages start\n")
# use this type of command to easily see dataset loading time in RStudio  
# currently 3 seconds from "start package load" to "finish gene_locations load"
library(shiny)
library(dplyr)
library(data.table)
library(forestplot) 
library(lattice)
library(stringr)
library(viridis) 
library(DT) 
library(Gviz)
source("utilities/meta.R")
source("utilities/comb_pval.R")
source("utilities/name_convert.R")

# load dataset descriptions
Dataset_Info <- readRDS("databases/microarray_data_infosheet_R.RDS")

#load and name GEO microarray and RNA-Seq datasets
for (i in Dataset_Info$Unique_ID) {assign(i, readRDS(paste0("databases/microarray_results/", i, ".RDS")))} 

Dataset_Info$PMID <- as.character(Dataset_Info$PMID) #else next line does not work
Dataset_Info[is.na(Dataset_Info$PMID), "PMID"] <- ""

#load info for gene tracks: gene locations, TFBS, SNPs, etc.
tfbs <- readRDS("databases/tfbs_for_app.RDS") #TFBS data from ENCODE - matched to gene ids using bedtools
snp <- readRDS("databases/grasp_output_for_app.RDS") #SNP data from GRASP - matched to gene ids using bedtools
snp_eve <- readRDS("databases/eve_data_realgar.RDS") #SNP data from EVE - was already in hg19 - matched to gene ids using bedtools 
snp_gabriel <- readRDS("databases/gabriel_data_realgar.RDS") #SNP data from GABRIEL - lifted over from hg17 to hg19 - matched to gene ids using bedtools 
gene_locations <- fread("databases/gene_positions.txt", header = TRUE, stringsAsFactors = FALSE) #gene location & transcript data from GENCODE
chrom_bands <- readRDS("databases/chrom_bands.RDS") #chromosome band info for ideogram - makes ideogram load 25 seconds faster
#unlike all other files, gene_locations is faster with fread than with readRDS (2s load, vs 4s)

#compute -log10 for SNPs -- used for SNP colors
snp <- dplyr::mutate(snp, neg_log_p = -log10(p))
snp_eve <- dplyr::mutate(snp_eve, neg_log_meta_p = -log10(meta_P),
                         neg_log_meta_p_aa = -log10(meta_P_AA),
                         neg_log_meta_p_ea = -log10(meta_P_EA),
                         neg_log_meta_p_lat = -log10(meta_P_LAT))
snp_gabriel <- dplyr::mutate(snp_gabriel, neg_log_p = -log10(P_ran))

#color tfbs based on binding score - used in tracks
#create color scheme based on encode binding score & snp p-values
tfbs$color <- inferno(50)[as.numeric(cut(tfbs$score,breaks = 50))]

breaks <- c(seq(0,8,by=0.001), Inf) # this sets max universally at 8 (else highest one OF THE SUBSET would be the max)

snp$color <- inferno(8002)[as.numeric(cut(snp$neg_log_p, breaks = breaks))]
snp_gabriel$color <- inferno(8002)[as.numeric(cut(snp_gabriel$neg_log_p, breaks = breaks))]

snp_eve$color_meta_P <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p,breaks = breaks))]
snp_eve$color_meta_P_AA <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p_aa,breaks = breaks))]
snp_eve$color_meta_P_EA <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p_ea,breaks = breaks))]
snp_eve$color_meta_P_LAT <- inferno(8002)[as.numeric(cut(snp_eve$neg_log_meta_p_lat,breaks = breaks))]

# make a list of gene symbols in all datasets for checking whether gene symbol entered is valid - used for GeneSymbol later on
genes_avail <- vector()
for (i in ls()[grep("GSE", ls())]) {
    gene_names <- as.character(levels(get(i)$SYMBOL))
    genes_avail <- unique(c(genes_avail, gene_names))
}

output.table <- data.frame() # initiate output table - used later in output.tableforplot()
heatmap_colors <-  inferno # heatmap colors - used in p-value plot

# server
shinyServer(function(input, output, session) {
    
    curr_gene <- reactive({
        if (gsub(" ", "", tolower(input$curr_gene), fixed=TRUE) %in% c(snp$snp, snp_eve$snp, snp_gabriel$snp)) { #if SNP ID is entered, convert internally to nearest gene symbol  
            all_matches <- rbind(rbind(snp[which(snp$snp==gsub(" ", "", tolower(input$curr_gene), fixed=TRUE)), c("snp", "end", "symbol")], 
                                       snp_eve[which(snp_eve$snp==gsub(" ", "", tolower(input$curr_gene), fixed=TRUE)), c("snp", "end", "symbol")]), 
                                 snp_gabriel[which(snp_gabriel$snp==gsub(" ", "", tolower(input$curr_gene), fixed=TRUE)), c("snp", "end", "symbol")])
            gene_locations_unique <- gene_locations[which(!duplicated(gene_locations$symbol)),]
            all_matches <- merge(all_matches, gene_locations_unique[,c("symbol", "start")], by="symbol")
            all_matches$dist <- abs(all_matches$start - all_matches$end) # here, "end" is snp position, "start" is gene start 
            unique(all_matches$symbol[which(all_matches$dist==min(all_matches$dist))]) # choose the gene symbol whose start is the smallest absolute distance away
        } else { 
            # if it is not in the list of snps, it is a gene id OR a snp that is not associated with asthma
            # in the latter case it will not show up in the list of genes & user gets an "enter valid gene/snp id" message
            gsub(" ", "", toupper(input$curr_gene), fixed = TRUE) #make uppercase, remove spaces
        }
    })  
    
    GeneSymbol <- reactive({if (curr_gene() %in% genes_avail) {TRUE} else {FALSE}})  #used later to generate error message when a wrong gene symbol is input
    
    #####################################################################
    ## "Select all" buttons for tissue, asthma and treatment selection ##
    #####################################################################
    
    #Tissue
    tissue_choices <- c("Airway smooth muscle"="ASM", "Bronchial epithelium"="BE", 
                        "Bronchoalveolar lavage"="BAL", "CD4"="CD4", "CD8"="CD8",
                        "Lens epithelium" = "LEC","Lymphoblastic leukemia cell" = "chALL", 
                        "Lymphoblastoid cell" = "LCL","Macrophage" = "MACRO", "MCF10A-Myc" = "MCF10A-Myc",
                        "Nasal epithelium"="NE","Osteosarcoma U2OS cell" = "U2O", 
                        "Peripheral blood mononuclear cell"="PBMC","Small airway epithelium"="SAE",
                        "White blood cell"="WBC","Whole lung"="Lung")
    observe({
        if(input$selectall_tissue == 0) return(NULL) # don't do anything if action button has been clicked 0 times
        else if (input$selectall_tissue%%2 == 0) { # %% means "modulus" - i.e. here you're testing if button has been clicked a multiple of 2 times
            updateCheckboxGroupInput(session,"Tissue","Tissue",choices=tissue_choices, inline = TRUE)
            updateActionButton(session, "selectall_tissue", label="Select all") # change action button label based on user input
        } else { # else is 1, 3, 5 etc.
            updateCheckboxGroupInput(session,"Tissue","Tissue",choices=tissue_choices,selected=c("BE", "LEC", "NE", "CD4", "CD8", "PBMC", "WBC", "ASM", "BAL", "Lung",
                                                                                                 "chALL", "MCF10A-Myc", "MACRO", "U2O", "LCL", "SAE"), inline = TRUE)
            updateActionButton(session, "selectall_tissue", label="Unselect all")
        }
    })
    
    
    #Asthma
    asthma_choices <- c("Allergic asthma"="allergic_asthma", "Asthma"="asthma", "Asthma and rhinitis"="asthma_and_rhinitis",
                        "Fatal asthma"="fatal_asthma", "Mild asthma"="mild_asthma", "Non-allergic asthma"="non_allergic_asthma",
                        "Severe asthma"="severe_asthma")
    observe({
        if(input$selectall_asthma == 0) return(NULL) 
        else if (input$selectall_asthma%%2 == 0) {
            updateCheckboxGroupInput(session,"Asthma","Asthma",choices=asthma_choices)
            updateActionButton(session, "selectall_asthma", label="Select all")
        } else {
            updateCheckboxGroupInput(session,"Asthma","Asthma",
                                     choices=asthma_choices,
                                     selected=c("allergic_asthma", "asthma", "asthma_and_rhinitis", "fatal_asthma", "mild_asthma", "non_allergic_asthma", "severe_asthma"))
            updateActionButton(session, "selectall_asthma", label="Unselect all")
        }})
    
    
    #Treatment
    treatment_choices <- c("Beta-agonist treatment"="BA", "Glucocorticoid treatment" = "GC", "Smoking"="smoking", "Vitamin D treatment"="vitD")
    
    observe({
        if(input$selectall_treatment == 0) return(NULL) 
        else if (input$selectall_treatment%%2 == 0) {
            updateCheckboxGroupInput(session,"Treatment","Treatment",choices=treatment_choices)
            updateActionButton(session, "selectall_treatment", label="Select all")
        }
        else {
            updateCheckboxGroupInput(session,"Treatment","Treatment",choices=treatment_choices,selected=c("BA", "GC", "smoking","vitD"))
            updateActionButton(session, "selectall_treatment", label="Unselect all")
        }})
    
    #########################################
    ## reactive UI for EVE p-value options ##
    #########################################
    
    #if EVE SNPs selected, display option to choose population
    output$eve_options <- reactive({if("snp_eve_subs" %in% input$which_SNPs) {"GWAS display options:"} else {""}})
    
    #######################
    ## GEO studies table ##
    #######################
    #select GEO studies matching desired conditions;
    #Jessica's initial app had an "and" condition here; Maya changed it to "or"
    # Mengyuan changed it to: if only select tissue or asthma/treatment, will use all the available studies; else use the intersection
    UserDataset_Info <- reactive({
        #Dataset_Info1 = subset(Dataset_Info,(((Dataset_Info$Tissue %in% input$Tissue) | (Dataset_Info$Asthma %in% input$Asthma)) & Dataset_Info$App == "asthma")) 
        #Dataset_Info2 = subset(Dataset_Info,(((Dataset_Info$Tissue %in% input$Tissue) | (Dataset_Info$Asthma %in% input$Treatment)) & Dataset_Info$App == "GC")) 
        ## Dataset_Info2 = subset(Dataset_Info, (((Dataset_Info$Tissue %in% input$Tissue) & ((Dataset_Info$Asthma %in% input$Treatment) | (Dataset_Info$App %in% input$Treatment)))))
        #Dataset_Info = rbind(Dataset_Info1, Dataset_Info2) # this separates GC and asthma data
        Dataset_Info_Tissue = subset(Dataset_Info, Dataset_Info$Tissue %in% input$Tissue)
        Dataset_Info_Asthma = subset(Dataset_Info, Dataset_Info$Asthma %in% input$Asthma | Dataset_Info$Asthma %in% input$Treatment)
        if ((nrow(Dataset_Info_Tissue)==0)|(nrow(Dataset_Info_Asthma)==0)) {Dataset_Info1=rbind(Dataset_Info_Tissue,Dataset_Info_Asthma)}
        else {Dataset_Info1=subset(Dataset_Info_Tissue,Dataset_Info_Tissue$Unique_ID%in%Dataset_Info_Asthma$Unique_ID)}
        Dataset_Info1
    })
    
    # print out further available options based on the current selected items
    
    avail_text <- reactive({
        UserDataset_Info <- UserDataset_Info()
        avail_asthma <- NULL
        avail_GC <- NULL
        avail_tissue <- NULL
        convname_func <- function(x){
            sort(sapply(as.character(x),function(x){nameconvert(x)}))
        }
        if (is.null(input$Tissue)&is.null(input$Asthma)&is.null(input$Treatment)) {text="Please select at least one data type"}
        else if (is.null(input$Tissue)&(!(is.null(input$Asthma)&is.null(input$Treatment)))) {
            avail_tissue=unique(Dataset_Info$Tissue[Dataset_Info$Asthma %in% input$Asthma|Dataset_Info$Asthma %in% input$Treatment])
            avail_tissue_fullname <- convname_func(avail_tissue)
            text=paste0("Based on the asthma endotype(s) and/or treatment(s) selected, these tissue(s) are available: ", paste(avail_tissue_fullname,collapse=", "),".")
        }
        else if ((!is.null(input$Tissue))&(is.null(input$Asthma)|is.null(input$Treatment))) {
            avail_asthma=unique(Dataset_Info$Asthma[(Dataset_Info$Tissue %in% input$Tissue)&(Dataset_Info$App == "asthma")])
            avail_GC=unique(Dataset_Info$Asthma[(Dataset_Info$Tissue %in% input$Tissue)&(Dataset_Info$App == "GC")])
            if (length(avail_asthma)>0) {avail_asthma_fullname <- convname_func(avail_asthma)} else {avail_asthma_fullname<-NULL}
            if (length(avail_GC)>0) {avail_GC_fullname <- convname_func(avail_GC)} else {avail_GC_fullname <- NULL}
            text=paste0("Based on the tissue(s) selected, these asthma endotype(s) and/or treatment(s) are available: ", paste(c(avail_asthma_fullname,avail_GC_fullname),collapse=", "),".")
        } else { # specific tissue(s), asthma type(s) and treatment(s) have been selected
            text=""
        }
        text})
    output$avail_choice = renderText(avail_text())
    
    # disable the notavailable function
    #  notavail_text <- reactive({
    #      UserDataset_Info <- UserDataset_Info()
    #      if (!is.null(input$Tissue)&(!is.null(input$Asthma)|!is.null(input$Treatment))) {
    #          not_avail_tissue <- input$Tissue[!input$Tissue%in%UserDataset_Info$Tissue]
    #          asthma_GC <- c(input$Asthma, input$Treatment)
    #          not_avail_asthma_GC <-asthma_GC[!asthma_GC%in%UserDataset_Info$Asthma]
    #	  if (length(not_avail_tissue)>0) {
    #	      not_avail_tissue_fullname <- sapply(not_avail_tissue,function(x){nameconvert(x)})
    #	      text_tissue <- paste0("These tissue(s) are not available in selected asthma type(s) and/or treatment(s): ",paste(not_avail_tissue_fullname,collapse=", "),". ")
    #	  } else {text_tissue <- ""}
    #	  if (length(not_avail_asthma_GC)>0) {
    #	      not_avail_asthma_GC_fullname <- sapply(not_avail_asthma_GC,function(x){nameconvert(x)})
    #	      text_asthma_GC <- paste0("These asthma type(s) and/or treatment(s) are not available in selected tissue(s): ",paste(not_avail_asthma_GC_fullname,collapse=", "),".")
    #	  } else {text_asthma_GC <- ""}
    #	  text=paste0(text_tissue,text_asthma_GC)
    #      } else {text=""}
    #      text})
    notavail_text <- reactive({text=""})
    output$notavail_choice = renderText(notavail_text())
    
    
    #add links for GEO_ID and PMID
    GEO_data <- reactive({
        validate(need(nrow(UserDataset_Info()) != 0, "Please choose at least one dataset.")) #Generate a error message when no data is loaded.
        
        UserDataset_Info() %>%
            dplyr::mutate(GEO_ID_link = ifelse(grepl("SRP", GEO_ID), #GEO link is conditional on whether GEO_ID is an "SRP" or "GSE"
                                               paste0("http://www.ncbi.nlm.nih.gov/sra/?term=", GEO_ID), 
                                               paste0("http://www.ncbi.nlm.nih.gov/gquery/?term=", GEO_ID)),
                          PMID_link = paste0("http://www.ncbi.nlm.nih.gov/pubmed/?term=", PMID))})
    
    
    
    Dataset <- reactive({paste0("<a href='",  GEO_data()$GEO_ID_link, "' target='_blank'>",GEO_data()$GEO_ID,"</a>")})
    PMID <- reactive({paste0("<a href='",  GEO_data()$PMID_link, "' target='_blank'>",GEO_data()$PMID,"</a>")})
    Description <- reactive({GEO_data()$Description})
    
    GEO_links <- reactive({
        df <- data.frame(Dataset(), PMID(), Description())
        colnames(df) <- c("Dataset", "PMID", "Description")
        df})
    
    output$GEO_table <- DT::renderDataTable(GEO_links(),  
                                            class = 'cell-border stripe', 
                                            rownames = FALSE, 
                                            options = list(paging = FALSE, searching = FALSE),
                                            escape=FALSE)
    #########################################
    ## Select GEO data for plots and table ##
    #########################################
    
    #select and modify data used for plots and accompanying table
    output.tableforplot <- reactive({
        validate(need(nrow(UserDataset_Info()) != 0, "Please choose at least one dataset.")) #Generate a error message when no data is loaded.
        validate(need(curr_gene() != "", "Please enter a gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
        
        #select data for the gene currently selected
        data_filter <- function(x){
            x <- x %>%
                dplyr::filter(SYMBOL==curr_gene()) %>%
                dplyr::select(logFC, P.Value, adj.P.Val, SD) %>% 
                dplyr::filter(P.Value==min(P.Value)) %>%
                dplyr::mutate(lower = logFC - 2*SD, upper = logFC + 2*SD)
        }
        
        #get data for given gene for each study selected
        for (i in UserDataset_Info()$Unique_ID){
            curr.gene.data=get(i,environment())
            Total = UserDataset_Info() %>% dplyr::filter(Unique_ID == i) %>% select(Total) # This 'Total' is total sample size and can be used to combine p-values
            data_type = UserDataset_Info() %>% dplyr::filter(Unique_ID == i) %>% select(App) #This 'data_type' can be used to separate asthma and GC data. 
            
            if(any(tbl_vars(curr.gene.data)=="qValuesBatch")) {
                curr.gene.data <- (curr.gene.data %>%
                                       dplyr::select(-P.Value,-adj.P.Val) %>%
                                       dplyr::rename(P.Value=pValuesBatch) %>%
                                       dplyr::rename(adj.P.Val=qValuesBatch))}
            
            #use data_filter function from above to filter curr.gene.data
            if (any(GeneSymbol())) {
                
                curr.gene.data <- data_filter(curr.gene.data)
                
                if(nrow(curr.gene.data) > 0) {
                    curr.gene.data <- curr.gene.data[order(curr.gene.data$P.Value,-abs(curr.gene.data$logFC)),][1,] # if multiple probes have the same smallest p-values, select the one with largest effect
                    curr.gene.data <- cbind(data_type, Unique_ID=i, curr.gene.data, Total)
                    #append curr.gene.data to all the other data that needs to be included in the levelplots
                    output.table <- rbind(output.table, curr.gene.data)}}}
        
        #preparing the data for levelplots
        #calculate the fold change, order by fold change for levelplots
        validate(need(GeneSymbol() != FALSE, "Please enter a valid gene symbol or SNP ID.")) # Generate error message if the gene symbol is not right.
        output.table <- dplyr::mutate(output.table, Fold_Change=2^(logFC), neglogofP=(-log10(adj.P.Val)), Lower_bound_CI = 2^(lower), Upper_bound_CI = 2^(upper)) #note that this is taking -log10 of adjusted p-value
        # row.names(output.table) <- output.table$Unique_ID #crucial for plot labels on levelplot
        output.table <- output.table[order(output.table$Fold_Change, output.table$Upper_bound_CI),]
        #print(output.table)
    })
    
    
    ###################################
    ## Data table for meta-analysis  ##
    ###################################
    
    # asthma
    data_Asthma <- reactive({ output.tableforplot_asthma = output.tableforplot() 
    output.tableforplot_asthma = output.tableforplot_asthma[output.tableforplot_asthma$App == "asthma",]
    output.tableforplot_asthma[rev(rownames(output.tableforplot_asthma)),]})
    
    # GC
    data_GC <- reactive({ output.tableforplot_GC = output.tableforplot()
    output.tableforplot_GC = output.tableforplot_GC[output.tableforplot_GC$App %in% c("GC", "BA", "smoking", "vitD"),]
    output.tableforplot_GC[rev(rownames(output.tableforplot_GC)),]})
    ###################################
    ##        Combined p-values      ##
    ###################################
    
    # asthma
    asthma_pcomb <- reactive({
        dat <- data_Asthma()
        if (length(dat$adj.P.Val)>1) {
            #write.table(dat,paste0("Asthma_",curr_gene(),".txt"),col.names=T,row.names=F,sep="\t",quote=F)
            asthma_liptak_pcomb <- liptak_stat(dat)
            asthma_sumlog_pcomb <- sumlog_stat(dat)
            pcomb_text=paste0("Combined p-values by Liptak's method =", asthma_liptak_pcomb, "; by Fisher's method = ", asthma_sumlog_pcomb)
        }
        else {pcomb_text=""}
        pcomb_text
    })
    
    output$asthma_pcomb_text <- renderText({asthma_pcomb()})
    
    # GC
    GC_pcomb <- reactive({
        dat <- data_GC()
        if (length(dat$adj.P.Val)>1) {
            #write.table(dat,paste0("GC_",curr_gene(),".txt"),col.names=T,row.names=F,sep="\t",quote=F)
            GC_liptak_pcomb <- liptak_stat(dat)
            GC_sumlog_pcomb <- sumlog_stat(dat)
            pcomb_text=paste0("Combined p-values by Liptak's method =", GC_liptak_pcomb, "; by Fisher's method = ", GC_sumlog_pcomb)
        }
        else {pcomb_text=""}
        pcomb_text
    })
    
    output$GC_pcomb_text <- renderText({GC_pcomb()})
    ###################################
    ##          Meta-analysis        ##
    ###################################
    
    # p-values, effect size and 95% CI will be used for forestplot
    # asthma
    meta_Asthma <- reactive({
        dat <- data_Asthma()
        if (length(dat$adj.P.Val)>1) {
            res <- meta_stat(dat)
        }
        else {res <- list(meta_pval=1,meta_fc=1,meta_lower=1,meta_upper=1)} # assign any values to trick the app if no "treatment" is selected
        res
    })  
    
    meta_GC <- reactive({
        dat <- data_GC()
        if (length(dat$adj.P.Val)>1) {
            res <- meta_stat(dat)
        }
        else {res <- list(meta_pval=NA,meta_fc=NA,meta_lower=NA,meta_upper=NA)}
        res
    })  
    
    ###################################
    ## Data table accompanying plots ##
    ###################################
    
    # asthma
    
    data2_Asthma <- reactive({ # dataset without meta-analysis results
        data_Asthma()%>%
            dplyr::select(Unique_ID, adj.P.Val, P.Value,Fold_Change, neglogofP, Lower_bound_CI, Upper_bound_CI) %>%
	    arrange(desc(Fold_Change),desc(Upper_bound_CI)) %>% # sort by first effect size (fold change) and then by upper CI in a descending order
            dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), P.Value =format(P.Value, scientific=TRUE, digits=3), 
                          Lower_bound_CI = round(Lower_bound_CI, digits = 2), Upper_bound_CI = round(Upper_bound_CI, digits = 2), Comparison = "Asthma vs. non-asthma")%>%
            dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Fold Change`=Fold_Change)})

    # Function: "interdata_func"
    # obtain intermediate data with meta-analysis results added for table output and forest plots
    interdata_func <- function(dat,meta_dat){
    	if (nrow(dat)==0) {return(data.frame())} # define an empty dataset if nothing is selected
	# add meta-analysis results if more than one study were observed
        if (nrow(dat)>1) {
            meta_pval=format(meta_dat[["meta_pval"]],scientific=TRUE, digits=3)
            meta_fc=round(meta_dat[["meta_fc"]],2)
            meta_lower=round(meta_dat[["meta_lower"]],2)
            meta_upper=round(meta_dat[["meta_upper"]],2)
            meta_neglogofP=-log10(meta_dat[["meta_pval"]])
            
            dat <- rbind(dat,rep(NA,ncol(dat)))
            dat$`Q Value`[nrow(dat)] <- meta_pval
            dat$`Fold Change`[nrow(dat)] <- meta_fc
            dat$Lower_bound_CI[nrow(dat)] <- meta_lower
            dat$Upper_bound_CI[nrow(dat)] <- meta_upper
            dat$neglogofP[nrow(dat)] <- meta_neglogofP
        }
	
	# merge information from info sheet
        text_temp <- merge(dat, as.matrix(Dataset_Info[which(Dataset_Info$Unique_ID %in% dat$`Study ID`),]), by.x="Study ID", by.y="Unique_ID", all.x=T)
        
	# sort by effect size (fold change) by individual study as the order changes after merging
	if (nrow(text_temp)>1) {
	    meta_row <- text_temp[nrow(text_temp),]
	    study_row <- text_temp[1:(nrow(text_temp)-1),] # sort exclude the meta-analysis result
	    study_row <- study_row[order(-study_row$`Fold Change`, -study_row$Upper_bound_CI),]
	    text_temp <- rbind(study_row,meta_row)
	}
	# select columns used for table output and forest plots
	text_temp <- text_temp[,c(names(dat),"GEO_ID","Asthma","Treatment","Long_tissue_name","Total","App")] # include total sample size
	# modify asthma endotypes or treatment conditions
	text_temp[,"Asthma"] <- gsub("_", " ", text_temp[,"Asthma"])
	text_temp[,"Treatment"] <- gsub("_", " ", text_temp[,"Treatment"])
	# compute total sample size for meta-analysis results
	text_temp[,"Total"] <- as.character(text_temp[,"Total"])
	text_temp[nrow(text_temp),"Total"] <- sum(as.numeric(text_temp[-nrow(text_temp),"Total"]))
	return(text_temp)
    }
  
    # Function: "tableforgraph_func" 
    # Obtain table output
    tableforgraph_func <- function(dat){
	empt_dat <- data.frame(matrix(ncol=6, nrow = 0)) # create an empty dataset if nothing is selected
	names(empt_dat) <- c("Study ID", "Tissue", "Comparison", "P Value", "Q Value", "Fold Change(95% CI)")
   	if (nrow(dat)==0) {return(empt_dat)}
	data4table <- dat %>% 
            dplyr::mutate(`Study ID`=GEO_ID, `Tissue`=Long_tissue_name, `Fold Change(95% CI)` = paste(`Fold Change`, " ","(", Lower_bound_CI, ",", Upper_bound_CI, ")", sep = "")) %>%
            dplyr::select(`Study ID`, `Tissue`, `Comparison`, `P Value`, `Q Value`, `Fold Change(95% CI)`) %>%
	    mutate_all(as.character)
    	if (nrow(data4table)>1) {
	    data4table[nrow(data4table),c("Study ID")] <- "Combined"
	    data4table[nrow(data4table),c("Comparison")] <- data4table[1,c("Comparison")]
	}
        return(data4table)
    }

    # obtain intermediate asthma data for table output and forest plots
    data3_Asthma <- reactive({
        interdata_func(data2_Asthma(),meta_Asthma())
    })

    # output table for asthma
    tableforgraph_Asthma <- reactive({
	tableforgraph_func(data3_Asthma())
    })
    
    # GC
    data2_GC <- reactive({
        data_GC()%>%
            dplyr::select(Unique_ID, adj.P.Val, P.Value,Fold_Change, neglogofP, Lower_bound_CI, Upper_bound_CI) %>%
            dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), P.Value =format(P.Value, scientific=TRUE, digits=3), 
                          Lower_bound_CI = round(Lower_bound_CI, digits = 2), Upper_bound_CI = round(Upper_bound_CI, digits = 2), Comparison = "Stimulation vs. at baseline")%>%
            dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Fold Change`=Fold_Change)})

    # obtain intermediate GC data for table output and forest plots
    data3_GC <- reactive({
        interdata_func(data2_GC(),meta_GC())
    })

    # output table for GC
    tableforgraph_GC <- reactive({
	tableforgraph_func(data3_GC())
    })


    #combine asthma & GC into one
    output$tableforgraph <- DT::renderDataTable(rbind(tableforgraph_Asthma(), tableforgraph_GC()),
                                                class = 'cell-border stripe', 
                                                rownames = FALSE, 
                                                options = list(paging = FALSE, searching = FALSE),
                                                width = "100%")
    
    
    #################
    ## Forestplots ##
    #################
    
    # Function: "forestplot_func" 
    # Forest plots

    forestplot_func <- function(dat) {
	
        validate(need(nrow(dat) != 0, "Please choose a dataset."))

	# select columns for forest plot text
        if ("asthma"%in%dat$App) {
            text_temp <- dat[,c("GEO_ID","Long_tissue_name","Asthma","Q Value")]
	} else {
	    text_temp <- dat[,c("GEO_ID","Long_tissue_name","Treatment","Q Value")]
	}
	# create an empty dataset for forestplot text
        tabletext <- data.frame(matrix(nrow=1, ncol=4))
	# assign column names same as the original data
	names(tabletext) <- names(text_temp)
	# create first line for forestplot text
	if ("asthma"%in%dat$App) {
            tabletext[1,] <- c("GEO ID", "Tissue", "Endotype", "Q Value")
	} else {
	    tabletext[1,] <- c("GEO ID", "Tissue", "Treatment", "Q Value")
	}
	# add text for individual study result
        tabletext <- rbind(tabletext,text_temp)

	# convert all columns into character
	tabletext <- tabletext %>% mutate_all(as.character)

	# assign variables to meta-analysis result row
        if (nrow(dat)>1) {
	    tabletext[nrow(tabletext),c("GEO_ID")] <- "Combined"
	    tabletext[nrow(tabletext),c("Long_tissue_name")] <- ""
            tabletext[nrow(tabletext),3] <- "" # "Asthma" or "Treatment" is in column 3
        }

	# remove double quote
        options(useFancyQuotes = FALSE)
        tabletext <- gsub('"', '', sapply(tabletext, dQuote))
	
	# table with fold changes for plot
        tableplot <- rbind(c(NA,NA,NA,NA),dat[,c("Fold Change","Lower_bound_CI","Upper_bound_CI")])	
        
        # function to color forestplot lines and boxes by -log10 of adjusted pvalue - always relative to the max of 8
        color_fn <- local({
            i <- 0
            breaks <- c(seq(0,8,by=0.001), Inf) # this sets max universally at 8 (else highest one OF THE SUBSET would be the max)
            b_clrs <- l_clrs <- inferno(8002)[as.numeric(cut(dat$neglogofP, breaks = breaks))] #8002 is length(breaks) - ensures there are enough colors
            
            function(..., clr.line, clr.marker){
                i <<- i + 1
                fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
            }
        })
        
        xticks = seq(from = min(0.9, min(dat$Lower_bound_CI)), to = max(max(dat$Upper_bound_CI),1.2), length.out = 5)
        
        #hrzl_lines are borders between rows... made wide enough to be a background
        size_par <- max(6, nrow(dat)) #else plot scaling messed up when fewer than 5 datasets selected
        hrzl_lines <- vector("list", nrow(tabletext)+1)
        # hrzl_lines[[2]] <- gpar(lwd=350/size_par, lineend="butt", columns=5, col="#BDB6B0")
        for (i in setdiff(c(3:(length(hrzl_lines)-1)),c(1))) {hrzl_lines[[i]]  <- if(nrow(dat)==1) {
            gpar(lwd=100, lineend="butt", columns=5, col="#BDB6B0")
        } else {
            gpar(lwd=1000/size_par, lineend="butt", columns=5, col="#BDB6B0")}
        }
        hrzl_lines[[1]] <- gpar(lwd=1000/size_par, lineend="butt", columns=5, col="#ffffff")
        hrzl_lines[[length(hrzl_lines)]] <- gpar(lwd=150/size_par, lineend="butt", columns=5, col="#ffffff")
        
        # adjust dot size based on sample size
        if (nrow(dat)>1) {
            total <- as.numeric(dat$Total)
            boxsize=c(0,0.1*log10(total)) # 0 for the header line
        } else {boxsize=0.2} # boxsize = 0.2 default

        if ("asthma"%in%dat$App) {
            title <- "Asthma vs. Non-asthma"
	} else {
	    title <- "Treatment vs. Control"
	}

        forestplot(tabletext, title = title, tableplot, zero = 1, 
                   xlab = "Fold Change", boxsize = boxsize, col = fpColors(zero="black"), 
                   lwd.ci = 2, xticks = xticks, is.summary=c(TRUE,rep(FALSE,nrow(dat))), 
                   lineheight = unit(19.7/size_par, "cm"), mar = unit(c(5,0,0,5),"mm"), fn.ci_norm = color_fn,
                   txt_gp = fpTxtGp(cex = 1.2, xlab = gpar(cex = 1.35), ticks = gpar(cex = 1.2), title = gpar(cex = 1.45)),
                   hrzl_lines=hrzl_lines)
    }
    
    # asthma forestplot

    forestplot_asthma <- reactive({
	forestplot_func(data3_Asthma())
    })

    # GC forestplot

    forestplot_GC <- reactive({
	forestplot_func(data3_GC())
    })
    
    output$forestplot_asthma = renderPlot({forestplot_asthma()}, height=650)
    output$forestplot_GC = renderPlot({forestplot_GC()}, height=650)
    
    output$color_scale1 <- output$color_scale2 <- renderImage({ #need two separate output names - else it fails (can't output same thing twice?)
        return(list(
            src = "databases/www/color_scale_vertical.png",
            height=550,
            width=59,
            filetype = "image/png",
            alt = "color_scale"))}, deleteFile = FALSE)
    
    
    ###############################
    ## Gene, SNP and TFBS tracks ##
    ###############################
    
    #horizontal color scale for gene tracks
    output$color_scale3 <- renderImage({ 
        return(list(
            src = "databases/www/color_scale_horizontal.png",
            height=109*1.05,
            width=1015*1.05,
            filetype = "image/png",
            alt = "color_scale"))}, deleteFile = FALSE)
    
    #filter data for selected gene
    gene_subs <- reactive({
        gene_subs_temp <- unique(filter(gene_locations, symbol==curr_gene()))
        gene_subs_temp <- gene_subs_temp[!(duplicated(gene_subs_temp$exon)),]})
    tfbs_subs <- reactive({unique(filter(tfbs, symbol==curr_gene()))})
    snp_subs <- reactive({
        if(("snp_subs" %in% input$which_SNPs)) {unique(filter(snp, symbol==curr_gene()))} 
        else {data.frame(matrix(nrow = 0, ncol = 0))}}) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
    snp_eve_subs <- reactive({
        if(("snp_eve_subs" %in% input$which_SNPs)) {
            snp_eve_temp <- snp_eve[which((snp_eve$symbol==curr_gene()) & (!is.na(unlist(snp_eve[,paste0("color_", input$which_eve_pvals)])))),] 
            #need the second filter criterion because otherwise will output snp names & otherwise blank if NA pvalues
            if (nrow(snp_eve_temp) > 0) {snp_eve_temp} else {data.frame(matrix(nrow = 0, ncol = 0))} #since pval_selector might remove all rows 
        } else {data.frame(matrix(nrow = 0, ncol = 0))}
    }) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
    
    snp_gabriel_subs <- reactive({
        if(("snp_gabriel_subs" %in% input$which_SNPs)) {unique(filter(snp_gabriel, symbol==curr_gene()))}
        else {data.frame(matrix(nrow = 0, ncol = 0))}}) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
    
    gene_tracks <- function() {
        validate(need(curr_gene() != "", "Please enter a gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
        validate(need(GeneSymbol() != FALSE, "Please enter a valid gene symbol or SNP ID.")) # Generate error message if the gene symbol is not right.
        validate(need(nrow(UserDataset_Info()) != 0, "Please choose at least one dataset.")) #Generate a error message when no data is loaded.
        
        gene_subs <- gene_subs()
        tfbs_subs <- tfbs_subs()
        snp_subs <- snp_subs()
        snp_eve_subs <- snp_eve_subs()
        snp_gabriel_subs <- snp_gabriel_subs()
        
        #for better visibility, increase tfbs and snp widths -- need scaling factor b/c different genes take up different amounts of space
        smallest_start <- min(gene_subs$start, tfbs_subs$start, snp_subs$start, snp_eve_subs$start, snp_gabriel_subs$start)
        largest_end <- max(gene_subs$end, tfbs_subs$end, snp_subs$end, snp_eve_subs$end, snp_gabriel_subs$end)
        scaling_factor <- (largest_end - smallest_start)/120000
        
        tfbs_subs$end <- tfbs_subs$end + 500*scaling_factor
        snp_subs$end <- snp_subs$end + 300*scaling_factor
        snp_eve_subs$end <- snp_eve_subs$end + 300*scaling_factor
        snp_gabriel_subs$end <- snp_gabriel_subs$end + 300*scaling_factor
        
        #constant for all tracks
        gen <- "hg19"
        chr <- unique(gene_subs$chromosome)
        
        #chromosome, axis and gene - these tracks show up for all genes
        #note that "col" refers to outline colors, whereas "fill" is the fill color
        bands <- chrom_bands[which(chrom_bands$chrom==chr),]
        chrom_track <- IdeogramTrack(genome = gen, bands = bands, fontcolor="black", fontsize=16) # formerly slow b/c of chromosome=chr; see https://support.bioconductor.org/p/78881/
        axis_track <- GenomeAxisTrack(col="black", fontcolor="black", fontsize=16)
        gene_track <- Gviz::GeneRegionTrack(gene_subs, genome = gen, chromosome = chr, name = "Transcripts", transcriptAnnotation="transcript", fill = "royalblue", col=NULL, grid=TRUE, col.grid="darkgrey") #add stacking="dense" if want transcript stacks combined into one
        
        #tfbs and snp tracks - only present for some genes
        
        #TFBS 
        if (nrow(tfbs_subs) > 0) {tfbs_track <- Gviz::AnnotationTrack(tfbs_subs, name="GR binding", fill = tfbs_subs$color, col=NULL, feature = tfbs_subs$score, grid=TRUE, col.grid="darkgrey")}
        
        # GRASP SNPs track
        if (nrow(snp_subs) > 0) { 
            snp_track <- Gviz::AnnotationTrack(snp_subs, name="SNPs (GRASP)", fill = snp_subs$color, col=NULL, feature=snp_subs$snp, grid=TRUE, col.grid="darkgrey")
            
            #rough estimate of number of stacks there will be in SNP track - for track scaling
            #note this stuff needs the SNPs to be ordered by position (smallest to largest)
            if (nrow(snp_subs) > 1) {
                snp_subs_temp <- snp_subs
                snp_range <- max(snp_subs_temp$start) - min(snp_subs_temp$start)
                snp_subs_temp$start_prev <- c(0, snp_subs_temp$start[1:(nrow(snp_subs_temp)-1)])
                snp_subs_temp$dist <- as.numeric(snp_subs_temp$start) - as.numeric(snp_subs_temp$start_prev)
                snp_size_init <- 0.9 + 3*as.numeric(nrow(snp_subs[which(snp_subs$dist < snp_range/10),])) + 0.3*length(unique(gene_subs$transcript))
            } else {snp_size_init <- 0.9 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs)}
        } else {snp_size_init <- 0.9 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs)}
        
        # EVE SNPs track
        if (nrow(snp_eve_subs) > 0) {
            pval_choice <- reactive({input$which_eve_pvals})  #pval_choice is responsible for dynamically coloring snps based on user selection of population
            snp_eve_track <- Gviz::AnnotationTrack(snp_eve_subs, name="SNPs (EVE)", fill = unlist(snp_eve_subs[, paste0("color_", pval_choice())]), col=NULL, feature=snp_eve_subs$snp, grid=TRUE, col.grid="darkgrey")
            
            #rough estimate of number of stacks there will be in SNP track - for track scaling
            if (nrow(snp_eve_subs) > 1) {
                snp_eve_subs_temp <- snp_eve_subs
                snp_eve_range <- max(snp_eve_subs_temp$start) - min(snp_eve_subs_temp$start)
                snp_eve_subs_temp$start_prev <- c(0, snp_eve_subs_temp$start[1:(nrow(snp_eve_subs_temp)-1)])
                snp_eve_subs_temp$dist <- as.numeric(snp_eve_subs_temp$start) - as.numeric(snp_eve_subs_temp$start_prev)
                snp_eve_size_init <- 1.5 + as.numeric(nrow(snp_eve_subs[which(snp_eve_subs$dist < snp_eve_range/10),])) + 0.3*length(unique(gene_subs$transcript))
            } else {snp_eve_size_init <- 1.4 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_eve_subs)}
        } else {snp_eve_size_init <- 1.4 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_eve_subs)}
        
        # GABRIEL SNPs track
        if (nrow(snp_gabriel_subs) > 0) {
            snp_gabriel_track <- Gviz::AnnotationTrack(snp_gabriel_subs, name="SNPs (GABRIEL)", fill = snp_gabriel_subs$color, col=NULL, feature=snp_gabriel_subs$snp, grid=TRUE, col.grid="darkgrey")
            
            #rough estimate of number of stacks there will be in SNP track - for track scaling
            if (nrow(snp_gabriel_subs) > 1) {
                snp_gabriel_subs_temp <- snp_gabriel_subs
                snp_gabriel_range <- max(snp_gabriel_subs_temp$start) - min(snp_gabriel_subs_temp$start)
                snp_gabriel_subs_temp$start_prev <- c(0, snp_gabriel_subs_temp$start[1:(nrow(snp_gabriel_subs_temp)-1)])
                snp_gabriel_subs_temp$dist <- as.numeric(snp_gabriel_subs_temp$start) - as.numeric(snp_gabriel_subs_temp$start_prev)
                snp_gabriel_size_init <- 1 + as.numeric(nrow(snp_gabriel_subs[which(snp_gabriel_subs$dist < snp_gabriel_range/10),])/4) + 0.12*length(unique(gene_subs$transcript))
            } else {snp_gabriel_size_init <- 1.4 + 0.1*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_gabriel_subs)}
        } else {snp_gabriel_size_init <- 1.4 + 0.1*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_gabriel_subs)}
        
        
        #track sizes - defaults throw off scaling as more tracks are added
        chrom_size <- 1.2 + 0.01*length(unique(gene_subs$transcript)) + 0.01*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.01*nrow(snp_gabriel_subs)
        axis_size <- 1 + 0.05*length(unique(gene_subs$transcript)) + 0.01*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.01*nrow(snp_gabriel_subs)
        gene_size <- 2 + 0.6*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs) + 0.05*nrow(snp_eve_subs) + 0.015*nrow(snp_gabriel_subs)
        tfbs_size <- 2 + 0.075*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs) + 0.05*nrow(snp_eve_subs) + 0.015*nrow(snp_gabriel_subs)
        snp_size <- snp_size_init #from above
        snp_eve_size <- snp_eve_size_init #from above
        snp_gabriel_size <- snp_gabriel_size_init #from above
        
        # #track sizes - defaults throw off scaling as more tracks are added --- if transcript stacking="dense"
        # chrom_size <- 1 + 0.005*length(unique(gene_subs$transcript)) + 0.005*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.005*nrow(snp_gabriel_subs)
        # axis_size <- 1.2+ 0.005*length(unique(gene_subs$transcript)) + 0.005*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.005*nrow(snp_gabriel_subs)
        # gene_size <- 1.75 + 0.005*length(unique(gene_subs$transcript)) + 0.005*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.005*nrow(snp_gabriel_subs)
        # tfbs_size <- 5 + 0.01*length(unique(gene_subs$transcript)) + 0.01*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.01*nrow(snp_gabriel_subs)
        # snp_size <- snp_size_init-1 #from above
        # snp_eve_size <- snp_eve_size_init-3 #from above
        # snp_gabriel_size <- snp_gabriel_size_init #from above
        
        #select the non-empty tracks to output -- output depends on whether there are TFBS and/or SNPs for a given gene
        subset_size <- sapply(c("tfbs_subs", "snp_subs", "snp_eve_subs", "snp_gabriel_subs"), function(x) {nrow(get(x))}) #size of each subset
        non_zeros <- names(subset_size)[which(!(subset_size==0))] #which subsets have non-zero size
        
        df_extract <- function(x,y) { #gives name of track and track size variable for non-zero subsets (y is "track" or "size")
            if (length(non_zeros) > 0) {
                get(paste0(strsplit(x, 'subs'),y)) #trim off "subs" and append either "track" or "size"
            } else {NULL} #to avoid meltdown if no subsets were non-zero
        }
        
        #use df_extract function to get track & track size corresponding to all non-zero subsets
        #note chrom_track, axis_track and gene_track are present for all
        selected_tracks <- list(chrom_track, axis_track, gene_track, sapply(non_zeros, df_extract, y="track")$tfbs_subs, sapply(non_zeros, df_extract, y="track")$snp_subs, sapply(non_zeros, df_extract, y="track")$snp_eve_subs, sapply(non_zeros, df_extract, y="track")$snp_gabriel_subs)
        selected_tracks <- Filter(Negate(function(x) is.null(unlist(x))), selected_tracks) #remove null elements from list
        
        selected_sizes <- na.omit(c(chrom_size,axis_size,gene_size, sapply(non_zeros, df_extract, y="size")[1], sapply(non_zeros, df_extract, y="size")[2], sapply(non_zeros, df_extract, y="size")[3], sapply(non_zeros, df_extract, y="size")[4]))
        selected_sizes <- Filter(Negate(function(x) is.null(unlist(x))), selected_sizes) #remove null elements from list - else run into trouble in conditions when no TFBS & no SNP tracks selected
        #note: use names to extract from selected_tracks b/c it is a list vs. index to extract from selected_sizes, since this is numeric
        
        #plot tracks 
        plotTracks(selected_tracks, sizes=selected_sizes, background.panel = "#BDB6B0", background.title = "firebrick4", col.border.title = "firebrick4", groupAnnotation = "feature", fontcolor.group = "darkblue", cex.group=0.75, just.group="below", cex.title=1.1)
        
    }
    
    #plot height increases if more tracks are displayed
    observe({output$gene_tracks_outp2 <- renderPlot({gene_tracks()}, width=1055,
                                                    height=400 + 15*length(unique(gene_subs()$transcript)) + 30*nrow(snp_eve_subs()) + 10*(nrow(snp_subs())+nrow(snp_gabriel_subs())))})
    
    
    #################################
    ## SNP data table for download ##
    #################################
    snp_subs_temp <- reactive({
        if (nrow(snp_subs()) > 0) {
            snp_subs()%>%
                dplyr::rename(position=start, meta_P=p) %>%
                dplyr::mutate(source = "GRASP") %>%
                dplyr::select(chromosome, snp, symbol, position, pmid, source, meta_P)
        }
    })
    
    snp_eve_subs_temp <- reactive({
        if (nrow(snp_eve_subs()) > 0) {
            snp_eve_subs() %>%
                dplyr::rename(position=start) %>%
                dplyr::mutate(source = "EVE") %>%
                dplyr::select(-c(end, color_meta_P, color_meta_P_EA, color_meta_P_AA, color_meta_P_LAT))
        }
    })
    
    snp_gabriel_subs_temp <- reactive({
        if (nrow(snp_gabriel_subs()) > 0) {
            snp_gabriel_subs() %>%
                dplyr::rename(position=start, meta_P=P_ran) %>%
                dplyr::mutate(source = "GABRIEL") %>%
                dplyr::select(-c(end, color))
        }
    })
    
    
    snp_data <- reactive({dplyr::bind_rows(snp_subs_temp(), snp_eve_subs_temp(), snp_gabriel_subs_temp())})
    
    ######################
    ## Download buttons ##
    ######################
    graphgene=reactive({curr_gene()})
    
    output$asthma_fc_download <- downloadHandler(
        filename= function(){paste0("fold_change_asthma_", graphgene(), "_", Sys.Date(), ".png")},
        content=function(file){
            png(file, width=12, height=9, units="in", res=300)
            forestplot_asthma()
            dev.off()})
    
    output$GC_fc_download <- downloadHandler(
        filename= function(){paste0("fold_change_GC_", graphgene(), "_", Sys.Date(), ".png")},
        content=function(file){
            png(file, width=12, height=9, units="in", res=300)
            forestplot_GC()
            dev.off()})
    
    # output$pval_download <- downloadHandler(
    #     filename= function(){paste0("-log(pval)_heatmap_", graphgene(), "_", Sys.Date(), ".png")},
    #     content=function(file){
    #         png(file, width=6, height=9, units="in", res=600)
    #         print(pval_plot()) # note that for this one, unlike other plot downloads, had to use print(). 
    #         dev.off()})        # else the download is a blank file. this seems to be b/c pval_plot() creates a graph 
    #                            # object but doesn't draw the plot, as per 
    #                            # http://stackoverflow.com/questions/27008434/downloading-png-from-shiny-r-pt-2
    
    output$gene_tracks_download <- downloadHandler(
        filename= function(){paste0("gene_tracks_", graphgene(), "_", Sys.Date(), ".png")},
        content=function(file){
            png(file, width=16, height=12, units="in", res=300)
            gene_tracks()
            dev.off()})
    
    output$table_download <- downloadHandler(filename = function() {paste0('Asthma&GC_data_summary_table_',graphgene(),"_", Sys.Date(), '.csv')},
                                             content = function(file) {write.csv(rbind(tableforgraph_Asthma(), tableforgraph_GC()), file, row.names=FALSE)})
    
    output$SNP_data_download <- downloadHandler(filename = function() {paste0('SNP_results_for_',graphgene(), Sys.Date(), '.csv')},
                                                content = function(file) {write.csv(snp_data(), file, row.names=FALSE)})
})
