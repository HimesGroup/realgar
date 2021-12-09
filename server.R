#detach("package:Gviz", unload=TRUE) # this is to keep RStudio happy - run if loading app more than once in same session - keep commented out otherwise
# if load Gviz 2x in same session (i.e. close & re-run app), get "object of type 'closure' is not subsettable" error
# should not be an issue when running app from the website
# cat(file=stderr(), as.character(Sys.time()),"packages start\n")
# use this type of command to easily see dataset loading time in RStudio  
# currently 3 seconds from "start package load" to "finish gene_locations load"
#rlang version 0.2.1

#detach("package:Gviz", unload=TRUE)
#library(Gviz, quietly = T)

# server
server <- shinyServer(function(input, output, session) {
  
  ## Loading proxy
  output$loadProxy <- renderUI({NULL})
  
  ## input values
  input_current <- eventReactive(input$go, {
    input$current
  }, ignoreNULL = FALSE
  )

  input_tissue <- eventReactive(input$go, {
    input$Tissue
  }, ignoreNULL = FALSE
  )   
  
  input_asthma <- eventReactive(input$go, {
    input$Asthma
  }, ignoreNULL = FALSE
  )
  
  input_treatment <- eventReactive(input$go, {
    input$Treatment
  }, ignoreNULL = FALSE
  ) 
  
  input_which_SNPs <- eventReactive(input$go, {
    input$which_SNPs
  }, ignoreNULL = FALSE
  ) 
  
  
  ## Get current gene or gene near selected SNP ID
   curr_gene <- reactive({
     current <- toString(input_current())
     if (gsub(" ", "", tolower(current), fixed=TRUE) %in% c(get_snp("snp"),get_snp("snp_gabriel"),get_snp("snp_fer"),
                                                            get_snp("snp_eve_all"),get_snp("snp_eve_ea"),get_snp("snp_eve_aa"),get_snp("snp_eve_la"),
                                                            get_snp("snp_TAGC_multi"), get_snp("snp_TAGC_euro"))){ #if SNP ID is entered, convert internally to nearest gene symbol
       rsid <- gsub(" ", "", tolower(current), fixed=TRUE)
       all_matches <- rbind(rbind(get_matches(rsid,"snp"),get_matches(rsid,"snp_gabriel"),get_matches(rsid,"snp_fer"),
                                  get_matches(rsid,"snp_eve_all"),get_matches(rsid,"snp_eve_ea"),get_matches(rsid,"snp_eve_aa"),get_matches(rsid,"snp_eve_la"),
                                  get_matches(rsid,"snp_TAGC_multi"), get_matches(rsid,"snp_TAGC_euro")))
       join_gene_snp(all_matches)
     } else {
       # if it is not in the list of snps, it is a gene id OR a snp that is not associated with asthma
       # in the latter case it will not show up in the list of genes & user gets an "enter valid gene/snp id" message
       gsub(" ", "", toupper(current), fixed = TRUE) #make uppercase, remove spaces
     }
   })

  
  GeneSymbol <- reactive({if (curr_gene() %in% gene_list) {TRUE} else {FALSE}})  #used later to generate error message when a wrong gene symbol is input
  
  ################################################
  ## reactive UI for EVE & TAGC p-value options ##
  ################################################
  
  # #if either EVE or TAGC SNPs selected, display the phrase "GWAS display options:"
  # output$GWAS_text <- reactive({if("snp_eve_subs" %in% input_which_SNPs() | "snp_TAGC_subs" %in% input_which_SNPs()){"GWAS display options:"} else {" "}})
  # 
  # #if EVE SNPs selected, display option to choose population
  # output$eve_options <- renderUI({if("snp_eve_subs" %in% input_which_SNPs()) {selectInput("which_eve_pvals", "Which EVE p-values to use?", 
  #                                                                                       list("All subjects"="meta_P", "African American"="meta_P_AA","European American"="meta_P_EA", "Latino"="meta_P_LAT"), 
  #                                                                                       selected="meta_P")} else {NULL}})
  # 
  # #if TAGC SNPs selected, display option to choose population
  # output$TAGC_options <- renderUI({if("snp_TAGC_subs" %in% input_which_SNPs()) {selectInput("which_TAGC_pvals", "Which TAGC p-values to use?", 
  #                                                                                         list("Multiancestry"="p_ran_multi", "European ancestry"="p_ran_euro"), 
  #                                                                                         selected="p_ran_multi")} else {NULL}})
  # #if UKBiobank SNPs selected, display phenotype of interest for association results
  # output$UKBB_options <- renderUI({if("snp_UKBB_subs" %in% input_which_SNPs()) {selectInput("which_UKBB_pheno", "Which UKBiobank phenotype to use?", 
  #                                                                                      list("Asthma"="Asthma", "COPD"="COPD","ACO"="ACO"), 
  #                                                                                      selected="pheno_asthma")} else {NULL}})
    
  #################################################################################
  ## "Select all" buttons for tissue, asthma, treatment and GWAS study selection ##
  #################################################################################
  
  #Tissue
  # observe({
  #   if(input$selectall_tissue == 0) return(NULL) # don't do anything if action button has been clicked 0 times
  #   else if (input$selectall_tissue%%2 == 0) { # %% means "modulus" - i.e. here you're testing if button has been clicked a multiple of 2 times
  #     #updateCheckboxGroupInput(session,"Tissue","Tissue:",choices=tissue_choices, selected = tissue_choices) 
  #     #updateActionButton(session, "selectall_tissue", label="Unselect all") # change action button label based on user input
  #   } else { # else is 1, 3, 5 etc.
  #     #updateCheckboxGroupInput(session,"Tissue","Tissue:",choices=tissue_choices)
  #     #updateActionButton(session, "selectall_tissue", label="Select all")
  #     updateSelectInput(session,"Tissue", label="Tissue:",choices=tissue_choices, selected = tissue_choices) # change mk
  #     updateActionButton(session, "selectall_tissue", label="Unselect all")
  #   }
  # })
  
  #Disease
  # observe({
  #   if(input$selectall_asthma == 0) return(NULL) 
  #   else if (input$selectall_asthma%%2 == 0) {
  #     #updateCheckboxGroupInput(session, "Asthma", label="Condition vs Healthy:", choices=asthma_choices)
  #     updateSelectInput(session, "Asthma", label="Condition vs Healthy:", choices=asthma_choices) # change mk
  #     updateActionButton(session, "selectall_asthma", label="Select all")
  #   }
  #   else {
  #     #updateCheckboxGroupInput(session,"Asthma",label="Condition vs Healthy:", choices=asthma_choices, selected = asthma_choices)
  #     updateSelectInput(session, "Asthma", label="Condition vs Healthy:", choices=asthma_choices, selected = asthma_choices) # change mk
  #     updateActionButton(session, "selectall_asthma", label="Unselect all")
  #   }})
  
  
  #Treatment
  # observe({
  #     if(input$selectall_treatment == 0) return(NULL) 
  #     else if (input$selectall_treatment%%2 == 0) {
  #       #updateCheckboxGroupInput(session, "Treatment", "Treatment:", choices = treatment_choices)
  #       #updateActionButton(session, "selectall_treatment", label="Select all")
  #       updateSelectInput(session, "Treatment", "Treatment:", choices = treatment_choices) # change mk
  #       updateActionButton(session, "selectall_treatment", label="Select all")
  #     }
  #     else {
  #       #updateCheckboxGroupInput(session, "Treatment", "Treatment:", choices = treatment_choices, selected = treatment_choices)
  #       #updateActionButton(session, "selectall_treatment", label="Unselect all")
  #       updateSelectInput(session, "Treatment", "Treatment:", choices = treatment_choices, selected = treatment_choices) # change mk
  #       updateActionButton(session, "selectall_treatment", label="Unselect all")
  #     }})
  
  #Smoking
  #observe({
  #  if(input$selectall_smoking == 0) return(NULL) 
  #  else if (input$selectall_smoking%%2 == 0) {
  #    updateCheckboxGroupInput(session, "Smoking", "Smoking:", choices = smoking_choices)
  #    updateActionButton(session, "selectall_smoking", label="Select all")
  #  }
  #  else {
  #    updateCheckboxGroupInput(session, "Smoking", "Smoking:", choices = smoking_choices, selected = smoking_choices)
  #    updateActionButton(session, "selectall_smoking", label="Unselect all")
  #  }})
  
  #GWAS
  # observe({
  #     if(input$selectall_GWAS == 0) return(NULL) 
  #     else if (input$selectall_GWAS%%2 == 0) {
  #         #updateCheckboxGroupInput(session,"which_SNPs","GWAS Results",choices=gwas_choices)
  #         #updateActionButton(session, "selectall_GWAS", label="Select all")
  #         updateSelectInput(session,"which_SNPs","GWAS Results",choices=gwas_choices) # change mk
  #         updateActionButton(session, "selectall_GWAS", label="Select all")
  #     }
  #     else {
  #         #updateCheckboxGroupInput(session,"which_SNPs","GWAS Results", choices=gwas_choices, selected=gwas_choices)
  #         #updateActionButton(session, "selectall_GWAS", label="Unselect all")
  #         updateSelectInput(session,"which_SNPs","GWAS Results", choices=gwas_choices, selected=gwas_choices) # change mk
  #         updateActionButton(session, "selectall_GWAS", label="Unselect all")
  #     }})
  
    
  #######################
  ## GEO studies table ##
  #######################
  #select GEO studies matching desired conditions;
  #Jessica's initial app had an "and" condition here; Maya changed it to "or"
  # Mengyuan changed it to: if only select tissue or asthma/treatment, will use all the available studies; else use the intersection
  UserDataset_Info <- reactive({
    Dataset_Info_Tissue <- subset(Dataset_Info, Dataset_Info$Tissue %in% input_tissue())
    
    if(is.null(input_asthma())){Dataset_Info_A <- subset(Dataset_Info, Dataset_Info$Asthma %in% input_asthma())} # change_mk
    else {Dataset_Info_A <- subset(Dataset_Info, Dataset_Info$Asthma %in% input_asthma())}

    if(is.null(input_treatment())){Dataset_Info_B <- subset(Dataset_Info, Dataset_Info$Asthma %in% input_treatment())}
    else {Dataset_Info_B <- subset(Dataset_Info, Dataset_Info$Asthma %in% input_treatment())}
    
    Dataset_Info_Asthma <- rbind(Dataset_Info_A,Dataset_Info_B)
    
    if ((nrow(Dataset_Info_Tissue)==0)|(nrow(Dataset_Info_Asthma)==0)) {Dataset_Info1 <- rbind(Dataset_Info_Tissue,Dataset_Info_Asthma)}
    else {Dataset_Info1 <- subset(Dataset_Info_Tissue,Dataset_Info_Tissue$Unique_ID%in%Dataset_Info_Asthma$Unique_ID)}
    
    #BA_PDE
    if(length(setdiff(c("BA","PDE"),input_treatment()))==0){Dataset_Info1 <- rbind(Dataset_Info1,BA_PDE_Info)}
    
    #Return
    Dataset_Info1
    
  }) %>% debounce(1000)
  
  # gene expression (GEO) studies table
  #add links for GEO_ID and PMID
  GEO_data <- reactive({
      validate(need(nrow(UserDataset_Info()) != 0, "No gene expression datasets selected")) #Generate a error message when no data is loaded.
      
      UserDataset_Info() %>%
          dplyr::mutate(GEO_ID_link = ifelse(grepl("SRP", GEO_ID), #GEO link is conditional on whether GEO_ID is an "SRP" or "GSE"
                                             paste0("https://www.ncbi.nlm.nih.gov/sra/?term=", GEO_ID), 
                                             paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GEO_ID)),
                        PMID_link = paste0("https://www.ncbi.nlm.nih.gov/pubmed/?term=", PMID),
                        QC_link = ifelse(grepl("SRP", GEO_ID), #QC link is conditional on whether GEO_ID is an "SRP" or "GSE"
                               # paste0("http://public.himeslab.org/realgar_qc/",GEO_ID,"_QC_RnaSeqReport.html"),
                               # paste0("http://public.himeslab.org/realgar_qc/",GEO_ID,"_QC_report.html")))})
                               paste0(GEO_ID,"_QC_RnaSeqReport.html"), 
                               paste0(GEO_ID,"_QC_report.html")))})

  
  GEO_links <- reactive({
      GEO_Dataset <- paste0("<a href='",  GEO_data()$GEO_ID_link, "' target='_blank'>",GEO_data()$GEO_ID,"</a>")
      GEO_PMID <- paste0("<a href='",  GEO_data()$PMID_link, "' target='_blank'>",GEO_data()$PMID,"</a>")
      GEO_Report <- paste0("<a href='",  GEO_data()$QC_link, "' target='_blank'>",GEO_data()$Report,"</a>")
      
      df <- data.frame(GEO_Dataset, GEO_PMID, GEO_Report, GEO_data()$Description)
      colnames(df) <- c("Dataset", "PMID", "Report","Description")
      df
  })
  
  # gwas studies table

  GWAS_data <- reactive({
    sub_convert = c("snp_subs", "snp_gabriel_subs", "snp_fer_sub",
                    "snp_eve_subs", "snp_eve_subs", "snp_eve_subs", "snp_eve_subs",
                    "snp_TAGC_subs", "snp_TAGC_subs")
    names(sub_convert) = c(
    "snp_subs", "snp_gabriel_subs", "snp_fer_sub",
    "snp_eve_all_subs" , "snp_eve_ea_subs", "snp_eve_aa_subs", "snp_eve_la_subs",
    "snp_TAGC_multi_subs", "snp_TAGC_euro_subs")
    sub_selected <- unname(sapply(input_which_SNPs(), function(x){sub_convert[x]}))
    df <- GWAS_Dataset_Info[which(GWAS_Dataset_Info$Tissue %in% sub_selected),c("GEO_ID", "PMID","Description")]
    validate(need(nrow(df) != 0, "No GWAS datasets selected")) #Generate a error message when no data is loaded.
    colnames(df) <- c("Dataset", "Link","Description") # I put the link for the study into the PMID column of the spreadsheet for convenience - change later?
    df
  })
  
  GWAS_links <- reactive ({
    GWAS_Dataset <- paste0("<a href='",  GWAS_data()$Link, "' target='_blank'>",GWAS_data()$Dataset,"</a>")
    df <- data.frame(GWAS_Dataset, GWAS_data()$Description)
    colnames(df) <- c("Dataset", "Description")
    df
  })

  
  #table output for "Datasets Loaded" tab
  output$GEO_table <- DT::renderDataTable(GEO_links(),  
                                          class = 'cell-border stripe', 
                                          rownames = FALSE, 
                                          options = list(paging = FALSE, searching = FALSE),
                                          escape = FALSE)
  
  output$GWAS_table <- DT::renderDataTable(GWAS_links(),
                                           class = 'cell-border stripe',
                                           rownames = FALSE,
                                           options = list(paging = FALSE, searching = FALSE), 
                                           escape = FALSE)
  
  #########################################
  ## Select GEO data for plots and table ##
  #########################################
  
  #select and modify data used for plots and accompanying table
  output.tableforplot <- reactive({
    
    validate(need(nrow(UserDataset_Info()) != 0, "Please choose at least one dataset.")) #Generate a error message when no data is loaded.
    validate(need(curr_gene() != "", "Please enter a gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
    
    gene <- paste0('"',curr_gene(),'"')
    query <- paste0("SELECT * FROM REALGAR WHERE Gene = ",gene)
    res <- dbSendQuery(ngs_db, query)
    out.table <- dbFetch(res)
    dbClearResult(res)
    
    # Generate error message if the gene symbol is not right.
    validate(need(GeneSymbol() != FALSE, "Please enter a valid gene symbol or SNP ID."))
    out.table <- out.table %>%
      dplyr::filter(Unique_ID %in% UserDataset_Info()$Unique_ID) %>%
      dplyr::rename(adj.P.Val = adjPVal,P.Value = PValue) %>%
      arrange(desc(Fold_Change),desc(Upper_bound_CI))
  })
  
  
  ###################################
  ## Data table for meta-analysis  ##
  ###################################
  
  # asthma
  data_Asthma <- reactive({ 
  output.tableforplot_asthma <- output.tableforplot() %>% dplyr::filter(App == "asthma")
  })
  
  # GC
  data_GC <- reactive({ 
  output.tableforplot_GC <- output.tableforplot() %>% dplyr::filter(App == "GC")
  })
  
  # Smoking
  data_cig <- reactive({ 
    output.tableforplot_cig <- output.tableforplot() %>% dplyr::filter(App == "Smoking"|App == "Pollutant")
  })

  
  ###################################
  ##        Combined p-values      ##
  ###################################
  
  # asthma
  asthma_pcomb <- reactive({
    dat <- data_Asthma()
    if (length(dat$adj.P.Val)>1) {
      #write.table(dat,paste0("Asthma_",curr_gene(),".txt"),col.names=T,row.names=F,sep="\t",quote=F)
      asthma_rankprod_pcomb <- rankprod_stat(dat)
      asthma_sumlog_pcomb <- sumlog_stat(dat)
      pcomb_text=paste0("P-value-based integration = ", asthma_sumlog_pcomb, "; Rank-based integration = ", asthma_rankprod_pcomb)
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
      GC_rankprod_pcomb <- rankprod_stat(dat)
      GC_sumlog_pcomb <- sumlog_stat(dat)
      pcomb_text=paste0("P-value-based integration = ", GC_sumlog_pcomb, "; Rank-based integration = ", GC_rankprod_pcomb)
    }
    else {pcomb_text=""}
    pcomb_text
  })
  
  output$GC_pcomb_text <- renderText({GC_pcomb()})
  
  # Smoking
  cig_pcomb <- reactive({
    dat <- data_cig()
    if (length(dat$adj.P.Val)>1) {
      #write.table(dat,paste0("GC_",curr_gene(),".txt"),col.names=T,row.names=F,sep="\t",quote=F)
      cig_rankprod_pcomb <- rankprod_stat(dat)
      cig_sumlog_pcomb <- sumlog_stat(dat)
      pcomb_text=paste0("P-value-based integration = ", cig_sumlog_pcomb, "; Rank-based integration = ", cig_rankprod_pcomb)
    }
    else {pcomb_text=""}
    pcomb_text
  })
  
  output$cig_pcomb_text <- renderText({cig_pcomb()})
  
  # pollutant
  
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
  
  meta_cig <- reactive({
    dat <- data_cig()
    if (length(dat$adj.P.Val)>1) {
      res <- meta_stat(dat)
    }
    else {res <- list(meta_pval=NA,meta_fc=NA,meta_lower=NA,meta_upper=NA)}
    res
  })  
  
  
  ###################################
  ## Data table accompanying plots ##
  ###################################
  
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
  
  # asthma
  data2_Asthma <- reactive({ # dataset without meta-analysis results
    data_Asthma()%>%
      dplyr::select(Unique_ID, adj.P.Val, P.Value,Fold_Change, neglogofP, Lower_bound_CI, Upper_bound_CI) %>%
      dplyr::arrange(desc(Fold_Change),desc(Upper_bound_CI)) %>% # sort by first effect size (fold change) and then by upper CI in a descending order
      dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),
                    adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), 
                    P.Value =format(P.Value, scientific=TRUE, digits=3), 
                    Lower_bound_CI = round(Lower_bound_CI, digits = 2), 
                    Upper_bound_CI = round(Upper_bound_CI, digits = 2), 
                    Comparison = "Asthma vs. non-asthma")%>%
      dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Fold Change`=Fold_Change)})
  
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
          dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),
                        adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), 
                        P.Value =format(P.Value, scientific=TRUE, digits=3), 
                        Lower_bound_CI = round(Lower_bound_CI, digits = 2), 
                        Upper_bound_CI = round(Upper_bound_CI, digits = 2), 
                        Comparison = "Treatment vs. at baseline")%>%
          dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Fold Change`=Fold_Change)})
  
  # obtain intermediate GC data for table output and forest plots
  data3_GC <- reactive({
    interdata_func(data2_GC(),meta_GC())
  }) 
  
  # output table for GC
  tableforgraph_GC <- reactive({
      tableforgraph_func(data3_GC())
  })
  
  
  # Smoking
  data2_cig <- reactive({
    data_cig()%>%
      dplyr::select(Unique_ID, adj.P.Val, P.Value,Fold_Change, neglogofP, Lower_bound_CI, Upper_bound_CI) %>%
      dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),
                    adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), 
                    P.Value =format(P.Value, scientific=TRUE, digits=3), 
                    Lower_bound_CI = round(Lower_bound_CI, digits = 2), 
                    Upper_bound_CI = round(Upper_bound_CI, digits = 2), 
                    Comparison = "Pollutant vs. at baseline")%>%
      dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Fold Change`=Fold_Change)})
  
  # obtain intermediate smoking data for table output and forest plots
  data3_cig <- reactive({
    interdata_func(data2_cig(),meta_cig())
  }) 
  
  # output table for Smoking
  tableforgraph_cig <- reactive({
    tableforgraph_func(data3_cig())
  })
  
  
  #combine asthma & GC & Smoking into one
  output$tableforgraph <- DT::renderDataTable(rbind(rbind(tableforgraph_Asthma(), tableforgraph_GC()),tableforgraph_cig()),
                                              class = 'cell-border stripe', 
                                              rownames = FALSE, 
                                              options = list(paging = FALSE, searching = FALSE),
                                              width = "100%")
  
  ################################################
  ## reactive UI for forestplot download options ##
  ################################################
  
  #if EVE SNPs selected, display option to choose population
  #output$asthma_fp_download <- renderUI({if(nrow(data3_Asthma())!= 0) {downloadButton(outputId="asthma_fc_download",label="Download asthma forest plot")} else {NULL}})
  
  #if TAGC SNPs selected, display option to choose population
  #output$GC_fp <- reactive({if(nrow(data3_GC())!= 0) {" "} else {NULL}})
  
  
  #################
  ## Forestplots ##
  #################
  
  # asthma forestplot
  
  getHeightAsthma <- reactive({
    getHeight_func(data3_Asthma())
  })
  
  forestplot_asthma <- reactive({
    forestplot_func(data3_Asthma(),"Asthma Transcriptomic Results for ",curr_gene())
  }) 
  
  # GC forestplot
  
  getHeightGC <- reactive({
    getHeight_func(data3_GC())
  })
  
  forestplot_GC <- reactive({
    forestplot_func(data3_GC(),"Exposure Transcriptomic Results for ",curr_gene())
  }) 
  
  # Smoking forestplot
  
  getHeightcig <- reactive({
    getHeight_func(data3_cig())
  })
  
  forestplot_cig <- reactive({
    forestplot_func(data3_cig(),"Smoking Transcriptomic Results for ",curr_gene())
  }) 
  
  
  output$forestplot_asthma = renderPlot({forestplot_asthma()}, height=getHeightAsthma)
  
  output$forestplot_GC = renderPlot({forestplot_GC()}, height=getHeightGC)
  
  output$forestplot_cig = renderPlot({forestplot_cig()}, height=getHeightcig)
  
  ###############################
  ## Gene, SNP and TFBS tracks ##
  ###############################
 
  observe({
    updateSelectizeInput(
      session, "pval_thr", server = TRUE,
      choices = pval_for_select,
      selected = pval_select,
      
      options = list(render = I(
        '{ option: function(item, escape) {
            return "<div>" + item.html + "</div>";
            }
         }'
      ))
    )
  })
   
  #horizontal color scale for gene tracks
  output$color_scale3 <- renderImage({ 
    return(list(
      src = "realgar_data/color_scale_horizontal.png",
      height=109*1.05,
      #width=1015*1.05,
      width=818*1.05,
      filetype = "image/png",
      alt = "color_scale"))}, deleteFile = FALSE)
  
  #filter data for selected gene
  # gene_subs <- reactive({
  #   gene_subs_temp <- get_query_db("gene_locations",curr_gene())
  #   gene_subs_temp <- gene_subs_temp[!(duplicated(gene_subs_temp$exon)),]
  #   })
  
  # tfbs_subs <- reactive({
  #   tfbs_db <- get_query_db("tfbs",curr_gene())
  #   unique(tfbs_db)
  #   })
  
  snp_subs <- reactive({
    print(input$pval_thr)
    if(("snp_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_db <- get_query_db("snp",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_db <- get_query_nominal_db("snp",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_db <- get_query_genomewide_db("snp",curr_gene())
      } 
      if (nrow(snp_db) > 0){
        snp_db %>% dplyr::mutate(start=start+1) %>% unique()
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    }
    else {data.frame(matrix(nrow = 0, ncol = 0))}}) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  
  snp_gabriel_subs <- reactive({
    print("select gabriel")
    if(("snp_gabriel_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_gabriel_db <- get_query_db("snp_gabriel",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_gabriel_db <- get_query_nominal_db("snp_gabriel",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_gabriel_db <- get_query_genomewide_db("snp_gabriel",curr_gene())
      }
      if (nrow(snp_gabriel_db) > 0){
        snp_gabriel_db %>% dplyr::mutate(start=start+1, pmid=20860503) %>% 
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle) %>% unique()
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
       }
    else {data.frame(matrix(nrow = 0, ncol = 0))}}) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  
  snp_fer_subs <- reactive({
    print("select ferreira")
    if(("snp_fer_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_fer_db <- get_query_db("snp_fer",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_fer_db <- get_query_nominal_db("snp_fer",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_fer_db <- get_query_genomewide_db("snp_fer",curr_gene())
      } 
      if (nrow(snp_fer_db) > 0){
        snp_fer_db %>% dplyr::mutate(start=start+1, pmid=29083406) %>% 
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle)
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    }
    else {data.frame(matrix(nrow = 0, ncol = 0))}}) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  
  snp_eve_all_subs <- reactive({
    print("select eve all")
    if(("snp_eve_all_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_eve_all_db <- get_query_db("snp_eve_all",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_eve_all_db <- get_query_nominal_db("snp_eve_all",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_eve_all_db <- get_query_genomewide_db("snp_eve_all",curr_gene())
      } 
      if (nrow(snp_eve_all_db)>0){
        snp_eve_all_db %>%
            dplyr::mutate(start=start+1, pmid=21804549)%>%
            dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle)
      }
      else {data.frame(matrix(nrow = 0, ncol = 0))} #since pval_selector might remove all rows
    }else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_eve_pvals is null
  }) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  

  snp_eve_ea_subs <- reactive({
    print("select eve ea")
    if(("snp_eve_ea_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_eve_ea_db <- get_query_db("snp_eve_ea",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_eve_ea_db <- get_query_nominal_db("snp_eve_ea",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_eve_ea_db <- get_query_genomewide_db("snp_eve_ea",curr_gene())
      } 
      if (nrow(snp_eve_ea_db)>0){
        snp_eve_ea_db %>%
          dplyr::mutate(start=start+1, pmid=21804549)%>%
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle)
      }
      else {data.frame(matrix(nrow = 0, ncol = 0))} #since pval_selector might remove all rows
    }else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_eve_pvals is null
  }) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  
  snp_eve_aa_subs <- reactive({
    print("select eve aa")
    if(("snp_eve_aa_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_eve_aa_db <- get_query_db("snp_eve_aa",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_eve_aa_db <- get_query_nominal_db("snp_eve_aa",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_eve_aa_db <- get_query_genomewide_db("snp_eve_aa",curr_gene())
      } 
      if (nrow(snp_eve_aa_db)>0){
        snp_eve_aa_db %>%
          dplyr::mutate(start=start+1, pmid=21804549)%>%
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle)
      }
      else {data.frame(matrix(nrow = 0, ncol = 0))} #since pval_selector might remove all rows
    }else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_eve_pvals is null
  }) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  
  
  snp_eve_la_subs <- reactive({
    print("select eve la")
    if(("snp_eve_la_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_eve_la_db <- get_query_db("snp_eve_la",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_eve_la_db <- get_query_nominal_db("snp_eve_la",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_eve_la_db <- get_query_genomewide_db("snp_eve_la",curr_gene())
      } 
      if (nrow(snp_eve_la_db)>0){
        snp_eve_la_db %>%
          dplyr::mutate(start=start+1, pmid=21804549)%>%
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle)
      }
      else {data.frame(matrix(nrow = 0, ncol = 0))} #since pval_selector might remove all rows
    }else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_eve_pvals is null
  }) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  

  snp_TAGC_multi_subs <- reactive({
    print("select tagc multi")
    if(("snp_TAGC_multi_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_TAGC_multi_db <- get_query_db("snp_TAGC_multi",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_TAGC_multi_db <- get_query_nominal_db("snp_TAGC_multi",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_TAGC_multi_db <- get_query_genomewide_db("snp_TAGC_multi",curr_gene())
      } 
      
      if (nrow(snp_TAGC_multi_db)>0){
        snp_TAGC_multi_db %>% 
            dplyr::mutate(start=start+1, pmid=29273806) %>% 
            dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle)
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    } else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_TAGC_pvals in NULL
  }) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error

  snp_TAGC_euro_subs <- reactive({
    print("select tagc euro")
    if(("snp_TAGC_euro_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_TAGC_euro_db <- get_query_db("snp_TAGC_euro",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_TAGC_euro_db <- get_query_nominal_db("snp_TAGC_euro",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_TAGC_euro_db <- get_query_genomewide_db("snp_TAGC_euro",curr_gene())
      } 
      
      if (nrow(snp_TAGC_euro_db)>0){
        snp_TAGC_euro_db %>% 
          dplyr::mutate(start=start+1, pmid=29273806) %>% 
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle)
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    } else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_TAGC_pvals in NULL
  }) #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error  

  snp_UKBB_asthma_subs <- reactive({
    print("select ukbb asthma")
    if(("snp_UKBB_asthma_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_UKBB_asthma_db <- get_query_db("snp_UKBB_asthma",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_UKBB_asthma_db <- get_query_nominal_db("snp_UKBB_asthma",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_UKBB_asthma_db <- get_query_genomewide_db("snp_UKBB_asthma",curr_gene())
      }

      if (nrow(snp_UKBB_asthma_db)>0){
        snp_UKBB_asthma_db %>% 
            dplyr::mutate(start=start+1, pmid=NA) %>% 
            dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle) 
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    } else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_UKBB_pheno in NULL
  }) 
  
  snp_UKBB_copd_subs <- reactive({
    if(("snp_UKBB_copd_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_UKBB_copd_db <- get_query_db("snp_UKBB_copd",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_UKBB_copd_db <- get_query_nominal_db("snp_UKBB_copd",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_UKBB_copd_db <- get_query_genomewide_db("snp_UKBB_copd",curr_gene())
      }
      
      if (nrow(snp_UKBB_copd_db)>0){
        print("select ukbb copd")
        snp_UKBB_copd_db %>% 
          dplyr::mutate(start=start+1, pmid=NA) %>% 
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle) 
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    } else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_UKBB_pheno in NULL
  }) 
  
  
  snp_UKBB_aco_subs <- reactive({
    if(("snp_UKBB_aco_subs" %in% input_which_SNPs())) {
      if (input$pval_thr=="normal") {
        snp_UKBB_aco_db <- get_query_db("snp_UKBB_aco",curr_gene())
      } else if (input$pval_thr=="norminal") {
        snp_UKBB_aco_db <- get_query_nominal_db("snp_UKBB_aco",curr_gene())
      } else if (input$pval_thr=="genomewide") {
        snp_UKBB_aco_db <- get_query_genomewide_db("snp_UKBB_aco",curr_gene())
      }
      
      if (nrow(snp_UKBB_aco_db)>0){
        print("select ukbb aco")
        snp_UKBB_aco_db %>% 
          dplyr::mutate(start=start+1, pmid=NA) %>% 
          dplyr::select(chromosome, start, end, snp, symbol, pmid, meta_P, neg_log_p, color, eqtl_lung, eqtl_blood, eqtl_muscle) 
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    } else {data.frame(matrix(nrow = 0, ncol = 0))} #since which_UKBB_pheno in NULL
  }) 
  
    
  #plot height increases if more tracks are displayed
  # observe({output$gene_tracks_outp2 <- renderPlot({gene_tracks()}, width=1055,
  #                                                 height=400 + 15*length(unique(gene_subs()$transcript)) + 30*nrow(snp_eve_subs()) + 10*(nrow(snp_subs())+nrow(snp_gabriel_subs())+nrow(snp_fer_subs())+nrow(snp_TAGC_subs())))})
  
  #################################
  ## SNP data table for download ##
  #################################
  snp_subs_temp <- reactive({
    if (nrow(snp_subs()) > 0) {
      snp_subs()%>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "GRASP") %>%
        dplyr::select(chromosome, snp, symbol, position, pmid, source, meta_P, eqtl_lung, eqtl_blood, eqtl_muscle)
    }
  })
  
  
  snp_gabriel_subs_temp <- reactive({
    if (nrow(snp_gabriel_subs()) > 0) {
      snp_gabriel_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "GABRIEL") %>%
        dplyr::select(-c(end, color, neg_log_p))
    }
  })
  
  snp_fer_subs_temp <- reactive({
    if (nrow(snp_fer_subs()) > 0) {
      snp_fer_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "Ferreira") %>%
        dplyr::select(-c(end, color, neg_log_p))
    }
  })
  
  snp_eve_all_subs_temp <- reactive({
    if (nrow(snp_eve_all_subs()) > 0) {
      snp_eve_all_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "EVE_ALL") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  })

  
  snp_eve_ea_subs_temp <- reactive({
    if (nrow(snp_eve_ea_subs()) > 0) {
      snp_eve_ea_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "EVE_EA") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  })  
 
  
  snp_eve_aa_subs_temp <- reactive({
    if (nrow(snp_eve_aa_subs()) > 0) {
      snp_eve_aa_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "EVE_AA") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  })
  
  snp_eve_la_subs_temp <- reactive({
    if (nrow(snp_eve_la_subs()) > 0) {
      snp_eve_la_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "EVE_LA") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  })
  
  snp_TAGC_multi_subs_temp <- reactive({
    if (nrow(snp_TAGC_multi_subs()) > 0) {
      snp_TAGC_multi_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "TAGC_Multi") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  }) 
  
  snp_TAGC_euro_subs_temp <- reactive({
    if (nrow(snp_TAGC_euro_subs()) > 0) {
      snp_TAGC_euro_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "TAGC_EUR") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  }) 
  
  snp_UKBB_asthma_subs_temp <- reactive({
    if (nrow(snp_UKBB_asthma_subs()) > 0) {
      snp_UKBB_asthma_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "UKBB_asthma") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  }) 
  
  snp_UKBB_copd_subs_temp <- reactive({
    if (nrow(snp_UKBB_copd_subs()) > 0) {
      snp_UKBB_copd_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "UKBB_COPD") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  })
  
  snp_UKBB_aco_subs_temp <- reactive({
    if (nrow(snp_UKBB_aco_subs()) > 0) {
      snp_UKBB_aco_subs() %>%
        dplyr::rename(position=start) %>%
        dplyr::mutate(source = "UKBB_ACO") %>%
        dplyr::select(-c(end, neg_log_p, color))
    }
  })
  
  #snp_data <- reactive({dplyr::bind_rows(snp_subs_temp(), snp_eve_subs_temp(), snp_gabriel_subs_temp(), snp_fer_subs_temp(), snp_TAGC_subs_temp())})
  snp_data <- reactive({
    print("combine gwas datasets")
    dplyr::bind_rows(snp_subs_temp(), snp_gabriel_subs_temp(), snp_fer_subs_temp(),
                                         snp_eve_all_subs_temp(), snp_eve_ea_subs_temp(), snp_eve_aa_subs_temp(), snp_eve_la_subs_temp(), 
                                         snp_TAGC_multi_subs_temp(), snp_TAGC_euro_subs_temp(), 
                                         snp_UKBB_asthma_subs_temp(), snp_UKBB_copd_subs_temp(), snp_UKBB_aco_subs_temp())
    })
  
  #########################
  ## KARYOPLOTR TRACKS ##
  ########################
  
  gene.region <- reactive({
    validate(need(curr_gene() != "", "Please enter a gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
    region <- all_genes %>% dplyr::filter(symbol == curr_gene()) 
    #region_split <- strsplit(x=gsub("\\s","",region$Pos), split="-|:| ")[[1]]
    if (nrow(region) > 0){
      unique(region)
    } else {NULL}
  })
  
  ## Make Plot
  # output$karyoPlot <- renderPlot({
  #   validate(need(!is.null(gene.region()), "Please enter a valid gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
  #   make_karyoplot(gene.region(),snp_subs(), snp_eve_subs(), snp_gabriel_subs(), snp_fer_subs(), snp_TAGC_subs())
  # })
  
  ## Make Plot
  karyoplot <- reactive({
    validate(need(!is.null(gene.region()), "Please enter a valid gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
    print("generate karyoplot")
    make_karyoplot(gene.region(),snp_subs(), snp_gabriel_subs(), snp_fer_subs(),
                   snp_eve_all_subs(), snp_eve_ea_subs(), snp_eve_aa_subs(), snp_eve_la_subs(), 
                   snp_TAGC_multi_subs(), snp_TAGC_euro_subs(),
                   snp_UKBB_asthma_subs(), snp_UKBB_copd_subs(), snp_UKBB_aco_subs())
  })
  
  observe({output$karyoPlot <- renderPlot({karyoplot()})})
  
  # GR-binding sites within the selected region
  GRbinding_table <- reactive({
    if (!is.null(gene.region())) {
      gene.sel <- gene.region()
      # generate GRanges object for gene regions
      if (gene.sel$strand == -1){
        gene.sel.gr <- toGRanges(data.frame(chr=gene.sel$Chromosome, start=gene.sel$Start, end=gene.sel$End+20000,genome="hg38"))
      } else if (gene.sel$strand == 1){
        gene.sel.gr <- toGRanges(data.frame(chr=gene.sel$Chromosome, start=gene.sel$Start-20000, end=gene.sel$End,genome="hg38"))
      }
    # generate GRanges object for GR-binding sites
    GRbinding_gr <- makeGRangesFromDataFrame(GRbinding[,1:3])
    overlap_gr <- findOverlaps(GRbinding_gr, gene.sel.gr)
    row_numbers = queryHits(overlap_gr)
    GRbinding[row_numbers, ] %>% data.frame()
    }
  })
  
  ###############################
  ##   Gene Expression Levels  ##
  ###############################
  
  #Gene input
  output$gene_name_out <- renderPrint(curr_gene())
  curr_gene_te <- reactive(curr_gene())
  
  #used to make situation-specific error messages
  in_all <- reactive({if (!(curr_gene_te() %in% all_genes_te$name)) {FALSE} else {TRUE}})  # wrong symbol was input
  in_unfiltered <- reactive({if ((curr_gene() %in% all_genes_te$name) & !(curr_gene_te() %in% unfiltered_genes$x)) {FALSE} else {TRUE}}) # gene not in database
  in_deseq2_filtered <- reactive({if ((curr_gene_te() %in% unfiltered_genes$x) & !(curr_gene_te() %in% deseq2_filtered_genes)) {FALSE} else {TRUE}}) # didn't pass sleuth filter
  
  output$gene_te <- renderPrint({
    curr_gene_te()
  })
  
  # more information about the dataset selected
  output$studyText <- renderUI({
    # if (!is.null(input$debugcode) && (input$debugcode == "studyText")) {
    #   browser()
    # }
    sras %>% filter(SRA_ID == input$tissue_te) %$% 
      p("Data used is available in the SRA under accession ",
        a(paste0(SRA_ID, ","), href=paste0("https://www.ncbi.nlm.nih.gov/sra/?term=", SRA_ID), target="_blank"),
        "and corresponds to ",
        Description,
        "More details were published ",
        a("here.", href=paste0("https://www.ncbi.nlm.nih.gov/pubmed/?term=", PMID), target="_blank"))
  })
  
  #generate faceted boxplot for the gene selected, using kallisto TPMs
  curr_data_te_tb <- reactive({
    x <- sras %>% 
      filter(SRA_ID == input$tissue_te) %$% 
      SRA_ID
    # obtain count data
    count_df <- get_query_counts(paste0(x, "_count"), curr_gene_te()) %>%
      tidyr::gather(Sample, value, -Gene, -gene_symbol) %>%
      data.frame()
    # obtain phenotype data
    pheno_df <- get_query_pheno(paste0(x, "_pheno"))
    pheno_df$Status <- gsub("-", "_", pheno_df$Status)
    pheno_df$Status <- gsub("1R6Fcig", "cig1R6F", pheno_df$Status) # avoid having numbers at the beginning of the characters
    
    merge(count_df, pheno_df, by="Sample")
  })
  
  getGeneBoxPlot <-reactive({
    validate(need(curr_gene_te() != "", "Please enter a gene id")) # no gene symbol was input
    validate(need(in_all() != FALSE, "Please enter a valid gene id.")) # invalid gene symbol was input
    validate(need(in_unfiltered() != FALSE, "This gene is not in the reference database.")) # gene not in database
    validate(need(in_deseq2_filtered() != FALSE, "Gene did not pass our DESeq2 filter (total counts should be greater than 10).")) # gene did not pass sleuth filter
    
    # This may be because less than , or all of the gene .
    curr_data_te <- curr_data_te_tb()
    #curr_data_te <- tpms[[x]] %>% filter(gene_symbol == curr_gene_te()) 
    in_filtered <- reactive({if (nrow(curr_data_te) == 0) {FALSE} else {TRUE}}) # didn't pass sleuth filter
    
    #gene plot initialize
    if (nrow(curr_data_te) > 0 & length(unique(curr_data_te$Donor))<=10) { # this iteration of curr_data has already been filtered to only have average_tpm > 1, with #donor<=10
      #gene_plot <- ggplot(curr_data_te, aes(x = Status, y = value, fill=Status)) + 
      gene_plot <- ggplot(curr_data_te, aes(x = Status, y = value)) + 
        geom_boxplot(outlier.colour=NA, lwd=0.2, color="grey18",fill="#1B9E77") + 
        stat_boxplot(geom ='errorbar', color="grey18") + 
        geom_jitter(aes(shape=Donor),size=2, width=0.2) +
        scale_shape_manual(values=seq(0,length(curr_data_te$Donor))) + 
        facet_wrap(~Gene) + 
        guides(fill=FALSE) + 
        theme_bw() +  
        labs(title=curr_gene_te()) + 
        labs(x="condition") + labs(y="Normalized Read Count") + 
        theme(#text = element_text(size=14), 
              strip.background = element_rect(colour="black",fill="#fbf7f5"),
              strip.text.x = element_text(size = 16), 
              #axis.text.x = element_text(angle = 90, hjust = 1, size=12),
              axis.text.x = element_text(angle=45, hjust=1, size=16),
              axis.text.y = element_text(size=14),
              title = element_text(size=16, face="bold.italic"),
              axis.title.x = element_text(size=16,face="bold"),
              legend.text=element_text(size=12),
              axis.title.y = element_text(size=16,face="bold"))
    }
    else if (nrow(curr_data_te) >0 & length(curr_data_te$Donor)>10) { # this iteration of curr_data has already been filtered to only have average_tpm > 1
      #gene_plot <- ggplot(curr_data_te, aes(x = Status, y = value, fill=Status)) + 
      gene_plot <- ggplot(curr_data_te, aes(x = Status, y = value)) + 
        geom_boxplot(outlier.colour=NA, lwd=0.2, color="grey18",fill="#1B9E77") + 
        stat_boxplot(geom ='errorbar', color="grey18") + 
        geom_jitter(size=2, width=0.2) +
        facet_wrap(~Gene) + 
        guides(fill=FALSE) + 
        theme_bw() +  
        labs(title=curr_gene_te()) + 
        labs(x="condition") + labs(y="Normalized Read Count") + 
        theme(#text = element_text(size=14), 
          strip.background = element_rect(colour="black",fill="#fbf7f5"),
          strip.text.x = element_text(size = 16), 
          #axis.text.x = element_text(angle = 90, hjust = 1, size=12),
          axis.text.x = element_text(angle=45, hjust=1, size=16),
          axis.text.y = element_text(size=14),
          title = element_text(size=16, face="bold.italic"),
          axis.title.x = element_text(size=16,face="bold"),
          legend.text=element_text(size=12),
          axis.title.y = element_text(size=16,face="bold"))
      }
      if (nrow(curr_data_te) > 0) {gene_plot}
  })
  
  #output boxplot
  output$GeneBoxPlot <- renderPlot(width = 650, height = 500, {
    # if (!is.null(input$debugcode) && (input$debugcode == "geneBoxPlot")) {
    #   browser()
    # }
    gbp <- getGeneBoxPlot()
    gbp
  })
  
  #output accompanying sleuth results
  output$table_title <- renderUI({
    title_string <- paste0("Differential Expression Results for ", curr_gene_te())
    HTML(title_string)
  })
  
  output$diffResults <- renderDataTable({
    # if (!is.null(input$debugcode) && (input$debugcode == "diffResults")) {
    #   browser()
    # }
  
    # count and phenotype results
    curr_data_te <- curr_data_te_tb() # count and phenotype table
    # summarize mean counts
    mean_count_df <- curr_data_te %>%
      dplyr::group_by(Gene, Status) %>%
      dplyr::summarise(mean=mean(value)) %>%
      data.frame() %>%
      tidyr::unite(name, c("Gene", "Status"), sep="_")
    count_mean <- mean_count_df$mean
    names(count_mean) <- mean_count_df$name

    # DE results
    x <- sras %>% 
      filter(SRA_ID == input$tissue_te) %$% 
      SRA_ID   
    de_df <- data.frame()
    gene_ids <- unique(curr_data_te$Gene)
    for (gene_id in gene_ids) {
      de_df_temp <- get_query_DE(paste0(x, "_de"), gene_id)
      de_df <- rbind(de_df, de_df_temp)
    }
    
    # modify DE results from wide to long
    de_df <- de_df %>%
      tidyr::gather(stat, value, -Gene) %>%
      dplyr::mutate(stat=stri_replace_last_fixed(stat, '.', 'statstartfromhere')) %>% # replace the last "." with 'statstartfromhere'
      tidyr::separate(stat, c("Comparison", "name"), sep="statstartfromhere") %>%
      tidyr::spread(name, value) %>%
      dplyr::mutate(cond1=gsub("(.*)_vs_.*", "\\1", Comparison), cond2=gsub(".*_vs_(.*)", "\\1", Comparison)) %>%
      dplyr::mutate(cond1_temp=paste0(Gene, "_", cond1), cond2_temp=paste0(Gene, "_", cond2)) %>%
      data.frame()
    # add mean count to cond1 and cond2
    de_df$cond1_mean <- unname(sapply(de_df$cond1_temp, function(x){count_mean[x]}))
    de_df$cond2_mean <- unname(sapply(de_df$cond2_temp, function(x){count_mean[x]}))
    
    de_df <- de_df %>%
      arrange(-log2FoldChange, padj) %>%
      dplyr::select(Gene, Comparison, log2FoldChange, padj, cond1_mean, cond2_mean) %>%
      mutate(log2FoldChange=round(log2FoldChange, digits=2), padj=format(padj, scientific=TRUE, digits=3), cond1_mean=round(cond1_mean, digits=0), cond2_mean=round(cond2_mean, digits=0)) %>% 
      dplyr::rename(`Gene`= Gene, `LogFC`= log2FoldChange, `Q-value`= padj, `Count.mean1` = cond1_mean, `Count.mean2` = cond2_mean)
    row.names(de_df) <- NULL
    de_df
    
    # sras %>% filter(SRA_ID == input$tissue_te) %$% de[[SRA_ID]] %>% subset(gene_symbol %in% curr_gene_te()) %>%  
    #   arrange(-log2FoldChange, padj) %>%   
    #   mutate(log2FoldChange=round(log2FoldChange, digits=2), padj=format(padj, scientific=TRUE, digits=3)) %>% 
    #   dplyr::select(Gene, Comparison,log2FoldChange, padj) %>% #rearrange columns in desired order
    #   dplyr::rename(`Gene`= Gene, `LogFC`= log2FoldChange, `Q-value`= padj)
    #removed p-values from the table, since we had originally saved two of the datasets without p-values to save space
    
  }, options=list(paging=FALSE, searching=FALSE)
  )
  
  ######################
  ## Download buttons ##
  ######################
  
  graphgene=reactive({curr_gene()})
  px_to_inch = 0.0104166653543
  
  #Plot height for downloading plot
  plotHeight <- function(dat){
    if (nrow(dat) <= 3){
      height = 3
    }else if (nrow(dat) %in% seq(4,10)){
      height = 5
    } else {
      height = getHeight_func(dat)*px_to_inch
    }
    return(height)
  }
  
  output$asthma_fc_download <- downloadHandler(
    filename= function(){paste0("REALGAR_asthma_forestplot_", graphgene(), ".png")},
    content=function(file){
      png(file, width=16, height=plotHeight(data3_Asthma()), units="in", res=300)
      print(forestplot_func(data3_Asthma(),"Asthma Transcriptomic Results for ",curr_gene()))
      dev.off()})
  
  output$GC_fc_download <- downloadHandler(
    filename= function(){paste0("REALGAR_treatment_forestplot_", graphgene(), ".png")},
    content=function(file){
      png(file, width=16, height=plotHeight(data3_GC()), units="in", res=300)
      print(forestplot_func(data3_GC(),"Exposure Transcriptomic Results for ",curr_gene()))
      dev.off()})
  
  output$cig_fc_download <- downloadHandler(
    filename= function(){paste0("REALGAR_smoking_forestplot_", graphgene(), ".png")},
    content=function(file){
      png(file, width=16, height=plotHeight(data3_cig()), units="in", res=300)
      print(forestplot_func(data3_cig(),"Smoking Transcriptomic Results for ",curr_gene()))
      dev.off()})
  
  output$gene_tracks_download <- downloadHandler(
      filename= function(){paste0("REALGAR_gene_tracks_", graphgene(), ".png")},
      content=function(file){
          png(file, width=16, height=12, units="in", res=300)
          print(make_karyoplot(gene.region(),snp_subs(), snp_gabriel_subs(), snp_fer_subs(),
                               snp_eve_all_subs(), snp_eve_ea_subs(), snp_eve_aa_subs(), snp_eve_la_subs(), 
                               snp_TAGC_multi_subs(), snp_TAGC_euro_subs(),
                               snp_UKBB_asthma_subs(), snp_UKBB_copd_subs(), snp_UKBB_aco_subs()))
          dev.off()})
  
  output$table_download <- downloadHandler(filename = function() {paste0('REALGAR_expression_summary_table_',graphgene(), '.csv')},
                                           content = function(file) {write.csv(do.call(rbind, list(tableforgraph_Asthma(), tableforgraph_GC(), tableforgraph_cig())), file, row.names=FALSE)})
  
  output$SNP_data_download <- downloadHandler(filename = function() {paste0('REALGAR_SNP_results_',graphgene(), '.csv')},
                                              content = function(file) {write.csv(snp_data(), file, row.names=FALSE)})
  output$GRbinding_data_download <- downloadHandler(filename = function() {paste0('REALGAR_GR_binding_site_results_',graphgene(), '.csv')},
                                                    content = function(file) {write.csv(GRbinding_table(), file, row.names=FALSE)})
  
  #download gene boxplot
  output$downloadPic <- downloadHandler(
    filename = function() {paste(input$tissue_te, "_", curr_gene_te(), "_", Sys.Date(), '.png', sep='')},
    content = function(file) {
      png(file, width=10, height=6, units="in", res=600)
      print(getGeneBoxPlot())
      dev.off()
    },
    contentType = 'image/png'
  )
  
})
