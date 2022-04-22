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
    print(paste0("check gene/SNP ID: ", input_current()))
    current <- toString(input_current())
    current_low <- gsub(" ", "", tolower(current), fixed=TRUE) # convert input to lower case
    current_up <- gsub(" ", "", toupper(current), fixed=TRUE) # convert input to upper case
    if (current_up %in% gene_list) { # if input is a gene within the gene_list
      current <- current_up
    } else if (grepl("^rs", current_low)) {
      if (current_low %in% all_snps) { #if SNP ID is entered, convert internally to nearest gene symbol
        rsid <- current_low
        all_matches <- get_matches(rsid)
        print(all_matches)
        current <- as.character(all_matches$symbol)
        if (nrow(all_matches)>1) {
          current <- join_gene_snp(all_matches) # this function will select one gene that nearest by distance to selected snp
        }
        print(current)
      } else { #input is not a valid SNP ID
        current <- "not_in_snp_db"
      }
    } else {
      # if it is not in the list of snps or genes, it is an invalid gene/snp, user will get an "enter valid gene/snp id" message
      current <- "invalid"
    }
    current
  })

  
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


  
  ChIPSeq_links <- reactive({
    chipseq_dataset %>%
      dplyr::mutate(Study_link1=ifelse(grepl("^GSE", Dataset), paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", Dataset), ifelse(grepl("^SRP", Dataset), paste0("https://www.ncbi.nlm.nih.gov/sra/?term=", Dataset), NA))) %>%
      dplyr::mutate(Study_link=paste0("<a href='",  Study_link1, "' target='_blank'>", Dataset,"</a>")) %>%
      dplyr::mutate(PMID_link1=paste0("https://www.ncbi.nlm.nih.gov/pubmed/?term=", PMID)) %>%
      dplyr::mutate(PMID_link=paste0("<a href='",  PMID_link1, "' target='_blank'>", PMID,"</a>")) %>%
      dplyr::mutate(QC_link=paste0(Dataset,"_QC_ChIPSeqReport.html")) %>%
      dplyr::mutate(Report=paste0("<a href='",  QC_link, "' target='_blank'>", "QC","</a>")) %>%
      dplyr::select(Study_link, PMID_link, Report, Description) %>%
      dplyr::rename(Dataset=Study_link, PMID=PMID_link)
    
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

  output$chipseq_table <- DT::renderDataTable(ChIPSeq_links(),
                                           class = 'cell-border stripe',
                                           rownames = FALSE,
                                           options = list(paging = FALSE, searching = FALSE), 
                                           escape = FALSE)
  
  
    
  #########################################
  ## Select GEO data for plots and table ##
  #########################################
  
  #select and modify data used for plots and accompanying table

  check_input_message <- reactive({ # check input
    # move the check process out of the function output.tableforplot and add it to forest plot output chunck, otherwise, an error message will display
    print("generate check message")
    message <- ""
    if (nrow(UserDataset_Info()) == 0) {
      message <- "Please choose at least one dataset." #Generate a error message when no data is loaded.
    }
    if (curr_gene() == "") {
      message <- "Please enter a gene symbol or SNP ID." #Generate a error message when no gene id is input.
    }
    if (curr_gene() == "invalid") {
      message <- "Please enter a valid gene symbol or SNP ID."
    }
    if (curr_gene() == "not_in_snp_db") {
      message <- "The input SNP ID is not present in the current datasets"
    }
    message
  })  
  
  output.tableforplot <- reactive({
    gene <- paste0('"',curr_gene(),'"')
    query <- paste0("SELECT * FROM REALGAR WHERE Gene = ",gene)
    res <- dbSendQuery(ngs_db, query)
    out.table <- dbFetch(res)
    dbClearResult(res)

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
    print("generate asthma dataset")
    output.tableforplot_asthma <- output.tableforplot() %>% dplyr::filter(App == "asthma")
  })
  
  # GC
  data_GC <- reactive({
    print("generate GC dataset")
    output.tableforplot_GC <- output.tableforplot() %>% dplyr::filter(App == "GC")
  })
  
  # Smoking
  data_cig <- reactive({
    print("generate cig dataset")
    output.tableforplot_cig <- output.tableforplot() %>% dplyr::filter(App == "Smoking")
  })

  
  ###################################
  ##        Combined p-values      ##
  ###################################
  
  # asthma
  asthma_pcomb <- reactive({
    print("combine p-values for asthma datasets")
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
    print("combine p-values for GC datasets")
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
    print("combine p-values for cigarette datasets")
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
    print("meta-analysis for asthma datasets")
    dat <- data_Asthma()
    if (length(dat$adj.P.Val)>1) {
      res <- meta_stat(dat)
    }
    else {res <- list(meta_pval=1,meta_fc=1,meta_lower=1,meta_upper=1)} # assign any values to trick the app if no "treatment" is selected
    res
  })  
  
  meta_GC <- reactive({
    print("meta-analysis for GC datasets")
    dat <- data_GC()
    if (length(dat$adj.P.Val)>1) {
      res <- meta_stat(dat)
    }
    else {res <- list(meta_pval=NA,meta_fc=NA,meta_lower=NA,meta_upper=NA)}
    res
  })  
  
  meta_cig <- reactive({
    print("meta-analysis for cigarette datasets")
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
    print("generate intermediate asthma table")
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
    print("generate final asthma table")
    interdata_func(data2_Asthma(),meta_Asthma())
  }) 
  
  # output table for asthma
  tableforgraph_Asthma <- reactive({
    print("generate asthma forestplot table")
    tableforgraph_func(data3_Asthma())
  })
  
  # GC
  data2_GC <- reactive({
    print("generate intermediate GC table")
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
    print("generate final GC table")
    interdata_func(data2_GC(),meta_GC())
  }) 
  
  # output table for GC
  tableforgraph_GC <- reactive({
    print("generate GC forestplot table")
    tableforgraph_func(data3_GC())
  })
  
  
  # Smoking
  data2_cig <- reactive({
    print("generate intermediate asthmacigarette table")
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
    print("generate final cigarette table")
    interdata_func(data2_cig(),meta_cig())
  }) 
  
  # output table for Smoking
  tableforgraph_cig <- reactive({
    print("generate cigarette forestplot table")
    tableforgraph_func(data3_cig())
  })
  
  
  #combine asthma & GC & Smoking into one
  output$tableforgraph <- DT::renderDataTable(dplyr::bind_rows(tableforgraph_Asthma(), tableforgraph_GC(),tableforgraph_cig()),
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
    validate(need(check_input_message()=="", check_input_message()))
    print("present asthma forest plot")
    forestplot_func(data3_Asthma(),"Asthma Transcriptomic Results for ",curr_gene())
  }) 
  
  # GC forestplot
  
  getHeightGC <- reactive({
    getHeight_func(data3_GC())
  })
  
  forestplot_GC <- reactive({
    print("present GC forest plot")
    validate(need(check_input_message()=="", check_input_message()))
    forestplot_func(data3_GC(),"Exposure Transcriptomic Results for ",curr_gene())
  }) 
  
  # Smoking forestplot
  
  getHeightcig <- reactive({
    getHeight_func(data3_cig())
  })
  
  forestplot_cig <- reactive({
    print("present cigarette forest plot")
    validate(need(check_input_message()=="", check_input_message()))
    forestplot_func(data3_cig(),"Pollutant Transcriptomic Results for ",curr_gene())
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

  # Use this function to extract gene from GWAS tables
  select_gene_from_GWAStb_func <- function(select_gwas_tb, pval_thr, query_table) {
    if (select_gwas_tb %in% input_which_SNPs()) {
      if (pval_thr=="normal") {
        snp_db <- get_query_db(query_table, curr_gene())
      } else if (pval_thr=="nominal") {
        snp_db <- get_query_nominal_db(query_table, curr_gene())
      } else if (pval_thr=="genomewide") {
        snp_db <- get_query_genomewide_db(query_table, curr_gene())
      } 
      if (nrow(snp_db) > 0){
        snp_db %>% dplyr::mutate(start=end) %>% unique()
      } else {data.frame(matrix(nrow = 0, ncol = 0))}
    } else {data.frame(matrix(nrow = 0, ncol = 0))} #only non-zero if corresponding checkbox is selected - but can't have "NULL" - else get "argument is of length zero" error
  }

  # GRASP  
  snp_subs <- reactive({
    print("select GRASP")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_subs", pval_thr = input$pval_thr, query_table = "snp")
  })
  
  # GABRIEL
  snp_gabriel_subs <- reactive({
    print("select gabriel")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_gabriel_subs", pval_thr = input$pval_thr, query_table = "snp_gabriel")
  })
  
  # Ferreira
  snp_fer_subs <- reactive({
    print("select ferreira")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_fer_subs", pval_thr = input$pval_thr, query_table = "snp_fer")
  })
 
  # EVE ALL
  snp_eve_all_subs <- reactive({
    print("select eve all")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_eve_all_subs", pval_thr = input$pval_thr, query_table = "snp_eve_all")
  })

  # EVE EA
  snp_eve_ea_subs <- reactive({
    print("select eve EA")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_eve_ea_subs", pval_thr = input$pval_thr, query_table = "snp_eve_ea")
  })
  
  # EVE AA
  snp_eve_aa_subs <- reactive({
    print("select eve AA")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_eve_aa_subs", pval_thr = input$pval_thr, query_table = "snp_eve_aa")
  })
  

  # EVE LA  
  snp_eve_la_subs <- reactive({
    print("select eve LA")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_eve_la_subs", pval_thr = input$pval_thr, query_table = "snp_eve_la")
  })
  
  # TAGC Multi-ancestry
  snp_TAGC_multi_subs <- reactive({
    print("select tagc multi")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_TAGC_multi_subs", pval_thr = input$pval_thr, query_table = "snp_TAGC_multi")
  })

  # TAGC EURO
  snp_TAGC_euro_subs <- reactive({
    print("select tagc euro")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_TAGC_euro_subs", pval_thr = input$pval_thr, query_table = "snp_TAGC_euro")
  })

  # UKBB asthma
  snp_UKBB_asthma_subs <- reactive({
    print("select ukbb asthma")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_UKBB_asthma_subs", pval_thr = input$pval_thr, query_table = "snp_UKBB_asthma")
  })

  # UKBB COPD
  snp_UKBB_copd_subs <- reactive({
    print("select ukbb copd")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_UKBB_copd_subs", pval_thr = input$pval_thr, query_table = "snp_UKBB_copd")
  })
  
  # UKBB ACO 
  snp_UKBB_aco_subs <- reactive({
    print("select ukbb aco")
    select_gene_from_GWAStb_func(select_gwas_tb = "snp_UKBB_aco_subs", pval_thr = input$pval_thr, query_table = "snp_UKBB_aco")
  })
  
  #plot height increases if more tracks are displayed
  # observe({output$gene_tracks_outp2 <- renderPlot({gene_tracks()}, width=1055,
  #                                                 height=400 + 15*length(unique(gene_subs()$transcript)) + 30*nrow(snp_eve_subs()) + 10*(nrow(snp_subs())+nrow(snp_gabriel_subs())+nrow(snp_fer_subs())+nrow(snp_TAGC_subs())))})
  
  #################################
  ## SNP data table for download ##
  #################################
  
  outp_gwas_tb_func <- function(dat, pmid, source) {
    if (nrow(dat)>0) {
      if ("pmid"%in%names(dat)) {
        dat %>%
          dplyr::rename(position=end) %>%
          dplyr::mutate(source=source) %>%
          dplyr::select(chromosome, SNP, symbol, position, meta_P, pmid, source) %>%
          dplyr::mutate(pmid=as.character(pmid))
      } else {
        dat %>%
          dplyr::rename(position=end) %>%
          dplyr::mutate(source=source, pmid=pmid) %>%
          dplyr::select(chromosome, SNP, symbol, position, meta_P, pmid, source)
      }
    } else {dat}
  }
  
  # GRASP
  snp_subs_temp <- reactive({
    print("create GRASP table")
    outp_gwas_tb_func(dat = snp_subs(), source = "GRASP", pmid = "")
  })
  
  # GABRIEL  
  snp_gabriel_subs_temp <- reactive({
    print("create GABREIL table")
    outp_gwas_tb_func(dat = snp_gabriel_subs(), source = "GABRIEL", pmid="20860503")
  })
  
  # Ferreira
  snp_fer_subs_temp <- reactive({
    print("create Ferreira table")
    outp_gwas_tb_func(dat = snp_fer_subs(), source = "Ferreira", pmid="29083406")
  })
  
  # EVE ALL
  snp_eve_all_subs_temp <- reactive({
    print("create EVE ALL table")
    outp_gwas_tb_func(dat = snp_eve_all_subs(), source = "EVE_ALL", pmid="21804549")
  })

  # EVE EA
  snp_eve_ea_subs_temp <- reactive({
    print("create EVE EA table")
    outp_gwas_tb_func(dat = snp_eve_ea_subs(), source = "EVE_EA", pmid="21804549")
  })  
 
  # EVE AA
  snp_eve_aa_subs_temp <- reactive({
    print("create EVE AA table")
    outp_gwas_tb_func(dat = snp_eve_aa_subs(), source = "EVE_AA", pmid="21804549")
  })
  
  # EVE LA
  snp_eve_la_subs_temp <- reactive({
    print("create EVE LA table")
    outp_gwas_tb_func(dat = snp_eve_la_subs(), source = "EVE_LA", pmid="21804549")
  })
  
  # TAGC Multi
  snp_TAGC_multi_subs_temp <- reactive({
    print("create TAGC multi table")
    outp_gwas_tb_func(dat = snp_TAGC_multi_subs(), source = "TAGC_Multi", pmid = "29273806")
  }) 
  
  # TAGC EUR
  snp_TAGC_euro_subs_temp <- reactive({
    print("create TAGC EUR table")
    outp_gwas_tb_func(dat = snp_TAGC_euro_subs(), source = "TAGC_EUR", pmid = "29273806")
  }) 
  
  # UKBB asthma
  snp_UKBB_asthma_subs_temp <- reactive({
    print("create UKBB asthma table")
    outp_gwas_tb_func(dat = snp_UKBB_asthma_subs(), source = "UKBB_asthma", pmid="NA")
  }) 
  
  # UKBB COPD
  snp_UKBB_copd_subs_temp <- reactive({
    print("create UKBB COPD table")
    outp_gwas_tb_func(dat = snp_UKBB_copd_subs(), source = "UKBB_COPD", pmid="NA")
  })
 
  # UKBB ACO 
  snp_UKBB_aco_subs_temp <- reactive({
    print("create UKBB ACO table")
    outp_gwas_tb_func(dat = snp_UKBB_aco_subs(), source = "UKBB_ACO", pmid="NA")
  })

  snp_data <- reactive({
    print("combine gwas datasets")
    dat <- dplyr::bind_rows(snp_subs_temp(), snp_gabriel_subs_temp(), snp_fer_subs_temp(),
                                         snp_eve_all_subs_temp(), snp_eve_ea_subs_temp(), snp_eve_aa_subs_temp(), snp_eve_la_subs_temp(), 
                                         snp_TAGC_multi_subs_temp(), snp_TAGC_euro_subs_temp(), 
                                         snp_UKBB_asthma_subs_temp(), snp_UKBB_copd_subs_temp(), snp_UKBB_aco_subs_temp())
    get_eqtl_annot_to_GWAS_tb(dat) # annotate eQTL info to GWAS SNPs
    })
  
  #########################
  ## KARYOPLOTR TRACKS ##
  ########################
 
  gene_to_GRanges <- function(gene_tb) {
    if (gene_tb$strand == -1){
      gene.gr <- toGRanges(data.frame(chr=gene_tb$Chromosome, start=gene_tb$Start, end=gene_tb$End+20000, genome="hg38"))
    } else if (gene_tb$strand == 1){
      gene.gr <- toGRanges(data.frame(chr=gene_tb$Chromosome, start=gene_tb$Start-20000, end=gene_tb$End, genome="hg38"))
    }
    return(gene.gr)
  }
   
  gene.region <- reactive({
    validate(need(curr_gene() != "", "Please enter a gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
    region <- all_genes %>% dplyr::filter(symbol == curr_gene()) 
    #region_split <- strsplit(x=gsub("\\s","",region$Pos), split="-|:| ")[[1]]
    if (nrow(region) > 0){
      gene_to_GRanges(unique(region))
    } else {NULL}
  })
  
  ## Make Plot
  karyoplot <- reactive({
    validate(need(!is.null(gene.region()), "Please enter a valid gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
    print(gene.region())
    print("generate karyoplot")
    make_karyoplot(gene.region(), snp_subs(), snp_gabriel_subs(), snp_fer_subs(),
                   snp_eve_all_subs(), snp_eve_ea_subs(), snp_eve_aa_subs(), snp_eve_la_subs(), 
                   snp_TAGC_multi_subs(), snp_TAGC_euro_subs(),
                   snp_UKBB_asthma_subs(), snp_UKBB_copd_subs(), snp_UKBB_aco_subs())
  })
  
  observe({output$karyoPlot <- renderPlot({karyoplot()})})
  
  # GR-binding sites within the selected region
  GRbinding_table <- reactive({
    if (!is.null(gene.region())) {
      # generate GRanges object for GR-binding sites
      GRbinding_gr <- makeGRangesFromDataFrame(GRbinding[,1:3])
      overlap_gr <- findOverlaps(GRbinding_gr, gene.region())
      row_numbers = queryHits(overlap_gr)
      GRbinding[row_numbers, ] %>% data.frame()
    }
  })

  # GRE motifs within the selected region
  GRE_table <- reactive({
    if (!is.null(gene.region())) {
      # find overlaps between GRanges object for GRE motifs and gene region
      overlap_gre <- findOverlaps(gre, gene.region())
      row_numbers = queryHits(overlap_gre)
      gre_tb <- gre[row_numbers, ] %>% data.frame()
      names(gre_tb)[1] <- "chromosome"
      gre_tb[,1:3]
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
    filename= function(){paste0("REALGAR_pollutant_forestplot_", graphgene(), ".png")},
    content=function(file){
      png(file, width=16, height=plotHeight(data3_cig()), units="in", res=300)
      print(forestplot_func(data3_cig(),"Pollutant Transcriptomic Results for ",curr_gene()))
      dev.off()})
  
#  output$gene_tracks_download <- downloadHandler(
#      filename= function(){paste0("REALGAR_gene_tracks_", graphgene(), ".png")},
#      content=function(file){
#          png(file, width=16, height=12, units="in", res=300)
#          print(make_karyoplot(gene.region(),snp_subs(), snp_gabriel_subs(), snp_fer_subs(),
#                               snp_eve_all_subs(), snp_eve_ea_subs(), snp_eve_aa_subs(), snp_eve_la_subs(), 
#                               snp_TAGC_multi_subs(), snp_TAGC_euro_subs(),
#                               snp_UKBB_asthma_subs(), snp_UKBB_copd_subs(), snp_UKBB_aco_subs()))
#          dev.off()})
  
  output$table_download <- downloadHandler(filename = function() {paste0('REALGAR_expression_summary_table_',graphgene(), '.csv')},
                                           content = function(file) {write.csv(do.call(rbind, list(tableforgraph_Asthma(), tableforgraph_GC(), tableforgraph_cig())), file, row.names=FALSE)})
  
  output$SNP_data_download <- downloadHandler(filename = function() {paste0('REALGAR_SNP_results_',graphgene(), '.csv')},
                                              content = function(file) {write.csv(snp_data(), file, row.names=FALSE)})
  output$GRbinding_data_download <- downloadHandler(filename = function() {paste0('REALGAR_GR_binding_site_results_',graphgene(), '.csv')},
                                                    content = function(file) {write.csv(GRbinding_table(), file, row.names=FALSE)})
  output$GRE_data_download <- downloadHandler(filename = function() {paste0('REALGAR_GRE_motif_results_',graphgene(), '.csv')},
                                                    content = function(file) {write.csv(GRE_table(), file, row.names=FALSE)})

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
