# detach("package:Gviz", unload=TRUE) # this is to keep RStudio happy - run if loading app more than once in same session - keep commented out otherwise
                                    # if load Gviz 2x in same session (i.e. close & re-run app), get "object of type 'closure' is not subsettable" error
                                    # should not be an issue when running app from the website
# cat(file=stderr(), as.character(Sys.time()),"packages start\n")
                                    # use this type command to easily see dataset loading time in RStudio  
                                    # currently 3 seconds from "start package load" to "finish gene_locations load"
library(shiny)
library(dplyr)
library(data.table)
library(forestplot) 
library(lattice)
library(stringr)
library(RColorBrewer)
library(viridis) 
library(DT) 
library(Gviz) 

# load dataset descriptions
Dataset_Info <- readRDS("databases/microarray_data_infosheet_R.RDS")
Dataset_Info$Unique_ID <- apply(Dataset_Info[, c("GEO_ID", "Tissue", "Asthma")], 1, paste, collapse = "_")

#load and name datasets
for (i in Dataset_Info$Unique_ID) {assign(i, readRDS(paste0("databases/microarray_results/", i, ".RDS")))}
Dataset_Info[is.na(Dataset_Info$PMID), "PMID"] <- ""

tfbs <-readRDS("databases/tfbs_for_app.RDS") #TFBS data from ENCODE - matched to gene ids using bedtools
snp <- readRDS("databases/grasp_output_for_app.RDS") #SNP data from GRASP - matched to gene ids using bedtools
snp_eve <- readRDS("databases/eve_data_realgar.RDS")#SNP data from EVE - matched to gene ids using bedtools #still have to do liftover
gene_locations <- fread("databases/gene_positions.txt", header = TRUE, stringsAsFactors = FALSE) #gene location & transcript data from GENCODE
chrom_bands <- readRDS("databases/chrom_bands.RDS") #chromosome band info for ideogram - makes ideogram load 25 seconds faster
#unlike all other files, gene_locations is faster with fread than with readRDS (2s load, vs 4s)

#color tfbs based on binding score - used in tracks
#create color scheme based on values tfbs binding score & snp p-values
getPalette = colorRampPalette(brewer.pal(9, "Blues"))
tfbs$color <- getPalette(50)[as.numeric(cut(tfbs$score,breaks = 50))]
snp$color <- getPalette(1024)[as.numeric(cut(-snp$p,breaks = 1024))]
snp_eve$color <- getPalette(1024)[as.numeric(cut(-snp_eve$meta_P,breaks = 1024))]

snp_eve$start <- snp_eve$end

# make a list of gene symbols in all datasets for checking whether gene symbol entered is valid - used for GeneSymbol later on
genes_avail <- vector()
for (i in ls()[grep("GSE", ls())]) {genes_avail <- unique(c(genes_avail, get(i)$SYMBOL))}

output.table <- data.frame() # initiate output table - used later in output.tableforplot()
heatmap_colors <-  inferno # heatmap colors - used in p-value plot

# server
shinyServer(function(input, output, session) {
  
    curr_gene <- reactive({
      if (tolower(input$curr_gene) %in% snp$snp) { #if SNP ID is entered, convert internally to corresponding gene symbol  
          snp$symbol[which(snp$snp==input$curr_gene)]
          } else {
          gsub(" ", "", toupper(input$curr_gene), fixed = TRUE) #make uppercase, remove spaces
          }
      }) 
    
  GeneSymbol <- reactive({if (curr_gene() %in% genes_avail) {TRUE} else {FALSE}})  #used later to generate error message when a wrong gene symbol is input
  
  ##############################################
  ## "Select all" button for tissue selection ##
  ##############################################
  checkbox_choices <- c("Airway smooth muscle"="ASM", "Bronchial epithelium"="BE", 
                        "Bronchoalveolar lavage"="BAL", "CD4"="CD4", "CD8"="CD8",
                        "Lens epithelium" = "LEC","Lymphoblastic leukemia cell" = "chALL", 
                        "Lymphoblastoid cell" = "LCL","Macrophage" = "MACRO", "MCF10A-Myc" = "MCF10A-Myc",
                        "Nasal epithelium"="NE","Osteosarcoma U2OS cell" = "U2O", 
                        "Peripheral blood mononuclear cell"="PBMC","White blood cell"="WBC","Whole lung"="Lung")

  observe({
      if(input$selectall == 0) return(NULL) 
      else if (input$selectall%%2 == 0)
      {updateCheckboxGroupInput(session,"Tissue","Tissue",choices=checkbox_choices, inline = TRUE)}
      else
      {updateCheckboxGroupInput(session,"Tissue","Tissue",choices=checkbox_choices,selected=c("BE", "LEC", "NE", "CD4", "CD8", "PBMC", "WBC", "ASM", "BAL", "Lung",
                                                                                                 "chALL", "MCF10A-Myc", "MACRO", "U2O", "LCL"), inline = TRUE)}})
  #######################
  ## GEO studies table ##
  #######################
  #select GEO studies matching desired conditions;
  #Jessica's initial app had an "and" condition here; I changed it to "or"
  UserDataset_Info <- reactive({
      Dataset_Info1 = subset(Dataset_Info,(((Dataset_Info$Tissue %in% input$Tissue) | (Dataset_Info$Asthma %in% input$Asthma)) & Dataset_Info$App == "asthma")) 
      Dataset_Info2 = subset(Dataset_Info, ((Dataset_Info$Tissue %in% input$Tissue) & (Dataset_Info$App %in% input$GC_included))) 
      Dataset_Info = rbind(Dataset_Info1, Dataset_Info2)}) # this separates GC and asthma data 
  
  #add links for GEO_ID and PMID
  GEO_data <- reactive({
      validate(need(nrow(UserDataset_Info()) != 0, "Please choose at least one dataset.")) #Generate a error message when no data is loaded.
     
       UserDataset_Info() %>%
          dplyr::mutate(GEO_ID_link = paste0("http://www.ncbi.nlm.nih.gov/gquery/?term=", GEO_ID),
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
      x %>%
          dplyr::filter(SYMBOL==curr_gene()) %>%
          dplyr::select(logFC, P.Value, adj.P.Val,t) %>% 
          dplyr::filter(P.Value==min(P.Value)) %>%
          dplyr::mutate(lower = logFC - 2* (logFC/t), upper = logFC + 2*(logFC/t))}
     
  #get data for given gene for each study selected
  for (i in UserDataset_Info()$Unique_ID){
      curr.gene.data=get(i,environment()) 
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
            curr.gene.data <- cbind(data_type, Unique_ID=i, curr.gene.data)
            #append curr.gene.data to all the other data that needs to be included in the levelplots
            output.table <- rbind(output.table, curr.gene.data)}}}
  
    #preparing the data for levelplots
    #calculate the fold change, order by fold change for levelplots
  validate(need(GeneSymbol() != FALSE, "Please enter a valid gene symbol or SNP ID.")) # Generate error message if the gene symbol is not right.
  output.table <- dplyr::mutate(output.table, Fold_Change=2^(logFC), neglogofP=(-log10(adj.P.Val)), Lower_bound_CI = 2^(lower), Upper_bound_CI = 2^(upper)) #note that this is taking -log10 of adjusted p-value
  row.names(output.table) <- output.table$Unique_ID #crucial for plot labels on levelplot
  output.table <- output.table[order(output.table$Fold_Change),]})
  

  ###################################
  ## Data table accompanying plots ##
  ###################################
  
  # asthma
  data_Asthma <- reactive({ output.tableforplot_asthma = output.tableforplot() 
  output.tableforplot_asthma = output.tableforplot_asthma[output.tableforplot_asthma$App == "asthma",]
  output.tableforplot_asthma[rev(rownames(output.tableforplot_asthma)),]})
  
  data2_Asthma <- reactive({
      data_Asthma()%>%
          dplyr::select(Unique_ID, adj.P.Val, P.Value,Fold_Change, neglogofP, Lower_bound_CI, Upper_bound_CI) %>%
          dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), P.Value =format(P.Value, scientific=TRUE, digits=3), 
                        Lower_bound_CI = round(Lower_bound_CI, digits = 2), Upper_bound_CI = round(Upper_bound_CI, digits = 2), Comparison = "Asthma vs. non-asthma")%>%
          dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Fold Change`=Fold_Change)})
  
  tableforgraph_Asthma <- reactive(data2_Asthma()%>% 
                                       dplyr::mutate(`Fold Change(95% CI)` = paste(`Fold Change`, " ","(", Lower_bound_CI, ",", Upper_bound_CI, ")", sep = "")) %>%
                                       dplyr::select(`Study ID`, `Comparison`, `P Value`, `Q Value`, `Fold Change(95% CI)`))
  # GC
  data_GC <- reactive({ output.tableforplot_GC = output.tableforplot()
  output.tableforplot_GC = output.tableforplot_GC[output.tableforplot_GC$App == "GC",]
  output.tableforplot_GC[rev(rownames(output.tableforplot_GC)),]})
  
  data2_GC <- reactive({
      data_GC()%>%
          dplyr::select(Unique_ID, adj.P.Val, P.Value, Fold_Change, neglogofP, Lower_bound_CI, Upper_bound_CI) %>%
          dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), P.Value =format(P.Value, scientific=TRUE, digits=3), 
                 Lower_bound_CI = round(Lower_bound_CI, digits = 2), Upper_bound_CI = round(Upper_bound_CI, digits = 2), Comparison = "Glucocorticoid vs. control")%>%
          dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Fold Change`=Fold_Change)})
  
  tableforgraph_GC <- reactive(data2_GC()%>% 
                                   dplyr::mutate(`Fold Change(95% CI)` = paste(`Fold Change`, " ","(", Lower_bound_CI, ",", Upper_bound_CI, ")", sep = "")) %>%
                                   dplyr::select(`Study ID`, `Comparison`, `P Value`, `Q Value`, `Fold Change(95% CI)`))
  #combine asthma & GC into one
  output$tableforgraph <- DT::renderDataTable(rbind(tableforgraph_Asthma(), tableforgraph_GC()),
                                                 class = 'cell-border stripe', 
                                                 rownames = FALSE, 
                                                 options = list(paging = FALSE, searching = FALSE),
                                                 width = "100%")
  #################
  ## Forestplots ##
  #################
  
  #asthma forestplot
  forestplot_asthma <- function() {
      
      data2_Asthma = data2_Asthma()
          
      validate(need(nrow(data2_Asthma) != 0, "Please choose a dataset.")) 
      
      # function to color forestplot lines and boxes by -log10 of adjusted pvalue - always relative to the max of 8
      color_fn <- local({
          i <- 0
          breaks <- c(seq(0,8,by=0.001), Inf) # this sets max universally at 8 (else highest one OF THE SUBSET would be the max)
          b_clrs <- l_clrs <- inferno(8002)[as.numeric(cut(data2_Asthma$neglogofP, breaks = breaks))] #8002 is length(breaks) - ensures there are enough colors

          function(..., clr.line, clr.marker){
              i <<- i + 1
              fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
          }
      })
      
      text_asthma = data2_Asthma$`Study ID`
      
      xticks = seq(from = min(0.9, min(data2_Asthma$Lower_bound_CI)), to = max(max(data2_Asthma$Upper_bound_CI),1.2), length.out = 5)
      
      
      tabletext <- matrix(nrow=1, ncol=4)
      tabletext[1,] <- c("GEO ID", "Tissue", "Endotype", "Q Value")
      text_temp <- cbind(as.matrix(Dataset_Info[which(Dataset_Info$Unique_ID %in% data2_Asthma$`Study ID`),c(1,11,4)]), data2_Asthma$`Q Value`)
      tabletext <- rbind(tabletext,text_temp)

      
      forestplot(tabletext, title = "Asthma vs. Non-asthma", rbind(c(NA,NA,NA,NA),data2_Asthma[,c("Fold Change","Lower_bound_CI","Upper_bound_CI")]), zero = 1, 
                 xlab = "Fold Change",boxsize = 0.2, col = fpColors(lines = "navyblue", box = "royalblue", zero = "lightgrey"), lwd.ci = 2, 
                 xticks = xticks, is.summary=c(TRUE,rep(FALSE,nrow(data2_Asthma))), lineheight = unit((19/nrow(data2_Asthma)), "cm"),mar = unit(c(5,0,0,5),"mm"), fn.ci_norm = color_fn,
                 txt_gp = fpTxtGp(cex = 1.2, xlab = gpar(cex = 1.35), ticks = gpar(cex = 1.2), title = gpar(cex = 1.45)))
  }
  
  #GC forestplot
  forestplot_GC <- function(){
      
      data2_GC = data2_GC()
      validate(need(nrow(data2_GC) != 0, "Please choose a dataset."))
      
      # function to color forestplot lines and boxes by -log10 of adjusted pvalue - always relative to the max of 8
      color_fn <- local({
          i <- 0
          breaks <- c(seq(0,8,by=0.001), Inf) # this sets max universally at 8 (else highest one OF THE SUBSET would be the max)
          b_clrs <- l_clrs <- inferno(8002)[as.numeric(cut(data2_GC$neglogofP, breaks = breaks))] # 8002 is length(breaks) - ensures there are enough colors
          
          function(..., clr.line, clr.marker){
              i <<- i + 1
              fpDrawNormalCI(..., clr.line = l_clrs[i], clr.marker = b_clrs[i])
          }
      })
      
      text_GC = data2_GC$`Study ID`
      
      xticks = seq(from = min(min(0.9, data2_GC$Lower_bound_CI)), to = max(max(data2_GC$Upper_bound_CI),1.2), length.out = 5)
      
      tabletext <- matrix(nrow=1, ncol=4)
      tabletext[1,] <- c("GEO ID", "Tissue", "Treatment", "Q Value")
      text_temp <- cbind(as.matrix(Dataset_Info[which(Dataset_Info$Unique_ID %in% data2_GC$`Study ID`),c(1,11,6)]), data2_GC$`Q Value`)
      tabletext <- rbind(tabletext,text_temp)
      
      par(bg = "thistle")
      forestplot(tabletext, title = "Glucocorticoid vs. Control", rbind(c(NA,NA,NA,NA),data2_GC[,c("Fold Change","Lower_bound_CI","Upper_bound_CI")]) ,zero = 1, 
                 xlab = "Fold Change",boxsize = 0.15, col = fpColors(lines = "navyblue", box = "royalblue", zero = "lightgrey"), lwd.ci = 2,
                 xticks = xticks, is.summary=c(TRUE,rep(FALSE,nrow(data2_GC))), lineheight = unit((15/(nrow(data2_GC))), "cm"),mar = unit(c(5,0,0,10),"mm"), fn.ci_norm = color_fn,
                 txt_gp = fpTxtGp(cex = 1.2, xlab = gpar(cex = 1.35), ticks = gpar(cex = 1.2), title = gpar(cex = 1.45)))}
  
  output$forestplot_asthma = renderPlot({forestplot_asthma()}, height=650)
  output$forestplot_GC = renderPlot({forestplot_GC()}, height=650)
  
  output$color_scale1 <- output$color_scale2 <- renderImage({ #need two separate output names - else it fails (can't output same thing twice?)
      return(list(
              src = "databases/www/color_scale.png",
              height=550,
              width=59,
              filetype = "image/png",
              alt = "color_scale"
          ))

  }, deleteFile = FALSE)
  
  # output$color_scale2 <- renderImage({ #need two separate output names - else it fails (can't output same thing twice?). also different size b/c different number of asthma & GC datasets
  #     return(list(
  #         src = "databases/www/color_scale.png",
  #         height=430,
  #         width=46,
  #         filetype = "image/png",
  #         alt = "color_scale"
  #     ))
  #     
  # }, deleteFile = FALSE)
  ###############################
  ## Gene, SNP and TFBS tracks ##
  ###############################
  
  #filter data for selected gene
  gene_subs <- reactive({
      gene_subs_temp <- unique(filter(gene_locations, symbol==curr_gene()))
      gene_subs_temp <- gene_subs_temp[!(duplicated(gene_subs_temp$exon)),]})
  tfbs_subs <- reactive({unique(filter(tfbs, symbol==curr_gene()))})
  snp_subs <- reactive({unique(filter(snp, symbol==curr_gene()))})
  snp_eve_subs <- reactive({unique(filter(snp_eve, symbol==curr_gene()))})
  
    gene_tracks <- function() {
      validate(need(curr_gene() != "", "Please enter a gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
      validate(need(GeneSymbol() != FALSE, "Please enter a valid gene symbol or SNP ID.")) # Generate error message if the gene symbol is not right.
      validate(need(nrow(UserDataset_Info()) != 0, "Please choose at least one dataset.")) #Generate a error message when no data is loaded.

      gene_subs <- gene_subs()
      tfbs_subs <- tfbs_subs()
      snp_subs <- snp_subs()
      snp_eve_subs <- snp_eve_subs()

      #constant for all tracks
      gen <- "hg19"
      chr <- unique(gene_subs$chromosome)

      #chromosome, axis and gene - these tracks show up for all genes
      bands <- chrom_bands[which(chrom_bands$chrom==chr),]
      chrom_track <- IdeogramTrack(genome = gen, bands = bands) # formerly slow b/c of chromosome=chr; see https://support.bioconductor.org/p/78881/
      axis_track <- GenomeAxisTrack()
      gene_track <- Gviz::GeneRegionTrack(gene_subs, genome = gen, chromosome = chr, name = "Transcripts", transcriptAnnotation="transcript", fill = "royalblue")
      
      #tfbs and snp tracks - only present for some genes
      
      #TFBS
      if (nrow(tfbs_subs) > 0) {tfbs_track <- Gviz::AnnotationTrack(tfbs_subs, name="GR binding", fill = tfbs_subs$color, group = " ")}
      
      # GRASP SNPs track
      if (nrow(snp_subs) > 0) { 
          snp_track <- Gviz::AnnotationTrack(snp_subs, name="SNPs (from GRASP)", fill = snp_subs$color, group=snp_subs$snp)
          
          #rough estimate of number of stacks there will be in SNP track - for track scaling
          if (nrow(snp_subs) > 1) {
              snp_subs_temp <- snp_subs
              snp_range <- max(snp_subs_temp$start) - min(snp_subs_temp$start)
              snp_subs_temp$start_prev <- c(0, snp_subs_temp$start[1:(nrow(snp_subs_temp)-1)])
              snp_subs_temp$dist <- as.numeric(snp_subs_temp$start) - as.numeric(snp_subs_temp$start_prev)
              snp_size_init <- 2 + as.numeric(nrow(snp_subs[which(snp_subs$dist < snp_range/10),])/2) + 0.8*length(unique(gene_subs$transcript))
          } else {snp_size_init <- 1.2 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs)}
      } else {snp_size_init <- 1.2 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs)}

      # EVE SNPs track
      if (nrow(snp_eve_subs) > 0) {
          snp_eve_track <- Gviz::AnnotationTrack(snp_eve_subs, name="SNPs (from EVE)", fill = snp_eve_subs$color, group=snp_eve_subs$snp_id)

          #rough estimate of number of stacks there will be in SNP track - for track scaling
          if (nrow(snp_eve_subs) > 1) {
              snp_eve_subs_temp <- snp_eve_subs
              snp_eve_range <- max(snp_eve_subs_temp$start) - min(snp_eve_subs_temp$start)
              snp_eve_subs_temp$start_prev <- c(0, snp_eve_subs_temp$start[1:(nrow(snp_eve_subs_temp)-1)])
              snp_eve_subs_temp$dist <- as.numeric(snp_eve_subs_temp$start) - as.numeric(snp_eve_subs_temp$start_prev)
              snp_eve_size_init <- 2 + as.numeric(nrow(snp_eve_subs[which(snp_eve_subs$dist < snp_eve_range/10),])/2) + 0.8*length(unique(gene_subs$transcript))
          } else {snp_eve_size_init <- 1.2 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_eve_subs)}
      } else {snp_eve_size_init <- 1.2 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_eve_subs)}

      #track sizes - defaults throw off scaling as more tracks are added
       chrom_size <- 1.2 + 0.01*length(unique(gene_subs$transcript)) + 0.01*nrow(snp_subs)
       axis_size <- 1 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs)
       gene_size <- 2 + 0.6*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs)
       tfbs_size <- 2 + 0.075*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs)
       snp_size <- snp_size_init #from above
       snp_eve_size <- snp_eve_size_init #from above

      #select the non-empty tracks to output -- output depends on whether there are TFBS and/or SNPs for a given gene
       subset_size <- sapply(c("tfbs_subs", "snp_subs", "snp_eve_subs"), function(x) {nrow(get(x))}) #size of each subset
       non_zeros <- names(subset_size)[which(!(subset_size==0))] #which subsets have non-zero size
       
       df_extract <- function(x,y) { #gives name of track and track size variable for non-zero subsets (y is "track" or "size")
           if (length(non_zeros) > 0) {
               get(paste0(strsplit(x, 'subs'),y)) #trim off "subs" and append either "track" or "size"
           } else {NULL} #to avoid meltdown if no subsets were non-zero
       }
       
       #use df_extract function to get track & track size corresponding to all non-zero subsets
       #note chrom_track, axis_track and gene_track are present for all
       selected_tracks <- list(chrom_track, axis_track, gene_track, sapply(non_zeros, df_extract, y="track")$tfbs_subs, sapply(non_zeros, df_extract, y="track")$snp_subs, sapply(non_zeros, df_extract, y="track")$snp_eve_subs)
       selected_tracks <- Filter(Negate(function(x) is.null(unlist(x))), selected_tracks) #remove null elements from list
       
       selected_sizes <- na.omit(c(chrom_size,axis_size,gene_size, sapply(non_zeros, df_extract, y="size")[1], sapply(non_zeros, df_extract, y="size")[2], sapply(non_zeros, df_extract, y="size")[3]))
       #note: use names to extract from selected_tracks b/c it is a list vs. index to extract from selected_sizes, since this is numeric
       
       #plot tracks 
       plotTracks(selected_tracks, sizes=selected_sizes, col=NULL, background.panel = "#d3cecc", background.title = "firebrick4", col.border.title = "firebrick4", groupAnnotation = "group", fontcolor.group = "darkblue", cex.group=0.75, just.group="below", cex.title=1.1)
       
    }
    
     #plot height increases if more tracks are displayed
     observe({output$gene_tracks_outp2 <- renderPlot({gene_tracks()}, height=400 + 15*length(unique(gene_subs()$transcript)) + 10*(nrow(snp_subs())), width=1055)})

  ######################
  ## Download buttons ##
  ######################
  graphgene=reactive({curr_gene()})
    
  output$asthma_fc_download <- downloadHandler(
    filename= function(){paste0("fold_change_asthma_", graphgene(), "_", Sys.Date(), ".png")},
    content=function(file){
      png(file, width=6, height=9, units="in", res=600)
      forestplot_asthma()
      dev.off()})
  
  output$GC_fc_download <- downloadHandler(
      filename= function(){paste0("fold_change_GC_", graphgene(), "_", Sys.Date(), ".png")},
      content=function(file){
          png(file, width=6, height=9, units="in", res=600)
          forestplot_GC()
          dev.off()})
  
  output$pval_download <- downloadHandler(
      filename= function(){paste0("-log(pval)_heatmap_", graphgene(), "_", Sys.Date(), ".png")},
      content=function(file){
          png(file, width=6, height=9, units="in", res=600)
          print(pval_plot()) # note that for this one, unlike other plot downloads, had to use print(). 
          dev.off()})        # else the download is a blank file. this seems to be b/c pval_plot() creates a graph 
                             # object but doesn't draw the plot, as per 
                             # http://stackoverflow.com/questions/27008434/downloading-png-from-shiny-r-pt-2
       
  output$gene_tracks_download <- downloadHandler(
      filename= function(){paste0("gene_tracks_", graphgene(), "_", Sys.Date(), ".png")},
      content=function(file){
          png(file, width=16, height=12, units="in", res=600)
          gene_tracks()
          dev.off()})
  
  output$table_download <- downloadHandler(filename = function() {paste0('Asthma&GC_data_summary_table_',graphgene(), Sys.Date(), '.csv')},
                                                  content = function(file) {write.csv(rbind(tableforgraph_Asthma(), tableforgraph_GC()), file, row.names=FALSE)})
})
