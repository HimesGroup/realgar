library(shiny)
library(dplyr)
library(data.table)
library(lattice)
library(stringr)
library(RColorBrewer)
library(DT)
library(GenomicRanges)
library(Gviz) # Gviz causes an odd problem... if run app second time within same session, get an "object of type 'closure' is not subsettable" error
              # this comes from loading Gviz, not from any code associated with it - and running first time is fine
              # hopefully will not affect the app when it is online

# Load data
Dataset_Info <- read.csv("databases/microarray_data_infosheet.csv")
Dataset_Info$Unique_ID <- apply(Dataset_Info[, c("GEO_ID", "Tissue", "Asthma")], 1, paste, collapse="_")
for (i in Dataset_Info$Unique_ID){
    assign(i, fread(paste0("databases/microarray_results/", i,".csv"), sep=" "))}
output.table <- data.frame()
Dataset_Info[is.na(Dataset_Info$PMID),"PMID"] <- ""
tfbs <- read.table("databases/tfbs_for_app.txt", header = TRUE, stringsAsFactors = FALSE) #TFBS data from ENCODE - matched to gene ids using bedtools
snp <- read.table("C:/Users/mayashum/Documents/realgar/databases/grasp_output_for_app.txt", header = TRUE, stringsAsFactors = FALSE) #SNP data from GRASP - matched to gene ids using bedtools

# heatmap colors
heatmap_colors <- colorRampPalette(c("navyblue","darkgoldenrod1","firebrick4"))

# Server
shinyServer(function(input,output) {
  
  curr_gene <- reactive({toupper(input$curr_gene)}) #can recognize gene names even if typed lowercase
  
  #######################
  ## GEO studies table ##
  #######################
  #select GEO studies matching desired conditions;
  #Jessica's initial app had an "and" condition here; I changed it to "or"
  UserDataset_Info <- reactive({
    Dataset_Info=subset(Dataset_Info,Dataset_Info$Tissue %in% input$Tissue | Dataset_Info$Asthma %in% input$Asthma) 
    Dataset_Info})
  
  #add links for GEO_ID and PMID
  GEO_data <- reactive({UserDataset_Info() %>%
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
  

  #select and modify data used for levelplots and accompanying table
  output.tableforplot <- reactive({
    
    #select data for the gene currently selected
    data_filter <- function(x){
      x %>%
        dplyr::filter(SYMBOL==curr_gene()) %>%
        dplyr::select(logFC, P.Value, adj.P.Val) %>%
        dplyr::filter(P.Value==min(P.Value))}
    
    #get data for given gene for each study selected
    for (i in UserDataset_Info()$Unique_ID){
      curr.gene.data=get(i,environment())
      
      if(any(tbl_vars(curr.gene.data)=="qValuesBatch")) {
        curr.gene.data <- (curr.gene.data %>%
                             dplyr::select(-P.Value,-adj.P.Val) %>%
                             dplyr::rename(P.Value=pValuesBatch) %>%
                             dplyr::rename(adj.P.Val=qValuesBatch))}
      
      #use data_filter function from above to filter curr.gene.data 
      curr.gene.data <- data_filter(curr.gene.data)
      if(nrow(curr.gene.data) > 0) {
        curr.gene.data <- cbind(Unique_ID=i, curr.gene.data)
        #append curr.gene.data to all the other data that needs to be included in the levelplots
        output.table <- rbind(output.table, curr.gene.data)}}
    
    #preparing the data for levelplots
    #calculate the Fold Change, order by Fold Change for levelplots
    output.table <- dplyr::mutate(output.table, Fold_Change=2^(logFC), neglogofP=(-log10(adj.P.Val))) #note that this is taking -log10 of adjusted p-value
    row.names(output.table) <- output.table$Unique_ID #crucial for plot labels on levelplot
    output.table <- output.table[order(output.table$Fold_Change),]})
  
  ########################################
  ## Data table accompanying levelplots ##
  ########################################
  levelplot_data <- reactive({dplyr::select(output.tableforplot(), Unique_ID, P.Value,Fold_Change,neglogofP)
    output.tableforplot()[rev(rownames(output.tableforplot())),]})
  
  levelplot_data_outp <- reactive({levelplot_data()%>%
      dplyr::select(Unique_ID,Fold_Change,adj.P.Val,P.Value)%>%
      dplyr::mutate(Fold_Change=round(Fold_Change,digits=2),adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), P.Value =format(P.Value, scientific=TRUE, digits=3))%>%
      dplyr::rename(`Study ID`=Unique_ID, `P Value`=P.Value, `Q Value`=adj.P.Val, `Log 2 Fold Change`=Fold_Change)})
  
  output$tableforgraph <- DT::renderDataTable(levelplot_data_outp(), class = 'cell-border stripe', rownames = FALSE, options = list(paging = FALSE, searching = FALSE))
  
  ################
  ## Levelplots ##
  ################
  output.tableforplot2 <- reactive({output.tableforplot() %>% dplyr::rename(' '=Fold_Change, ' '=neglogofP)})
  heatmapMAT <- reactive({output.tableforplot2()})
  
  fc_data <- reactive({as.matrix(t(heatmapMAT()[5]))})
  pval_data <- reactive({t(heatmapMAT()[6])})
 
  #set up min & max boundaries for levelplots
  minFC <- reactive({if(max(fc_data())<=1.5 & min(fc_data())>=-1.5){
    minFC=-1.5} else {minFC=min(fc_data())}})
  
  maxFC <- reactive({if(max(fc_data())<=1.5 & min(fc_data())>=-1.5){
    maxFC=1.5} else {maxFC=max(fc_data())}})
  
  # minNLOP <- reactive({0})
  # maxNLOP <- reactive({if(max(plot_data_pval())<=1.5 & min(plot_data_pval())>=-1.5){
  #   maxNLOP=1.5} else {maxNLOP=max(plot_data_pval())}})
  
  # levelplots output for fold change & log p-value
  fc_plot <- reactive({
      levelplot(fc_data(),
                col.regions=heatmap_colors,
                xlab =NULL,
                ylab="GEO ID",
                main = "Fold Change",
                aspect=2,
                par.settings=list(
                    layout.widths=list(left.padding=-4)),
                scales=list(x=list(cex=1, tck = c(0,0,0,0)),
                            y=list(cex=1, tck = c(1,0,0,0))),
                at=seq(minFC(), maxFC(), length.out=100))})
  
  output$fc_plot_outp <- renderPlot({fc_plot()})

  pval_plot <- reactive({
      levelplot(pval_data(),
                col.regions=heatmap_colors,
                xlab =NULL,
                ylab="GEO ID",
                main = "-log10(adjusted p-value)",
                aspect=2,
                scales=list(x=list(cex=1, tck = c(0,0,0,0)), 
                            y=list(cex=1, tck = c(1,0,0,0))),
                at=seq(0, 8,length.out=100))})
  
  output$pval_plot_outp <- renderPlot({pval_plot()})
  
  ###################
  ## SNPs and TFBS ##
  ###################

  # mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  # gene <- makeGene(id = "ENSG00000095203", type="ensembl_gene_id", biomart = mart)
  # 
  # plusStrand <- makeGeneRegion(chromosome = 19, start = 12050000, end = 12230000, strand = "+", biomart = mart)
  # genomeAxis <- makeGenomeAxis(add53 = TRUE)
  # output$GRASP <- renderPlot({gdPlot(list(plusStrand, genomeAxis))})
  # 

  ##########
  ## TFBS ##
  ##########
  
  tfbs_tracks <- reactive({
      tfbs_subs <- filter(tfbs, GENE3==curr_gene())

      if (nrow(tfbs_subs) > 0) {
          gr <- GRanges(seqnames = tfbs_subs$CHR, ranges = IRanges(start = tfbs_subs$START, end = tfbs_subs$STOP))
          chr <- as.character(unique(seqnames(gr)))
          gen <- "hg19"
          atrack <- Gviz::AnnotationTrack(gr, name="NR3C1 binding sites", stacking="dense")
          itrack <- IdeogramTrack(genome = gen, chromosome = chr) # this step really slows the app down...
          gtrack <- GenomeAxisTrack()
          plotTracks(list(gtrack, itrack, atrack), sizes=c(1,0.5,1.25))
      }
  })

    output$tfbs_plot <- renderPlot({tfbs_tracks()})

  ######################
  ## Download buttons ##
  ######################
  graphgene=reactive({curr_gene()})
    
  output$fc_download <- downloadHandler(
    filename= function(){paste0("Fold_change_heatmap_", graphgene(), "_", Sys.Date(), ".png")},
    content=function(file){
      png(file)
      print(fc_plot())
      dev.off()})
  
  output$pval_download <- downloadHandler(
    filename= function(){paste0("-log(pval)_heatmap_", graphgene(), "_", Sys.Date(), ".png")},
    content=function(file){
      png(file)
      print(pval_plot())
      dev.off()})

  output$table_download <- downloadHandler(filename = function() {paste0('Heatmap_summary_table_',graphgene(), Sys.Date(), '.csv')},
                                         content = function(file) {write.csv(data(), file)})
})