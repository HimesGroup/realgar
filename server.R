# detach("package:Gviz", unload=TRUE) # this is to keep RStudio happy - run if loading app more than once in same session - keep commented out otherwise
                                    # if load Gviz 2x in same session (i.e. close & re-run app), get "object of type 'closure' is not subsettable" error
                                    # should not be an issue when running app from the website
                                             
library(shiny)
library(dplyr)
library(data.table)
library(forestplot)
library(lattice)
library(stringr)
library(RColorBrewer)
library(DT)
library(GenomicRanges)
library(Gviz) 

# load dataset descriptions
Dataset_Info <- read.csv("databases/microarray_data_infosheet_R.csv")
Dataset_Info$Unique_ID <- apply(Dataset_Info[, c("GEO_ID", "Tissue", "Asthma")], 1, paste, collapse="_")

#load and name datasets
for (i in Dataset_Info$Unique_ID) {
    assign(i, fread(paste0("databases/microarray_results/", i,".csv"), sep=","))}
Dataset_Info[is.na(Dataset_Info$PMID),"PMID"] <- ""

tfbs <- read.table("databases/tfbs_for_app.txt", header = TRUE, stringsAsFactors = FALSE) #TFBS data from ENCODE - matched to gene ids using bedtools
snp <- read.table("databases/grasp_output_for_app.txt", header = TRUE, stringsAsFactors = FALSE) #SNP data from GRASP - matched to gene ids using bedtools
gene_locations <- read.table("databases/gene_positions.txt", header = TRUE, stringsAsFactors = FALSE) #gene location data from our hg19 gtf annotation file

#color tfbs based on binding score - used in tracks
#create color scheme based on values
getPalette = colorRampPalette(brewer.pal(9, "Blues"))
tfbs$color <- getPalette(50)[as.numeric(cut(tfbs$SCORE,breaks = 50))]

#for SNP annotations
# snp$PMID_link <- paste0("http://www.ncbi.nlm.nih.gov/pubmed/?term=", snp$PMID)
# snp$pval_annot <- paste0(snp$SNP, "\nPMID: \n", "<a href='",  snp$PMID_link, "' target='_blank'>",snp$PMID,"</a>", "\npval=\n",format(snp$P, scientific = TRUE, digits=2))
# snp$pval_annot <- paste0(snp$SNP, "\nPMID: \n", snp$PMID, "\npval=\n",format(snp$P, scientific = TRUE, digits=2))
# paste0("<a href='",  snp$PMID_link, "' target='_blank'>",snp$PMID,"</a>")
snp$pval_annot <- format(snp$P, scientific = TRUE, digits=2)

# make a list of gene symbols in all datasets for checking whether gene symbol entered is valid - used for GeneSymbol later on
genes_avail <- vector()
for (i in ls()[grep("GSE", ls())]) {genes_avail <- unique(rbind(genes_avail, get(i)$SYMBOL))}

output.table <- data.frame() # initiate output table - used later in output.tableforplot()
heatmap_colors <- colorRampPalette(c("navyblue","darkgoldenrod1","firebrick4")) # heatmap colors - used in p-value plot

# server
shinyServer(function(input,output) {
  
  curr_gene <- reactive({toupper(input$curr_gene)}) #can recognize gene names even if typed lowercase
  
  GeneSymbol <- reactive({
      if (curr_gene() %in% genes_avail) {TRUE} else {FALSE}  #For generating error message when a wrong gene symbol is input.
  })
      
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
      validate(need(curr_gene() != "", "Please enter a gene id")) #Generate a error message when no gend id is input.
      
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
  validate(need(GeneSymbol() != FALSE, "Please enter a valid gene id.")) # Generate error message if the gene symbol is not right.
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
                 Lower_bound_CI = round(Lower_bound_CI, digits = 2), Upper_bound_CI = round(Upper_bound_CI, digits = 2), Comparison = "Glucocorticoid treatment vs. placebo")%>%
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
  forestplot_asthma <- function(){
      data2_Asthma = data2_Asthma()
      validate(need(nrow(data2_Asthma) != 0, "Please choose a dataset.")) #Generate the user-friendly error message
      
      text_asthma = data2_Asthma$`Study ID`

      xticks = seq(from = min(0.9, min(data2_Asthma$Lower_bound_CI)), to = max(max(data2_Asthma$Upper_bound_CI),1.2), length.out = 5)
      forestplot(as.vector(text_asthma), title = "Asthma vs. non-asthma", data2_Asthma[,c("Fold Change","Lower_bound_CI","Upper_bound_CI")], zero = 1, 
                 xlab = "Fold Change",ylab = "Studies", boxsize = 0.2, col = fpColors(lines = "darkblue", box = "royalblue", zero = "lightgrey"), lwd.ci = 2, 
                 xticks = xticks, lineheight = unit((22.5/nrow(data2_Asthma)), "cm"), graphwidth = unit(4.5, "cm"),mar = unit(c(0,0,0,0),"mm"),
                 txt_gp = fpTxtGp(xlab = gpar(cex = 1.35), ticks = gpar(cex = 1.2), title = gpar(cex = 1.2)))}

  forestplot_GC <- function(){
      data2_GC = data2_GC()
      validate(need(nrow(data2_GC) != 0, "Please choose a dataset."))
          
      text_GC = data2_GC$`Study ID`

      xticks = seq(from = min(min(0.9, data2_GC$Lower_bound_CI)), to = max(max(data2_GC$Upper_bound_CI),1.2), length.out = 5)
      forestplot(as.vector(text_GC), title = "Glucocorticoid treatment vs. placebo", data2_GC[,c("Fold Change","Lower_bound_CI","Upper_bound_CI")] ,zero = 1, 
                       xlab = "Fold Change",ylab = "Studies", boxsize = 0.15, col = fpColors(lines = "darkblue", box = "royalblue", zero = "lightgrey"), lwd.ci = 2,
                       xticks = xticks, lineheight = unit((22.5/nrow(data2_GC)), "cm"), graphwidth = unit(4.5, "cm"),mar = unit(c(0,0,0,0),"mm"),
                       txt_gp = fpTxtGp(xlab = gpar(cex = 1.35), ticks = gpar(cex = 1.2), title = gpar(cex = 1.2)))}
  
  output$forestplot_asthma = renderPlot(forestplot_asthma())
  output$forestplot_GC = renderPlot(forestplot_GC())
  
  #######################
  ## p-value levelplot ##
  #######################
  output.tableforplot2 <- reactive({output.tableforplot() %>% dplyr::rename(' '=Fold_Change, ' '=neglogofP)})
  heatmapMAT <- reactive({output.tableforplot2()})
  pval_data <- reactive({t(heatmapMAT()[10])}) #It's 5 in Maya's script. I have add some columns in the table so this number is changed.
  
  #set up max boundary for levelplot (min is fixed at 0)
  maxNLOP <- reactive({if(max(pval_data())<=1.5 & min(pval_data())>=-1.5){maxNLOP=1.5} else {maxNLOP=max(pval_data())}})
  
  # levelplot for log p-value
  pval_plot <- function() {
      levelplot(pval_data(),
                col.regions=heatmap_colors,
                xlab = NULL,
                ylab = NULL,
                main = "-log10(adjusted p-value)",
                pretty = FALSE,
                aspect = 2,
                width = 3,
                scales=list(x=list(cex=1, tck = c(0,0,0,0)),
                            y=list(cex=1, tck = c(1,0,0,0))),
                at=seq(0,maxNLOP(),length.out=100))}
  
  output$pval_plot_outp <- renderPlot({pval_plot()})
  
  ###############################
  ## Gene, SNP and TFBS tracks ##
  ###############################
    gene_tracks <- function() {
      tfbs_subs <- unique(filter(tfbs, GENE3==curr_gene()))
      gene_subs <- unique(filter(gene_locations, GENE3==curr_gene()))
      snp_subs <- unique(filter(snp, GENE3==curr_gene()))
      
      gen <- "hg19"
      
      #gene - this track shows up for all genes
      gr_gene <- GRanges(seqnames = gene_subs$CHR, ranges = IRanges(start = gene_subs$START, end = gene_subs$STOP))
      chr <- as.character(unique(seqnames(gr_gene)))
      atrack_gene <- Gviz::GeneRegionTrack(gr_gene, name="Exons", stacking="dense", fill = "dodgerblue3")
      # atrack_gene <- Gviz::GeneRegionTrack(geneModels, genome = gen, chromosome = chr, name = "Gene Model")
      itrack <- IdeogramTrack(genome = gen, chromosome = chr) # this step really slows the app down...but only the first time?
      gtrack <- GenomeAxisTrack()
      
      #tfbs - if statement b/c many genes don't have one
      if (nrow(tfbs_subs) > 0) {
          gr_tfbs <- GRanges(seqnames = tfbs_subs$CHR, ranges = IRanges(start = tfbs_subs$START, end = tfbs_subs$STOP))
          atrack_tfbs <- Gviz::AnnotationTrack(gr_tfbs, name="NR3C1 binding sites", stacking="dense", fill = tfbs_subs$color)
          feature(atrack_tfbs) <- ""
          }
      
      #snp - if statement b/c many genes don't have one
      if (nrow(snp_subs) > 0) {
          gr_snp <- GRanges(seqnames = snp_subs$CHR, ranges = IRanges(start = snp_subs$START, end = snp_subs$STOP))
          atrack_snp <- Gviz::AnnotationTrack(gr_snp, name="SNPs", stacking="dense")
          feature(atrack_snp) <- snp_subs$pval_annot
      }

      #output depends on whether there is are TFBS & SNPs for a given gene    
      if ((nrow(tfbs_subs) > 0) & (nrow(snp_subs) > 0)) {
          plotTracks(list(gtrack, atrack_gene, atrack_tfbs, atrack_snp, itrack), featureAnnotation = "feature", fontcolor.feature = "darkblue", cex.feature=0.85, sizes=c(1,1.25,1.25,1.25,0.5), col=NULL,  background.panel = "#FFFEDB", background.title = "darkblue")
      } else if (nrow(tfbs_subs) > 0) {
          plotTracks(list(gtrack, atrack_gene, atrack_tfbs, itrack), sizes=c(1,1.25,1.25,0.5), col=NULL,  background.panel = "#FFFEDB", background.title = "darkblue")
      } else if (nrow(snp_subs) > 0) {
          plotTracks(list(gtrack, atrack_gene, atrack_snp, itrack), featureAnnotation = "feature", fontcolor.feature = "darkblue", cex.feature=0.85, sizes=c(1,1.25,1.25,0.5), col=NULL, background.panel = "#FFFEDB", background.title = "darkblue")
      } else {
          plotTracks(list(gtrack, atrack_gene, itrack), sizes=c(1,1.25,0.5), col=NULL)
      }
  }

    output$gene_tracks_outp <- renderPlot({gene_tracks()})

  ######################
  ## Download buttons ##
  ######################
  graphgene=reactive({curr_gene()})
    
  output$asthma_fc_download <- downloadHandler(
    filename= function(){paste0("fold_change_asthma_", graphgene(), "_", Sys.Date(), ".png")},
    content=function(file){
      png(file)
      forestplot_asthma()
      dev.off()})
  
  output$GC_fc_download <- downloadHandler(
      filename= function(){paste0("fold_change_GC_", graphgene(), "_", Sys.Date(), ".png")},
      content=function(file){
          png(file)
          forestplot_GC()
          dev.off()})
  
  output$pval_download <- downloadHandler(
      filename= function(){paste0("-log(pval)_heatmap_", graphgene(), "_", Sys.Date(), ".png")},
      content=function(file){
          png(file)
          pval_plot()
          dev.off()})
  
  output$gene_tracks_download <- downloadHandler(
      filename= function(){paste0("gene_tracks_", graphgene(), "_", Sys.Date(), ".png")},
      content=function(file){
          png(file)
          gene_tracks()
          dev.off()})
  
  output$table_download <- downloadHandler(filename = function() {paste0('Asthma&GC_data_summary_table_',graphgene(), Sys.Date(), '.csv')},
                                                  content = function(file) {write.csv(rbind(tableforgraph_Asthma(), tableforgraph_GC()), file, row.names=FALSE)})
})