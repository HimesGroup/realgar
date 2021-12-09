# Load libaries in global.R
#library(karyoploteR, quietly = T) #Need R > 3.5
#library(TxDb.Hsapiens.UCSC.hg38.knownGene, quietly = T)

# Data Loading

#ChromHMM
# K562.hmm <- toGRanges("realgar_data/chrhmm/wgEncodeBroadHmmNhlfHMM.bed.gz")
# 
# chromHMM <- as.vector(unique(K562.hmm$V4))
# chromHMM_sort <- str_sort(chromHMM, numeric = TRUE)
# K562.hmm$V4 <- factor(K562.hmm$V4, levels=chromHMM_sort)
# #colors = c("Red","Light Coral","Purple","Orange","Orange","Yellow","Yellow","Blue","Dark Green","Dark Green","Light Green","Slate Gray","Light Gray","Light Gray","Light Gray")
# colors = c("#FF0000","#CD5C5C","#800080","#FFA500","#FFA500","#FFFF00","#FFFF00","#00BFFF","#006400","#006400","#90EE90","#808080","#D3D3D3","#D3D3D3","#D3D3D3")
# K562.hmm$V6 <- as.character(factor(K562.hmm$V4, labels = colors))
# K562.hmm <- toGRanges(K562.hmm)

#Genes
# genes <- read.csv("/home/avantika/all_genes_hg19.csv")
# genes$Chromosome <- paste0("chr",genes$Chromosome)
# saveRDS(genes,"/home/avantika/gene_symbol_coords_hg19.RDS")

#GRE data
gre <- toGRanges(readRDS("realgar_data/gre.RDS")) 
#gre$V6 <- inferno(50)[as.numeric(cut(gre$V5,breaks = 50))]

#EQTL data
eqtl_lung <- toGRanges(readRDS("realgar_data/eQTL/Lung.RDS") %>% dplyr::select(-id))
eqtl_ms <- toGRanges(readRDS("realgar_data/eQTL/Muscle_Skeletal.RDS") %>% dplyr::select(-id))
eqtl_wb <- toGRanges(readRDS("realgar_data/eQTL/Whole_Blood.RDS") %>% dplyr::select(-id))

#Bigwig ChIP-Seq datasets
ASM.merge.GRs <- c( GR_dex="ASM.GR.dex.bw", GR_EtOH="ASM.GR.control.bw", RNAP2_dex="ASM.RNAP2.dex.bw", RNAP2_EtOH="ASM.RNAP2.control.bw")
BE.merge.GRs <- c( GR_dex="BE.GR.dex.bw", GR_EtOH="BE.GR.control.bw", RNAP2_dex="BE.RNAP2.dex.bw", RNAP2_EtOH="BE.RNAP2.control.bw")
A549.merge.GRs <- c( GR_dex="A549.GR.dex.bw", GR_EtOH="A549.GR.control.bw", RNAP2_dex="A549.RNAP2.dex.bw", RNAP2_EtOH="A549.RNAP2.control.bw")
total.tracks <- length(ASM.merge.GRs) + length(BE.merge.GRs) + length(A549.merge.GRs)
path="https://public.himeslab.dev/bigwig_merge/"

#helper functions

#ChIP-Seq ymax
bigwig_ymax_func <- function(file) {
  data <- rtracklayer::BigWigFile(file)
  wig.data <- rtracklayer::import(data, format = "bigWig", selection=gr)
  ymax <- max(0, max(mcols(wig.data)[,1]))
  return(ymax)
}

# add UCSC track
track_plot_func <- function(r0, histone.marks, base.url, start, label_name, colour="cadetblue2", track_ymax="visible.region", n_group=2,kp){
  end = start+length(histone.marks)-1
  out.at <- autotrack(start:end, total.tracks, margin = 0.3, r0=r0)
  bigwig.files <- paste0(base.url, histone.marks)
  
  kpAddLabels(kp, labels = label_name, r0 = out.at$r0, r1=out.at$r1, cex=1.5, srt=90, pos=1, label.margin = 0.14) # margin 0.3 -> label.margin 0.2
  
  ymax_vect <- rep(0, length(histone.marks)) # parameter vector for ymax
  if (length(track_ymax)==1) {
  ymax_vect[seq_len(length(histone.marks))] <- track_ymax
  } else {
    ymax_vect <- track_ymax
  }
  
  
  #RNAP2 tracks
  at1 <- autotrack(1, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
  kp1 <- kpPlotBigWig(kp, data=bigwig.files[1], ymax=ymax_vect[1],
                     r0=at1$r0, r1=at1$r1, col = colour)
  at2 <- autotrack(2, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
  kp2 <- kpPlotBigWig(kp1, data=bigwig.files[2], ymax=ymax_vect[2],
                      r0=at2$r0, r1=at2$r1, col = colour)
  computed.ymax <- max(ceiling(kp1$latest.plot$computed.values$ymax),ceiling(kp2$latest.plot$computed.values$ymax))
  kpAxis(kp1, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax, r0=at1$r0, r1=at1$r1, cex=1.1)
  kpAddLabels(kp1, labels = names(histone.marks)[1], r0=at1$r0, r1=at1$r1, cex=1.0, label.margin = 0.035)
  kpAxis(kp2, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax, r0=at2$r0, r1=at2$r1, cex=1.1)
  kpAddLabels(kp2, labels = names(histone.marks)[2], r0=at2$r0, r1=at2$r1, cex=1.0, label.margin = 0.035)
  
  #GR tracks
  at3 <- autotrack(3, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
  kp3 <- kpPlotBigWig(kp2, data=bigwig.files[3], ymax=ymax_vect[3],
                      r0=at3$r0, r1=at3$r1, col = colour)
  at4 <- autotrack(4, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
  kp4 <- kpPlotBigWig(kp3, data=bigwig.files[4], ymax=ymax_vect[4],
                      r0=at4$r0, r1=at4$r1, col = colour)
  computed.ymax <- max(ceiling(kp3$latest.plot$computed.values$ymax),ceiling(kp4$latest.plot$computed.values$ymax))
  kpAxis(kp3, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax, r0=at3$r0, r1=at3$r1, cex=1.1)
  kpAddLabels(kp3, labels = names(histone.marks)[3], r0=at3$r0, r1=at3$r1, cex=1.0, label.margin = 0.035)
  kpAxis(kp4, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax, r0=at4$r0, r1=at4$r1, cex=1.1)
  kpAddLabels(kp4, labels = names(histone.marks)[4], r0=at4$r0, r1=at4$r1, cex=1.0, label.margin = 0.035)
  
  
  return(end)
}

## KaryoplotR function
make_karyoplot <- function(region, ggwas_df, gabriel_gwas_df, fer_gwas_df,
                           eve_gwas_all_df, eve_gwas_ea_df, eve_gwas_aa_df, eve_gwas_la_df, 
                           tagc_gwas_multi_df, tagc_gwas_eur_df,
                           ukbb_gwas_asthma_df, ukbb_gwas_copd_df, ukbb_gwas_aco_df){
  
  #Plot margins
  ## Plot karyoplotR
  # generate plot parameters
  pp <- getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 1.34/5.5 # 4.5 in -> 0.3 margin
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0
  
  
  # need to adjust the proportion of delta_locus and delta_bw based on the number of locuszoom plots
  delta_r=0.02
  delta_locus=0.18 # proportion of individual locuszoom
  delta_bw=0.2 # proportion of all plots of Encode histone mark bigwig files
  delta_hmm=delta_gre=0.01
  
  #gene region
  if (region$strand == -1){
    gene.region <- toGRanges(data.frame(chr=region$Chromosome, start=region$Start, end=region$End+20000,genome="hg38"))
  } else if (region$strand == 1){
    gene.region <- toGRanges(data.frame(chr=region$Chromosome, start=region$Start-20000, end=region$End,genome="hg38"))
  }
  #gene.region <- toGRanges(data.frame(chr=region$Chromosome, start=region$Start-20000, end=region$End+20000,genome="hg19"))
  
  #Add genes data
  genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                      karyoplot=plotKaryotype(zoom = gene.region),
                                      plot.transcripts = TRUE,
                                      plot.transcripts.structure = TRUE)
  
  ## Merge transcripts and get gene symbols
  genes.data <- addGeneNames(genes.data)
  genes.data <- mergeTranscripts(genes.data)
  
  #KP initialize
  kp <- plotKaryotype(zoom = gene.region, cex=1.5, plot.params = pp)
  #kpAddMainTitle(kp, paste0("Asthma and Glucocorticoid Associations for ",region$symbol),cex=2,label.margin = 0.035)
  
  #Add ChromHMM track
  #Get plot
  #units
  unit <- 10^floor(log10((region$End+20000)-(region$Start-20000)))
  minor_unit <- unit/2
  kpAddBaseNumbers(kp, tick.dist = unit, minor.tick.dist = minor_unit,add.units = F, cex=1.0, tick.len = 1)
  
  #genes
  kpPlotGenes(kp, data=genes.data, r0=0, r1=0.05, gene.name.cex = 1.2)
  r1=0.1
  
  #r0=0.15
  #r1=r0+0.03
  #kpPlotRegions(kp, K562.hmm, col=K562.hmm$V6, r0=r0, r1=r1,avoid.overlapping=F)
  #kpAddLabels(kp, labels = "Chromatin\nState (HMM)", r0=r0, r1=r1, cex=1.0)
  
  #Add eQTL data
  # r0=r1
  # r1=r0+0.02
  # if (length(eqtl_lung)>0){
  # no_eqtl_lung <- IRanges::setdiff(gene.region, eqtl_lung, ignore.strand=TRUE)
  # kpPlotRegions(kp, eqtl_lung, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  # kpPlotRegions(kp, no_eqtl_lung, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  # } else {kpPlotRegions(kp, gene.region, col="#efefef", r0=r0, r1=r1)}
  # kpAddLabels(kp, labels = "Lung", r0=r0, r1=r1, cex=1.0, label.margin = 0.035)
  
  r0=r1
  r1=r0+0.02
  no_eqtl_lung <- setdiff(gene.region, eqtl_lung, ignore.strand=TRUE)
  kpPlotRegions(kp, no_eqtl_lung, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, eqtl_lung, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "Lung", r0=r0, r1=r1, cex=1.0, label.margin = 0.035)
  
  r0=r1
  r1=r0+0.02
  g0=r1+0.01 #for gwas title
  no_eqtl_ms <- IRanges::setdiff(gene.region, eqtl_ms, ignore.strand=TRUE)
  # gre_a549$itemRgb <- gre_a549$V6
  # no_gre_a549$itemRgb <- "#efefef"
  kpPlotRegions(kp, no_eqtl_ms, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, eqtl_ms, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "Muscle Skeletal", r0=r0, r1=r1, cex=1.0, label.margin = 0.035)
  
  r0=r1
  r1=r0+0.02
  g0=r1+0.01 #for gwas title
  no_eqtl_wb <- IRanges::setdiff(gene.region, eqtl_wb, ignore.strand=TRUE)
  kpPlotRegions(kp, no_eqtl_wb, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, eqtl_wb, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "Whole Blood", r0=r0, r1=r1, cex=1.0, label.margin = 0.035)
  kpAddLabels(kp, labels = "eQTL", r0=0.1, r1=r1, cex=1.5, srt=90, pos=1, label.margin = 0.14)
  
  #Get GWAS SNPS
  ## Add GWAS data
  #EVE ALL
  if(!is.null(eve_gwas_all_df) & length(eve_gwas_all_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    eve_gwas_all_df$color <- ifelse(is.na(eve_gwas_all_df$color),"#000000",eve_gwas_all_df$color)
    eve_gwas <- toGRanges(eve_gwas_all_df)
    
    #Lines
    kpPlotRegions(kp, data = eve_gwas,col=eve_gwas$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = eve_gwas, labels = as.character(eve_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color = eve_gwas$color)
    kpAddLabels(kp, labels = "EVE_ALL", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}
    
  #EVE EA
  if(!is.null(eve_gwas_ea_df) & length(eve_gwas_ea_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    eve_gwas_ea_df$color <- ifelse(is.na(eve_gwas_ea_df$color),"#000000",eve_gwas_ea_df$color)
    eve_gwas <- toGRanges(eve_gwas_ea_df)
    
    #Lines
    kpPlotRegions(kp, data = eve_gwas,col=eve_gwas$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = eve_gwas, labels = as.character(eve_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color = eve_gwas$color)
    kpAddLabels(kp, labels = "EVE_EA", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}
  
  #EVE AA
  if(!is.null(eve_gwas_aa_df) & length(eve_gwas_aa_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    eve_gwas_aa_df$color <- ifelse(is.na(eve_gwas_aa_df$color),"#000000",eve_gwas_aa_df$color)
    eve_gwas <- toGRanges(eve_gwas_aa_df)
    
    #Lines
    kpPlotRegions(kp, data = eve_gwas,col=eve_gwas$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = eve_gwas, labels = as.character(eve_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color = eve_gwas$color)
    kpAddLabels(kp, labels = "EVE_AA", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}

  #EVE LA
  if(!is.null(eve_gwas_la_df) & length(eve_gwas_la_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    eve_gwas_la_df$color <- ifelse(is.na(eve_gwas_la_df$color),"#000000",eve_gwas_la_df$color)
    eve_gwas <- toGRanges(eve_gwas_la_df)
    
    #Lines
    kpPlotRegions(kp, data = eve_gwas,col=eve_gwas$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = eve_gwas, labels = as.character(eve_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color = eve_gwas$color)
    kpAddLabels(kp, labels = "EVE_LA", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}
  
    
  #FER
  if(!is.null(fer_gwas_df) & length(fer_gwas_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    fer_gwas_df$color <- ifelse(is.na(fer_gwas_df$color),"#000000",fer_gwas_df$color)
    fer_gwas <- toGRanges(fer_gwas_df)
    
    #Lines
    kpPlotRegions(kp, data = fer_gwas,col=fer_gwas$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = fer_gwas, labels = as.character(fer_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color=fer_gwas$color)
    kpAddLabels(kp, labels = "FERREIRA", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}
  

  #GABRIEL
  if(!is.null(gabriel_gwas_df) & length(gabriel_gwas_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    gabriel_gwas_df$color <- ifelse(is.na(gabriel_gwas_df$color),"#000000",gabriel_gwas_df$color)
    gabriel_gwas <- toGRanges(gabriel_gwas_df)
    
    #Lines
    kpPlotRegions(kp, data = gabriel_gwas,col=gabriel_gwas$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = gabriel_gwas, labels = as.character(gabriel_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T, line.color=gabriel_gwas$color)
    kpAddLabels(kp, labels = "GABRIEL", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}
  
 
  #GRASP
  if(!is.null(ggwas_df) & length(ggwas_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    ggwas_df$color <- ifelse(is.na(ggwas_df$color),"#000000",ggwas_df$color)
    ggwas <- toGRanges(ggwas_df)
    
    #Lines
    kpPlotRegions(kp, data = ggwas,col=ggwas$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = ggwas, labels = as.character(ggwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color=ggwas$color)
    kpAddLabels(kp, labels = "GRASP", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}
  

  #TAGC multi-ethnic groups
  if(!is.null(tagc_gwas_multi_df) & length(tagc_gwas_multi_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    tagc_gwas_multi_df$color <- ifelse(is.na(tagc_gwas_multi_df$color),"#000000",tagc_gwas_multi_df$color)
    tagc_gwas <- toGRanges(tagc_gwas_multi_df)
    
    #Lines
    kpPlotRegions(kp, data = tagc_gwas,col=tagc_gwas$color, r0=r0, r1=r1-0.03, cex=5.0,srt=90)
    kpPlotMarkers(kp, data = tagc_gwas, labels = as.character(tagc_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color=tagc_gwas$color)
    kpAddLabels(kp, labels = "TAGC_Mutli", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}

  #TAGC European population
  if(!is.null(tagc_gwas_eur_df) & length(tagc_gwas_eur_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    tagc_gwas_eur_df$color <- ifelse(is.na(tagc_gwas_eur_df$color),"#000000",tagc_gwas_eur_df$color)
    tagc_gwas <- toGRanges(tagc_gwas_eur_df)
    
    #Lines
    kpPlotRegions(kp, data = tagc_gwas,col=tagc_gwas$color, r0=r0, r1=r1-0.03, cex=5.0,srt=90)
    kpPlotMarkers(kp, data = tagc_gwas, labels = as.character(tagc_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color=tagc_gwas$color)
    kpAddLabels(kp, labels = "TAGC_EUR", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)}
  
    
    #points
    # ymax = ceiling(max(tagc_gwas$neg_log_p))
    # if(length(tagc_gwas)==1){
    #   ymin = 0
    # } else {ymin = floor(min(tagc_gwas$neg_log_p))}
    # 
    # kpPoints(kp, chr=tagc_gwas$chromosome, x=tagc_gwas$start, y=tagc_gwas$neg_log_p, data.panel=1, ymin=ymin, ymax=ymax, col="black",bg=tagc_gwas$color,r0=r0, r1=r1, pch=21, cex=1.0)
    # kpText(kp, chr=tagc_gwas$chromosome, x=tagc_gwas$start, y=tagc_gwas$neg_log_p+(ymax-ymin)/5, ymin=ymin, ymax=ymax, r0=r0, r1=r1, labels=tagc_gwas$snp, col="blue", cex=1.0)
    # kpAddLabels(kp, labels = "TAGC\nGWAS", label.margin = 0.035, r0=r0, r1=r1, cex=1.0)}

    #Get ChIP-Seq tracks
    #ChIP-Seq tracks
    # total.tracks <- length(histone.marks)
    # out.at <- autotrack(1:length(histone.marks), total.tracks, margin = 0.3, r0=0.60)
    # kpAddLabels(kp, labels = "ASM", r0 = out.at$r0, r1=out.at$r1, cex=2.0,
    #             srt=90, pos=1, label.margin = 0.14)
    # ymax <- c(150,150,60,60)
    # for(i in seq_len(length(histone.marks))) {
    #   bigwig.file <- paste0(base.url, histone.marks[i])
    #   at <- autotrack(i, length(histone.marks), r0=out.at$r0, r1=out.at$r1, margin = 0.1)
    #   kp<- kpPlotBigWig(kp, data=bigwig.file, ymax=ymax[i],r0=at$r0, r1=at$r1, col = colors[i])
    #   computed.ymax <- kp$latest.plot$computed.values$ymax
    #   kpAxis(kp, ymin=0, ymax=computed.ymax, tick.pos = computed.ymax, 
    #          r0=at$r0, r1=at$r1, cex=1.1)
    #   kpAddLabels(kp, labels = names(histone.marks)[i], r0=at$r0, r1=at$r1, 
    #               cex=1.0, label.margin = 0.035)
    # }
  
  #Add GWAS title track
  # if (exists("g1")){
  #   kpAddLabels(kp, labels = "GWAS", r0=0.1, r1=r1, cex=1.5, srt=90, pos=1, label.margin = 0.14)
  # }
  
  #UKBB asthma
  if(!is.null(ukbb_gwas_asthma_df) & length(ukbb_gwas_asthma_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    ukbb_gwas_asthma_df$color <- ifelse(is.na(ukbb_gwas_asthma_df$color),"#000000",ukbb_gwas_asthma_df$color)
    ukbb_gwas <- toGRanges(ukbb_gwas_asthma_df)
    
    #Lines
    kpPlotRegions(kp, data = ukbb_gwas,col=ukbb_gwas$color, r0=r0, r1=r1-0.03, cex=5.0,srt=90)
    kpPlotMarkers(kp, data = ukbb_gwas, labels = as.character(ukbb_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color=ukbb_gwas$color)
    kpAddLabels(kp, labels = "UKBB_asthma", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)
    }

  #UKBB COPD
  if(!is.null(ukbb_gwas_copd_df) & length(ukbb_gwas_copd_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    ukbb_gwas_copd_df$color <- ifelse(is.na(ukbb_gwas_copd_df$color),"#000000",ukbb_gwas_copd_df$color)
    ukbb_gwas <- toGRanges(ukbb_gwas_copd_df)
    
    #Lines
    kpPlotRegions(kp, data = ukbb_gwas,col=ukbb_gwas$color, r0=r0, r1=r1-0.03, cex=5.0,srt=90)
    kpPlotMarkers(kp, data = ukbb_gwas, labels = as.character(ukbb_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color=ukbb_gwas$color)
    kpAddLabels(kp, labels = "UKBB_COPD", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)
  }

  #UKBB ACO
  if(!is.null(ukbb_gwas_aco_df) & length(ukbb_gwas_aco_df) > 0){
    r0=r1+0.01
    r1=r0+0.06
    
    #Convert data
    ukbb_gwas_aco_df$color <- ifelse(is.na(ukbb_gwas_aco_df$color),"#000000",ukbb_gwas_aco_df$color)
    ukbb_gwas <- toGRanges(ukbb_gwas_aco_df)
    
    #Lines
    kpPlotRegions(kp, data = ukbb_gwas,col=ukbb_gwas$color, r0=r0, r1=r1-0.03, cex=5.0,srt=90)
    kpPlotMarkers(kp, data = ukbb_gwas, labels = as.character(ukbb_gwas$snp), r0=r0, r1=r1, cex=0.8,
                  label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T,line.color=ukbb_gwas$color)
    kpAddLabels(kp, labels = "UKBB_ACO", label.margin = 0.035, r0=r0, r1=r1, cex=1.2)
  }
  
    
  #Get GR gre sites
  #Add gre Data
  r0 = r1 + 0.01
  r1= r0 + 0.03
  no_gre <- IRanges::setdiff(gene.region, gre, ignore.strand=FALSE)
  kpPlotRegions(kp, no_gre, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, gre, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "GRE", r0=r0, r1=r1, cex=1.2, label.margin = 0.035)
  
  #ChIP-Seq
  r0 = r1 + 0.01
  ASM_end <- track_plot_func(r0=r0, histone.marks = rev(ASM.merge.GRs), base.url = path, start = 1, label_name = "ASM", colour = "#1B9E77", track_ymax="visible.region", n_group=2,kp)
  BE_end <- track_plot_func(r0=r0+0.01, histone.marks = rev(BE.merge.GRs), base.url = path, start = ASM_end+1, label_name = "BEAS-2B", colour = "#D95F02", track_ymax="visible.region", n_group=2,kp)
  A549_end <- track_plot_func(r0=r0+0.01, histone.marks = rev(A549.merge.GRs), base.url = path, start = BE_end+1, label_name = "A549", colour = "#7570B3", track_ymax="visible.region", n_group=2,kp)

}

