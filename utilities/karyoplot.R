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

#EQTL data
eqtl_lung <- toGRanges(readRDS("realgar_data/eQTL/Lung.RDS"))
eqtl_ms <- toGRanges(readRDS("realgar_data/eQTL/Muscle.RDS"))
eqtl_wb <- toGRanges(readRDS("realgar_data/eQTL/Blood.RDS"))

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
  
  kpAddLabels(kp, labels = label_name, r0 = out.at$r0, r1=out.at$r1, cex=1, srt=90, pos=1, label.margin = 0.14) # margin 0.3 -> label.margin 0.2
  
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
make_karyoplot <- function(gene.region, ggwas_df, gabriel_gwas_df, fer_gwas_df,
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
  delta_gwas=0.04 # proportion of GWAS track. 0.06 used before
  delta_locus=0.18 # proportion of individual locuszoom
  delta_bw=0.2 # proportion of all plots of Encode histone mark bigwig files
  delta_hmm=delta_gre=0.01
  delta_gene=0.03 # proportion of gene track. 0.05 used before
  cex_gwas=1.0
  cex_gre=1.0
  cex_chr=cex_gene=1.0
  
  #gene region
  #if (region$strand == -1){
  #  gene.region <- toGRanges(data.frame(chr=region$Chromosome, start=region$Start, end=region$End+20000,genome="hg38"))
  #} else if (region$strand == 1){
  #  gene.region <- toGRanges(data.frame(chr=region$Chromosome, start=region$Start-20000, end=region$End,genome="hg38"))
  #}
  #gene.region <- toGRanges(data.frame(chr=region$Chromosome, start=region$Start-20000, end=region$End+20000,genome="hg19"))
  
  #Add genes data
  print("plot gene track")
  genes.data <- makeGenesDataFromTxDb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                                      karyoplot=plotKaryotype(zoom = gene.region),
                                      plot.transcripts = TRUE,
                                      plot.transcripts.structure = TRUE)
  
  ## Merge transcripts and get gene symbols
  genes.data <- addGeneNames(genes.data)
  genes.data <- mergeTranscripts(genes.data)
  
  #KP initialize
  kp <- plotKaryotype(zoom = gene.region, cex=cex_chr, plot.params = pp)
  #kpAddMainTitle(kp, paste0("Asthma and Glucocorticoid Associations for ",region$symbol),cex=2,label.margin = 0.035)
  
  #Add ChromHMM track
  #Get plot
  #units
  unit <- 10^floor(log10((end(gene.region)+20000)-(start(gene.region)-20000)))
  minor_unit <- unit/2
  kpAddBaseNumbers(kp, tick.dist = unit, minor.tick.dist = minor_unit,add.units = F, cex=1.0, tick.len = 1)
  
  #genes
  kpPlotGenes(kp, data=genes.data, r0=0, r1=delta_gene, gene.name.cex = cex_gene)
  r1=0.1
  
  #r0=0.15
  #r1=r0+0.03
  #kpPlotRegions(kp, K562.hmm, col=K562.hmm$V6, r0=r0, r1=r1,avoid.overlapping=F)
  #kpAddLabels(kp, labels = "Chromatin\nState (HMM)", r0=r0, r1=r1, cex=1.0)
  
  #Add eQTL data
  r0=r1
  r1=r0+0.02
  print("plot lung eQTL track")
  no_eqtl_lung <- setdiff(gene.region, eqtl_lung, ignore.strand=TRUE)
  kpPlotRegions(kp, no_eqtl_lung, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, eqtl_lung, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "Lung", r0=r0, r1=r1, cex=1.0, label.margin = 0.035)
  
  r0=r1
  r1=r0+0.02
  g0=r1+0.01 #for gwas title
  no_eqtl_ms <- IRanges::setdiff(gene.region, eqtl_ms, ignore.strand=TRUE)
  print("plot skeletal muscle eQTL track")
  kpPlotRegions(kp, no_eqtl_ms, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, eqtl_ms, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "Muscle Skeletal", r0=r0, r1=r1, cex=1.0, label.margin = 0.035)
  
  r0=r1
  r1=r0+0.02
  g0=r1+0.01 #for gwas title
  print("plot whole blood eQTL track")
  no_eqtl_wb <- IRanges::setdiff(gene.region, eqtl_wb, ignore.strand=TRUE)
  kpPlotRegions(kp, no_eqtl_wb, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, eqtl_wb, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "Whole Blood", r0=r0, r1=r1, cex=1.0, label.margin = 0.035)
  kpAddLabels(kp, labels = "eQTL", r0=0.1, r1=r1, cex=1, srt=90, pos=1, label.margin = 0.14)
  
  #Get GWAS SNPS
  ## GWAS plot function
  gwas_plot_func <- function(df, label_gwas) {
    #Convert data
    df <- df[,c("chromosome", "start", "end", "SNP", "color")]
    df$color <- ifelse(is.na(df$color),"#000000",df$color)
    df_gr <- toGRanges(df)
    
    #Lines
    kpPlotRegions(kp, data = df_gr, col=df_gr$color, r0=r0, r1=r1-0.03, cex=4.5,srt=90)
    kpPlotMarkers(kp, data = df_gr, labels = as.character(df_gr$SNP), r0=r0, r1=r1, cex=0.8, label.margin = 0.035,text.orientation = "horizontal",adjust.label.position=T, line.color = df_gr$color)
    kpAddLabels(kp, labels = label_gwas, label.margin = 0.035, r0=r0, r1=r1, cex=cex_gwas)
  }
  
  
  ## Add GWAS data
  #EVE ALL
  print("plot EVE_ALL track")
  if(!is.null(eve_gwas_all_df) & length(eve_gwas_all_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=eve_gwas_all_df, label_gwas="EVE_ALL")    
  }
    
  #EVE EA
  print("plot EVE_EA track")
  if(!is.null(eve_gwas_ea_df) & length(eve_gwas_ea_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=eve_gwas_ea_df, label_gwas="EVE_EA") 
  }
  
  #EVE AA
  print("plot EVE_AA track")
  if(!is.null(eve_gwas_aa_df) & length(eve_gwas_aa_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=eve_gwas_aa_df, label_gwas="EVE_AA") 
  }

  #EVE LA
  print("plot EVE_LA track")
  if(!is.null(eve_gwas_la_df) & length(eve_gwas_la_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=eve_gwas_la_df, label_gwas="EVE_LA") 
  }
  
    
  #FER
  print("plot ferreira track")
  if(!is.null(fer_gwas_df) & length(fer_gwas_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=fer_gwas_df, label_gwas="FERREIRA") 
  }
  

  #GABRIEL
  print("plot gabriel track")
  if(!is.null(gabriel_gwas_df) & length(gabriel_gwas_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=gabriel_gwas_df, label_gwas="GABRIEL") 
  }
  
 
  #GRASP
  print("plot GRASP track")
  if(!is.null(ggwas_df) & length(ggwas_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=ggwas_df, label_gwas="GRASP") 
  }
  

  #TAGC multi-ethnic groups
  print("plot TAGC_Multi track")
  if(!is.null(tagc_gwas_multi_df) & length(tagc_gwas_multi_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=tagc_gwas_multi_df, label_gwas="TAGC_Multi") 
  }

  #TAGC European population
  print("plot TAGC_EUR track")
  if(!is.null(tagc_gwas_eur_df) & length(tagc_gwas_eur_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=tagc_gwas_eur_df, label_gwas="TAGC_EUR") 
  }
  
  #UKBB asthma
  print("plot UKBB_asthma track")
  if(!is.null(ukbb_gwas_asthma_df) & length(ukbb_gwas_asthma_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=ukbb_gwas_asthma_df, label_gwas="UKBB_asthma") 
  }

  #UKBB COPD
  print("plot UKBB_COPD track")
  if(!is.null(ukbb_gwas_copd_df) & length(ukbb_gwas_copd_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=ukbb_gwas_copd_df, label_gwas="UKBB_COPD") 
  }

  #UKBB ACO
  print("plot UKBB_ACO track")
  if(!is.null(ukbb_gwas_aco_df) & length(ukbb_gwas_aco_df) > 0){
    r0=r1+0.01
    r1=r0+delta_gwas
    gwas_plot_func(df=ukbb_gwas_aco_df, label_gwas="UKBB_ACO")
  }
  
    
  #Get GR gre sites
  #Add gre Data
  print("plot GRE track")
  r0 = r1 + 0.01
  r1= r0 + 0.01
  no_gre <- IRanges::setdiff(gene.region, gre, ignore.strand=FALSE)
  kpPlotRegions(kp, no_gre, col="#efefef", r0=r0, r1=r1, avoid.overlapping=F)
  kpPlotRegions(kp, gre, col="#710C04", r0=r0, r1=r1, avoid.overlapping=F)
  kpAddLabels(kp, labels = "GRE", r0=r0, r1=r1, cex=cex_gre, label.margin = 0.035)
  
  #ChIP-Seq
  r0 = r1 + 0.01
  print("plot ASM GR track")
  ASM_end <- track_plot_func(r0=r0, histone.marks = rev(ASM.merge.GRs), base.url = path, start = 1, label_name = "ASM", colour = "#1B9E77", track_ymax="visible.region", n_group=2,kp)
  print("plot BE GR track")
  BE_end <- track_plot_func(r0=r0+0.01, histone.marks = rev(BE.merge.GRs), base.url = path, start = ASM_end+1, label_name = "BEAS-2B", colour = "#D95F02", track_ymax="visible.region", n_group=2,kp)
  print("plot A549 GR track")
  A549_end <- track_plot_func(r0=r0+0.01, histone.marks = rev(A549.merge.GRs), base.url = path, start = BE_end+1, label_name = "A549", colour = "#7570B3", track_ymax="visible.region", n_group=2,kp)

}

