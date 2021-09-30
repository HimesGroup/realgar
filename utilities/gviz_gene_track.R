gene_tracks <- function() {
  validate(need(curr_gene() != "", "Please enter a gene symbol or SNP ID.")) #Generate a error message when no gene id is input.
  validate(need(GeneSymbol() != FALSE, "Please enter a valid gene symbol or SNP ID.")) # Generate error message if the gene symbol is not right.
  validate(need(nrow(UserDataset_Info()) != 0, "Please choose at least one dataset.")) #Generate a error message when no data is loaded.
  
  gene_subs <- gene_subs()
  tfbs_subs <- tfbs_subs()
  snp_subs <- snp_subs()
  snp_eve_subs <- snp_eve_subs()
  snp_gabriel_subs <- snp_gabriel_subs()
  snp_fer_subs <- snp_fer_subs()
  snp_TAGC_subs <- snp_TAGC_subs()
  
  #for better visibility, increase tfbs and snp widths -- need scaling factor b/c different genes take up different amounts of space
  smallest_start <- min(gene_subs$start, tfbs_subs$start, snp_subs$start, snp_eve_subs$start, snp_gabriel_subs$start, snp_fer_subs$start, snp_TAGC_subs$start)
  largest_end <- max(gene_subs$end, tfbs_subs$end, snp_subs$end, snp_eve_subs$end, snp_gabriel_subs$end, snp_fer_subs$end, snp_TAGC_subs$end)
  scaling_factor <- (largest_end - smallest_start)/120000
  
  tfbs_subs$end <- tfbs_subs$end + 500*scaling_factor
  snp_subs$end <- snp_subs$end + 300*scaling_factor
  snp_eve_subs$end <- snp_eve_subs$end + 300*scaling_factor
  snp_gabriel_subs$end <- snp_gabriel_subs$end + 300*scaling_factor
  snp_fer_subs$end <- snp_fer_subs$end + 300*scaling_factor
  snp_TAGC_subs$end <- snp_TAGC_subs$end + 300*scaling_factor
  
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
  if (nrow(tfbs_subs) > 0) {tfbs_track <- Gviz::AnnotationTrack(tfbs_subs, name="TF binding", fill = tfbs_subs$color, col=NULL, feature = tfbs_subs$score, grid=TRUE, col.grid="darkgrey")}
  
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
    tval_choice <- reactive({input$which_eve_pvals})  #pval_choice is responsible for dynamically coloring snps based on user selection of population
    if (!is.null(tval_choice())){
      snp_eve_track <- Gviz::AnnotationTrack(snp_eve_subs, name="SNPs (EVE)", fill = unlist(snp_eve_subs[, paste0("color_", tval_choice())]), col=NULL, feature=snp_eve_subs$snp, grid=TRUE, col.grid="darkgrey")
    } else {snp_eve_track <- Gviz::AnnotationTrack(NULL)}
    
    #rough estimate of number of stacks there will be in SNP track - for track scaling
    if (nrow(snp_eve_subs) > 1) {
      snp_eve_subs_temp <- snp_eve_subs
      snp_eve_range <- max(snp_eve_subs_temp$start) - min(snp_eve_subs_temp$start)
      snp_eve_subs_temp$start_prev <- c(0, snp_eve_subs_temp$start[1:(nrow(snp_eve_subs_temp)-1)])
      snp_eve_subs_temp$dist <- as.numeric(snp_eve_subs_temp$start) - as.numeric(snp_eve_subs_temp$start_prev)
      snp_eve_size_init <- 1.5 + as.numeric(nrow(snp_eve_subs[which(snp_eve_subs$dist < snp_eve_range/10),])) + 0.3*length(unique(gene_subs$transcript))
    } else {snp_eve_size_init <- 1.4 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_eve_subs)}
  }
  
  else {snp_eve_size_init <- 1.4 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_eve_subs)}
  
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
  
  # Ferreira SNPs track
  if (nrow(snp_fer_subs) > 0) {
    snp_fer_track <- Gviz::AnnotationTrack(snp_fer_subs, name="SNPs (Ferreira)", fill = snp_fer_subs$color, col=NULL, feature=snp_fer_subs$snp, grid=TRUE, col.grid="darkgrey")
    
    #rough estimate of number of stacks there will be in SNP track - for track scaling
    if (nrow(snp_fer_subs) > 1) {
      snp_fer_subs_temp <- snp_fer_subs
      snp_fer_range <- max(snp_fer_subs_temp$start) - min(snp_fer_subs_temp$start)
      snp_fer_subs_temp$start_prev <- c(0, snp_fer_subs_temp$start[1:(nrow(snp_fer_subs_temp)-1)])
      snp_fer_subs_temp$dist <- as.numeric(snp_fer_subs_temp$start) - as.numeric(snp_fer_subs_temp$start_prev)
      snp_fer_size_init <- 10 + as.numeric(nrow(snp_fer_subs[which(snp_fer_subs$dist < snp_fer_range),])) + 0.12*length(unique(gene_subs$transcript))
    } else {snp_fer_size_init <- 1.4 + 0.1*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_fer_subs)}
  } else {snp_fer_size_init <- 1.4 + 0.1*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_fer_subs)}
  
  # TAGC SNPs track
  if (nrow(snp_TAGC_subs) > 0) {
    pval_choice <- reactive({input$which_TAGC_pvals})  #pval_choice is responsible for dynamically coloring snps based on user selection of population
    if (!is.null(pval_choice())){
      snp_TAGC_track <- Gviz::AnnotationTrack(snp_TAGC_subs, name="SNPs (TAGC)", fill = unlist(snp_TAGC_subs[, paste0("color_", pval_choice())]), col=NULL, feature=snp_TAGC_subs$snp, grid=TRUE, col.grid="darkgrey")
    } else {snp_TAGC_track <- Gviz::AnnotationTrack(NULL)}
    
    #rough estimate of number of stacks there will be in SNP track - for track scaling
    if (nrow(snp_TAGC_subs) > 1) {
      snp_TAGC_subs_temp <- snp_TAGC_subs
      snp_TAGC_range <- max(snp_TAGC_subs_temp$start) - min(snp_TAGC_subs_temp$start)
      snp_TAGC_subs_temp$start_prev <- c(0, snp_TAGC_subs_temp$start[1:(nrow(snp_TAGC_subs_temp)-1)])
      snp_TAGC_subs_temp$dist <- as.numeric(snp_TAGC_subs_temp$start) - as.numeric(snp_TAGC_subs_temp$start_prev)
      snp_TAGC_size_init <- 1.5 + as.numeric(nrow(snp_TAGC_subs[which(snp_TAGC_subs$dist < snp_TAGC_range/10),])) + 0.3*length(unique(gene_subs$transcript))
    } else {snp_TAGC_size_init <- 1.4 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_TAGC_subs)}
  } else {snp_TAGC_size_init <- 1.4 + 0.05*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_TAGC_subs)}
  
  #track sizes - defaults throw off scaling as more tracks are added
  chrom_size <- 1.2 + 0.01*length(unique(gene_subs$transcript)) + 0.01*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.01*nrow(snp_gabriel_subs)
  axis_size <- 1 + 0.05*length(unique(gene_subs$transcript)) + 0.01*nrow(snp_subs) + 0.005*nrow(snp_eve_subs) + 0.01*nrow(snp_gabriel_subs)
  gene_size <- 2 + 0.6*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs) + 0.05*nrow(snp_eve_subs) + 0.015*nrow(snp_gabriel_subs)
  tfbs_size <- 2 + 0.075*length(unique(gene_subs$transcript)) + 0.015*nrow(snp_subs) + 0.05*nrow(snp_eve_subs) + 0.015*nrow(snp_gabriel_subs)
  snp_size <- snp_size_init #from above
  snp_eve_size <- snp_eve_size_init #from above
  snp_gabriel_size <- snp_gabriel_size_init #from above
  snp_fer_size <- snp_fer_size_init #from above
  snp_TAGC_size <- snp_TAGC_size_init #from above
  
  #select the non-empty tracks to output -- output depends on whether there are TFBS and/or SNPs for a given gene
  subset_size <- sapply(c("tfbs_subs", "snp_subs", "snp_eve_subs", "snp_gabriel_subs", "snp_fer_subs", "snp_TAGC_subs"), function(x) {nrow(get(x))}) #size of each subset
  non_zeros <- names(subset_size)[which(!(subset_size==0))] #which subsets have non-zero size
  
  df_extract <- function(x,y) { #gives name of track and track size variable for non-zero subsets (y is "track" or "size")
    if (length(non_zeros) > 0) {
      get(paste0(strsplit(x, 'subs'),y)) #trim off "subs" and append either "track" or "size"
    } else {NULL} #to avoid meltdown if no subsets were non-zero
  }
  
  #use df_extract function to get track & track size corresponding to all non-zero subsets
  #note chrom_track, axis_track and gene_track are present for all
  selected_tracks <- list(chrom_track, axis_track, gene_track, 
                          sapply(non_zeros, df_extract, y="track")$tfbs_subs, 
                          sapply(non_zeros, df_extract, y="track")$snp_subs, 
                          sapply(non_zeros, df_extract, y="track")$snp_eve_subs, 
                          sapply(non_zeros, df_extract, y="track")$snp_gabriel_subs, 
                          sapply(non_zeros, df_extract, y="track")$snp_fer_subs, 
                          sapply(non_zeros, df_extract, y="track")$snp_TAGC_subs)
  
  selected_tracks <- Filter(Negate(function(x) is.null(unlist(x))), selected_tracks) #remove null elements from list
  
  selected_sizes <- na.omit(c(chrom_size,axis_size,gene_size,
                              sapply(non_zeros, df_extract, y="size")[1], 
                              sapply(non_zeros, df_extract, y="size")[2], 
                              sapply(non_zeros, df_extract, y="size")[3], 
                              sapply(non_zeros, df_extract, y="size")[4], 
                              sapply(non_zeros, df_extract, y="size")[5],
                              sapply(non_zeros, df_extract, y="size")[6]))
  
  selected_sizes <- Filter(Negate(function(x) is.null(unlist(x))), selected_sizes) #remove null elements from list - else run into trouble in conditions when no TFBS & no SNP tracks selected
  #note: use names to extract from selected_tracks b/c it is a list vs. index to extract from selected_sizes, since this is numeric
  
  #plot tracks 
  plotTracks(selected_tracks, sizes=selected_sizes, background.panel = "#BDB6B0", background.title = "firebrick4", col.border.title = "firebrick4", groupAnnotation = "feature", fontcolor.group = "darkblue", cex.group=0.75, just.group="below", cex.title=1.1)
  
}