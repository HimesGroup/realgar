##########################
## FORESTPLOT FUNCTION ##
##########################

# Function: "forestplot_func" 
# Forest plots
forestplot_func <- function(dat, ptitle,curr_gene) {
  
  validate(need(nrow(dat) != 0, "No entries available for your selections. Please choose other options.")) #Generate a error message when no data available for selected options.
  
  # create an empty dataset for forestplot text
  tabletext <- data.frame(matrix(nrow=1, ncol=4))
  
  if ("asthma"%in%dat$App) {
    text_temp <- dat[,c("GEO_ID","Long_tissue_name","Asthma","Q Value")]
  } else {
    text_temp <- dat[,c("GEO_ID","Long_tissue_name","Treatment","Q Value")]
  }
  
  # assign column names same as the original data
  names(tabletext) <- names(text_temp) <- c("GEO_ID","Long_tissue_name","Condition","Q Value")
  tabletext[1,] <- c("GEO ID", "Tissue", "Condition", "Q-Value")
  
  #Chnage to character strings
  text_temp$GEO_ID <- as.character(text_temp$GEO_ID)
  text_temp$Long_tissue_name <- as.character(text_temp$Long_tissue_name)
  
  # assign variables to meta-analysis result row
  if (nrow(dat)>1) {
    text_temp[nrow(dat),c("GEO_ID")] <- " "
    text_temp[nrow(dat),c("Long_tissue_name")] <- " "
    text_temp[nrow(dat),3] <- "Effect size-based integration =   " # "Asthma" or "Treatment" is in column 3
  }
  
  #Change to factor
  text_temp <- text_temp %>% mutate_all(as.factor)
  
  #List of alphabets - if forestplot gives overlap/height scaling error : add more combinations to list alphabets
  mix = paste(letters,LETTERS)
  alph1 = append(letters, LETTERS)
  alphabets <- append(alph1,mix)
  
  #Text table
  # add text for individual study result
  tablevector <- as.vector(as.matrix(rbind(tabletext,text_temp)))
  levels<-alphabets[1:(nrow(text_temp)+1)]
  x<-levels
  for (j in 1:(ncol(text_temp)-1)){x<-append(x,levels)}
  #Assign position for columns
  table <- data.frame(V0 = factor(x, rev(levels)), V05 = rep(c(0,0.16,0.50,0.80),each=nrow(text_temp)+1),V1=tablevector) #c(1,1.5,2.7,3.7) #0,0.2,0.7,1.1 #0,0.2,0.6,1.0
  
  # remove double quote
  options(useFancyQuotes = FALSE)
  
  
  ##G1 plot1 - datatable
  g1 <- ggplot(table, aes(x = V05, y = V0, label = format(V1, nsmall = 1))) +
    geom_text(size = 4.2, hjust=0, vjust=0.5, #4.4
              fontface = ifelse(table$V1 %in% c("Effect size-based integration =   ","GEO ID","Tissue","Condition","Q-Value")|table$V1 == table$V1[nrow(table)] & table$V0 == table$V0[nrow(table)], 2, 1)) + 
    theme_bw() +
    geom_segment(y=nrow(dat)+1.5,yend=nrow(dat)+1.5,x=0,xend=0.90) + #x=1, xend=4.0 #x=0,xend=1.2 #xend=1.1
    geom_segment(y=nrow(dat)+0.5,yend=nrow(dat)+0.5,x=0,xend=0.90) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "lines"),
          plot.title = element_text(size=20, hjust=0.5, vjust=0.5))+
    labs(x="",y="") + #+ xlim(-0.5, 4) 
    coord_cartesian(xlim=c(0, 0.85)) #1,4.5 #3.8 #0, 1.1
  
  #Forest plot
  # table with fold changes for plot
  tableplot <- rbind(c(NA,NA,NA,NA),dat[,c("Fold Change","Lower_bound_CI","Upper_bound_CI")])	
  no_of_values <- nrow(dat)+1
  tableplot <- data.frame(sapply(tableplot, as.numeric))
  tableplot$group <-table$V0[1:no_of_values]
  
  #X-ticks
  xticks = seq(from = min(0.9, min(dat$Lower_bound_CI)), to = max(max(dat$Upper_bound_CI),1.2), length.out = 5)
  xticks = round(xticks,2)
  
  #Colors: (inferno(8002) can be used)
  breaks <- c(seq(0,8,by=0.001), Inf) # this sets max universally at 8 (else highest one OF THE SUBSET would be the max)
  reds = rev(c(rep("#F58178",1001),rep("#F26256",1001),rep("#EF3B2C",1000),rep("#BF2F23",1000),rep("#99261C",1000),rep("#7A1E16",1000),rep("#621812",1000), rep("#4E130E",1000)))
  b_clrs  <- reds[findInterval(dat$neglogofP, breaks)] #8002 is length(breaks) - ensures there are enough colors
  
  #Colors
  colors <- append(NA,b_clrs)
  tableplot$neglogofP2 <- append(NA,dat$neglogofP)
  
  
  #Point : (22 can be used to fill outline color)
  point<-c(NA)
  if (nrow(dat)>1){
    for (i in 2:(nrow(dat))){point<-append(point,15)}
    point <- append(point,18)
  }else{ point<-append(point,15)}
  
  #Box size
  if (nrow(dat)>1) {
    total <- as.numeric(dat$Total)
    boxsize=c(0,2.5*log10(total)) # 0 for the header line
  } else {boxsize=2.5} # boxsize = 2 default
  
  
  #Theme
  theme_set(theme_bw())
  
  theme_update(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black"), #border around the plot
    panel.background = element_rect(fill="#E7E7E7"),  #edebe9 #DEDAD7 #F6F5F4 #f7f3ef #BDB6B0 ##ECECEC #E7E7E7 #C6C6C5
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size=16),
    plot.margin = unit(c(0, 0, 0, 0), "lines"),
    axis.text=element_text(size=10)
  )
  
  ##G2 Plot 2 - forestplot
  g2 <- ggplot(tableplot,aes(Fold.Change,group)) + 
    geom_errorbarh(aes(xmax = Upper_bound_CI, xmin = Lower_bound_CI), height = 0.40, 
                   colour= ifelse(tableplot$group == tableplot$group[nrow(tableplot)] & tableplot$neglogofP2 == tableplot$neglogofP2[nrow(tableplot)],"#E7E7E7",colors),
                   na.rm=TRUE, size=1.2) + #0.15 #removed error bar for the final integration
    geom_point(aes(fill=neglogofP2),size=boxsize, shape=point, colour = colors, na.rm=TRUE)+
    #geom_point(aes(fill=neglogofP2),size=boxsize, shape=point, colour= ifelse(tableplot$group == tableplot$group[nrow(tableplot)] & tableplot$neglogofP2 == tableplot$neglogofP2[nrow(tableplot)],colors,"black"), na.rm=TRUE)+
    geom_vline(xintercept = 1, linetype = "longdash") + scale_x_continuous(breaks = xticks) + labs(x= "Fold Change",y="",fill="-log(Q-Value)") +
    scale_fill_gradientn(colours=reds,limits=c(0,8),breaks = seq(0,8),oob=squish) + #squish from scales package manages values out of range
    guides(fill = guide_colourbar(barheight = 10)) 
  
  #Legend size
  if (nrow(dat) > 5) {g2 <- g2 + theme(legend.text = element_text(size=10),legend.title = element_text(size=14))}
  else {g2 <- g2 + theme(legend.text = element_text(size=8),legend.title = element_text(size=12))}
  
  #Add title if specified  
  title <- ggdraw() + draw_label(paste0(ptitle, curr_gene),fontface = 'bold') #Initialize title
  
  #Bind two plots together and then bind this with the title
  #fplot <- plot_grid(g1,g2,ncol=2,rel_widths = c(2.1,1.0),align = "h") #rel_widths #align="h",axis="tblr",
  fplot <- plot_grid(g1,g2,ncol=2,rel_widths = c(2.1,1.3),align = "h") #increased width of forestplot
  plot_grid(title,fplot,ncol = 1, rel_heights = c(0.1, 1),align = "v")
  
}

#1 PX = 0.0104166653543 in.
getHeight_func <- function(dat){
  height=25 #27 optimum
  height<- 150 + height*(abs(nrow(dat)-2)) #Title: 3x + pcomb object (5px) #110
  return(height)
}