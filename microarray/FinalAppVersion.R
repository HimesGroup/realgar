#User Interface

#Load Shiny 
#install.packages("shiny")
library(shiny)
library(dplyr)
library(lattice)
library(stringr)
library(RColorBrewer)

#Colors for the heatmap
cols2 <-colorRampPalette(c("white","mediumspringgreen"))


#Load data files
dataset.info <- read.csv("../databases/microarray_data_infosheet.csv")
dataset.info$Unique_ID <- apply(dataset.info[, c("GEO_ID", "Tissue", "Asthma")], 1, paste, collapse="_")
for (i in dataset.info$Unique_ID){
    assign(i, read.csv(paste0("../databases/microarray_results/", i,".csv"), sep=""))}
output.table <- data.frame()

#Template
ui=fluidPage(

titlePanel(h2("Asthma Gene Explorer", align="center")),    
br(),
p("This app is used to discover the significance of genes based on different studies from GEO. To use this app please
  choose a gene, tissue(s), and asthma type(s). Once you have selected your input and typed the official gene symbol, click the update button.
  This button will update the differential expression results as well as the selected gene data accordign to the Microarray Inclusion Criteria."),
hr(),

#User Input
fluidRow(
        column(12,class="well",h4("Microarray Inclusion Criteria", align="center"),
            column(4,checkboxGroupInput(inputId="Tissue", label="Tissue", choices=c("Broncial Epithelium"="BE",
                "Nasal Epithelium"="NE","CD4"="CD4","CD8"="CD8","PBMC"="PBMC","White Blood Cell"="WBC", "Airway Smooth Muscle"="ASM",
                "BAL"="BAL", "Lung"="Lung"),selected="BE" )),
            column(5,checkboxGroupInput(inputId="Asthma", label="Asthma Types according to GEO deposit", choices=c("Allergic Asthma"="allergic_asthma",
                 "Severe Asthma"="severe_asthma","Asthma"="asthma","Mild Asthma"="mild_asthma","Non Allergic Asthma"=
                 "non_allergic_asthma", "Asthma and Rhinitis"="asthma_and_rhinitis"),selected="asthma"),
                 h6("*Note that for this application the data was always analyzed with a reference patten of asthma vs [other type of asthma].")),
            column(3,radioButtons(inputId="allP.Values",label="Probes to Include", choices=c("Use all probes for each gene type"="TRUE",
                                                                                             "Probe with Minimum P.Value"="FALSE"),selected="FALSE"),
                textInput(inputId="curr_gene",label="Type the Official Gene Symbol", value= "GAPDH"),
                #insert an action Button
                actionButton(inputId="updategraph",label="Update"),
                hr(),
                fixedRow(img(src="http://shiny.rstudio.com/tutorial/lesson2/www/bigorb.png", height=42, width=48),
                "Created with ",
                a("RStudio's shiny", href="http://www.rstudio.com/shiny")))),
        hr(),
        tabsetPanel(
            tabPanel("Corresponding Number of Selected Studies from GEO ", fluidRow(h4("Description of Studies Selected from GEO that meet Criteria", align="center"),
                                                             br(),
                                                             fixedRow(align="center",downloadButton(outputId="dGT",label="Download Table")),
                                                             br(),
                                                             tableOutput(outputId="user.choice"))),
            tabPanel("Differential Expression Results",
                     fixedRow(align="center",
                         column(6, downloadButton(outputId="dplotFC",label="Download Fold Change Plot")),
                         column(6, downloadButton(outputId="dplotNLOP", label="Download P.Value Plot"))),
                     hr(),
                     
                     fluidRow(align="center",
                         column(6,plotOutput(outputId="heatmapMATFC",width="100%", height="800px")),
                         column(6,plotOutput(outputId="heatmapMATNLOP",width="100%", height="800px")))),
            tabPanel("Selected Gene Results",fluidRow(h4("Corresponding Differential Expression Results", align="center"),
                                                   tableOutput(outputId="tableforgraph"),align="center")))))
        

server=function(input,output){
    #Create the Table to Show the choices for the user
    Userdataset.info= reactive({
        dataset.info=subset(dataset.info,dataset.info$Tissue %in% input$Tissue & dataset.info$Asthma %in% input$Asthma)
        dataset.info
    })
    Userdataset.info2= reactive({rename(Userdataset.info(),'Study ID'=Unique_ID, 'Sample Size'=Total, 'Asthma Type'=Asthma)})
    output$user.choice=renderTable(Userdataset.info2())
    curr_gene=reactive({input$curr_gene})
    allP.Values= reactive({input$allP.Values})
    
    #Use the observe function to generate the table and graphs
    output.tableforplot=eventReactive(input$updategraph,{
        myfunction <- function(x){
            if(allP.Values()=="FALSE"){
                x %>%
                    filter(SYMBOL==curr_gene()) %>%
                    select(logFC, P.Value, adj.P.Val) %>%
                    filter(P.Value==min(P.Value))
            } else {
                x %>%
                    filter(SYMBOL==curr_gene()) %>%
                    select(logFC, P.Value, adj.P.Val)
            }
        }
        #In for loop create if then statement 
        for (i in Userdataset.info()$Unique_ID){
            curr.gene.data=get(i,environment())
            if(any(tbl_vars(curr.gene.data)=="qValuesBatch")){
                curr.gene.data=( curr.gene.data%>%
                                     select(-P.Value,-adj.P.Val)%>%
                                     rename(P.Value=pValuesBatch)%>%
                                     rename(adj.P.Val=qValuesBatch)
                )}
            curr.gene.data <- myfunction(curr.gene.data)
            if(nrow(curr.gene.data) > 0){
                curr.gene.data <- cbind(Unique_ID=i, curr.gene.data)
                output.table <- rbind(output.table, curr.gene.data)
            }}
        
    #preparing the data for the Image
    #calculate the Fold Change
    output.table=mutate(output.table, Fold_Change=2^(logFC) )
    output.table=mutate(output.table, neglogofP=(-log10(P.Value)) )
    #Sort the Data
    output.table=output.table[order(output.table$Fold_Change),]})
    output.tableforplot2=reactive({output.tableforplot()%>%
                                       rename(' '=Fold_Change, ' '=neglogofP)})
    
    ##
    
    data=reactive({select(output.tableforplot(), Unique_ID, P.Value,Fold_Change,neglogofP )
                                      output.tableforplot()[rev(rownames(output.tableforplot())),]})
    data2=reactive({data()%>%
                        select(Unique_ID,Fold_Change,adj.P.Val,P.Value )%>%
                        mutate(Fold_Change=round(Fold_Change,digits=2),adj.P.Val=format(adj.P.Val, scientific=TRUE, digits=3), P.Value =format(P.Value, scientific=TRUE, digits=3))%>%
                        rename('Study ID'=Unique_ID, 'P Value'=P.Value, 'Q Value'=adj.P.Val, 'Log 2 Fold Change'=Fold_Change)})
    
    
    output$tableforgraph=renderTable(data2())
    heatmapMAT=eventReactive(input$updategraph,{output.tableforplot2()})
    heatmapMATFC=reactive({data.matrix(t(heatmapMAT()[5]))})
    heatmapMATNLOP=reactive({data.matrix(t(heatmapMAT()[6]))})
    
    minFC=reactive({if(max(heatmapMATFC())<=1.5 & min(heatmapMATFC())>=-1.5){
        minFC=-1.5
    } else {
        minFC=min(heatmapMATFC())}})
    maxFC=reactive({if(max(heatmapMATFC())<=1.5 & min(heatmapMATFC())>=-1.5){
        maxFC=1.5
    } else {
        maxFC=max(heatmapMATFC())}})
    
    minNLOP=reactive({0})
    maxNLOP=reactive({if(max(heatmapMATNLOP())<=1.5 & min(heatmapMATNLOP())>=-1.5){
        maxNLOP=1.5
    } else {
        maxNLOP=max(heatmapMATNLOP())}})
      
    #ProducePlots
    graphgene=eventReactive(input$updategraph,{curr_gene()})
    
    imageFC= reactive({print(levelplot(heatmapMATFC(),
                        col.regions=cm.colors,
                        xlab =NULL,
                        ylab="Corresponging Number for Plot Data",
                        main = paste0("Fold Change in Gene ", graphgene()),
                        scales=list(x=list(cex=1,rot=35,tck = c(0,0,0,0)), 
                                    y=list(tck = c(1,0,0,0))),
                        at=seq(minFC(), maxFC(), length.out=100)))})
    
    output$heatmapMATFC= renderPlot(imageFC())
    imageNLOP=reactive({ print(levelplot(heatmapMATNLOP(),
                        col.regions=cols2,
                        xlab =NULL,
                        ylab="Corresponging Number for Plot Data",
                        main = paste0("Negative Log P-Value for Gene ", graphgene()),
                        scales=list(x=list(cex=1,rot=35,tck = c(0,0,0,0)), 
                                    y=list(tck = c(1,0,0,0))),
                        at=seq(minNLOP(), maxNLOP(),length.out=100)))})
    output$heatmapMATNLOP=renderPlot(imageNLOP())
    
    #Downloading things
   output$dplotFC=downloadHandler(
       filename= function(){paste("Fold Change Graph for", "_", graphgene(), "_", Sys.Date(), ".png", sep='')},
       content=function(file){
           png(file)
           print(imageFC())
           dev.off() 
       })
   output$dplotNLOP=downloadHandler(
       filename= function(){paste("Negative Log10 Graph for", "_", graphgene(), "_", Sys.Date(), ".png", sep='')},
       content=function(file){
           png(file)
           print(imageNLOP())
           dev.off() 
       })
   
   output$dGT=downloadHandler(
       filename = function() {
           paste('Graph Data for',graphgene(), Sys.Date(), '.csv', sep="")
       },
       content = function(file) {
           write.csv(data(), file)
       })
   
   
}
shinyApp(ui=ui, server=server)







