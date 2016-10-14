library(shiny)
library(dplyr)
library(lattice)
library(stringr)
library(RColorBrewer)

#Load data files
dataset.info <- read.csv("../databases/microarray_data_infosheet.csv")
dataset.info$Unique_ID <- apply(dataset.info[, c("GEO_ID", "Tissue", "Asthma")], 1, paste, collapse="_")
for (i in dataset.info$Unique_ID){
    assign(i, read.csv(paste0("../databases/microarray_results/", i,".csv"), sep=""))}
output.table <- data.frame()

#Colors for the heatmap
cols2 <-colorRampPalette(c("white","mediumspringgreen"))

# Define server logic required to display results of multiple GEO datasets
shinyServer(function(input,output) {
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
})
