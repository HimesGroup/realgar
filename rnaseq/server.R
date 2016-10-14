library(shiny)
library(ggplot2)
library(data.table)
library(magrittr)
library(dplyr)

#Enable Browser Mode on Errors for Debugging Purposes
#options(error = browser)

#Load data files - gene names and dataset info
sras <- fread("../databases/dataset_info.csv") %>% tbl_df
gene_names <- fread("../databases/gene_names.csv") %>% tbl_df

#sleuth output: fold changes, pvalues, conditions - for data table beneath plot
de <- list()
de[["SRP005411"]] <- fread("../databases/SRP005411_full_sleuth.txt", sep = " ") %>% tbl_df
de[["SRP043162"]] <- fread("../databases/SRP043162_full_sleuth.txt", sep = " ") %>% tbl_df
de[["SRP033351"]] <- fread("../databases/SRP033351_full_sleuth.txt", sep = " ") %>% tbl_df


#kallisto results - by-transcript TPMs, used for the plots
tpms <- list()
tpms[["SRP005411"]] <- fread("../databases/SRP005411_full_kallisto.txt", sep = " ") %>% tbl_df
tpms[["SRP043162"]] <- fread("../databases/SRP043162_full_kallisto.txt", sep = " ") %>% tbl_df
tpms[["SRP033351"]] <- fread("../databases/SRP033351_full_kallisto.txt", sep = " ") %>% tbl_df

print("Loaded datasets.")


server <- shinyServer(function(input, output, session) {
    
    updateSelectizeInput(session, "gene", choices=gene_names$name, selected="GAPDH", server=TRUE)
    
    output$gene <- renderPrint({
        input$gene
    })
    
    output$studyText <- renderUI({
        if (!is.null(input$debugcode) && (input$debugcode == "studyText")) {
            browser()
        }
        
        sras %>% filter(Tissue == input$tissue) %$% 
            p("Data used is available in the SRA under accession ",
              a(paste0(SRA_ID, ","), href=paste0("http://www.ncbi.nlm.nih.gov/sra/?term=", SRA_ID)),
              "and corresponds to ",
              Description,
              if (PMID != "-") {
                  HTML(paste0("More details were published <a href=http://www.ncbi.nlm.nih.gov/pubmed/?term=", PMID, ">here.</a>"))
              }
            )
    })
    
    getGeneBoxPlot <-reactive({
        
        x <- sras %>% 
            filter(Tissue == input$tissue) %$% 
            SRA_ID
        
        curr_data <- tpms[[x]] %>%
            filter(ext_gene == input$gene) 
        
        if (nrow(curr_data) > 0) {
            
            curr_data$average_tpm <- vector(length=nrow(curr_data))
            
            for (i in 1:nrow(curr_data)) {
                curr_data$average_tpm[i] <- mean(curr_data[which(curr_data$target_id == curr_data$target_id[i]),]$tpm)
            }
            
            curr_data <- curr_data %>%
                filter(average_tpm > 1)
            
            gene_plot <- ggplot(curr_data, aes(x = condition, y = tpm, fill=condition)) + 
                geom_boxplot(outlier.colour=NA, lwd=0.2, color="grey18") + 
                stat_boxplot(geom ='errorbar', color="grey18") + 
                geom_jitter(size=0.8, width=0.2) +
                facet_wrap(~target_id) + 
                guides(fill=FALSE) + 
                theme_bw() +  
                labs(title=input$gene) + 
                labs(x="condition") + labs(y="TPM") + 
                theme(text = element_text(size=9), 
                      strip.text.x = element_text(size = 10), 
                      axis.text.x = element_text(angle = 90, hjust = 1, size=12),
                      axis.text.y = element_text(size=9),
                      title = element_text(size=12),
                      axis.title.x = element_text(size=12),
                      axis.title.y = element_text(size=12))
            if (nrow(curr_data) > 0) {gene_plot}
        }
    })    
    
    output$GeneBoxPlot <- renderPlot(width = 650, height = 500, {
        if (!is.null(input$debugcode) && (input$debugcode == "geneBoxPlot")) {
            browser()
        }
        gbp <- getGeneBoxPlot()
        validate(need(gbp, message="Plot not available. Either the gene is not in reference database, less than 47% of samples from this dataset had 5 or more reads in the gene's transcript(s), or all of the gene transcripts had low expression (average TPM<1)."))
        gbp
    })
    
    output$diffResults <- renderDataTable({
        if (!is.null(input$debugcode) && (input$debugcode == "diffResults")) {
            browser()
        }
        
        sras %>% filter(Tissue == input$tissue) %$% de[[SRA_ID]] %>% subset(ext_gene %in% input$gene) %>%  
            mutate(b=round(b, digits=2), qval=format(qval, scientific=TRUE, digits=3)) %>% 
            arrange(-b, qval) %>% 
            select(ext_gene, target_id, Comparison, b, qval) %>% #rearrange columns in desired order
            rename(`Gene`=ext_gene, `Transcript`=target_id, `Beta value`=b, `Q Value`=qval) 
        
        
    }, options=list(paging=FALSE, searching=FALSE)
    )
    
    output$downloadPic <- downloadHandler(
        filename = function() {paste(input$tissue, "_", input$gene, "_", Sys.Date(), '.png', sep='')},
        content = function(file) {
            png(file, width=10, height=6, units="in", res=600)
            print(getGeneBoxPlot())
            dev.off()
        },
        contentType = 'image/png'
    )
    
})
