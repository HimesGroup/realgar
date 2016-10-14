library(shiny)

# Define UI for application that displays GEO results based on user choices
shinyUI(fluidPage(
    
    titlePanel(h2("Asthma Gene Expression Explorer", align="center")),    

    #Top row to gather input
    fluidRow(
        column(12, class="well",
               p("Differential expression results for asthma vs. non-asthma samples from individual publicly available microarray datasets in GEO. Select input categories, then click the update button."),
               
               column(4,checkboxGroupInput(inputId="Tissue", label="Tissue:", choices=c("Broncial Epithelium"="BE",
                                                                                       "Nasal Epithelium"="NE","CD4"="CD4","CD8"="CD8","PBMC"="PBMC","White Blood Cell"="WBC", "Airway Smooth Muscle"="ASM",
                                                                                       "BAL"="BAL", "Lung"="Lung"),selected="BE" )),
               column(5,checkboxGroupInput(inputId="Asthma", label="Asthma Type:", choices=c("Allergic Asthma"="allergic_asthma",
                                                                                                                      "Severe Asthma"="severe_asthma","Asthma"="asthma","Mild Asthma"="mild_asthma","Non-Allergic Asthma"=
                                                                                                                          "non_allergic_asthma", "Asthma with Rhinitis"="asthma_and_rhinitis"), selected="asthma")),
               column(3,radioButtons(inputId="allP.Values",label="Probes to Include:", choices=c("All available"="TRUE",
                                                                                                "Probe with minimum p-value"="FALSE"),selected="FALSE"),
                      textInput(inputId="curr_gene",label="Gene (official gene symbol):", value= "GAPDH"),
                      #insert an action Button
                      actionButton(inputId="updategraph",label="Update"),
                      hr(),
                      fixedRow(img(src="http://shiny.rstudio.com/tutorial/lesson2/www/bigorb.png", height=42, width=48),
                               "Created with ",
                               a("RStudio's shiny", href="http://www.rstudio.com/shiny")))),
        hr(),
        tabsetPanel(
            tabPanel("Study Information", 
                     fluidRow(align="center",
                              br(),
                              fixedRow(align="center", downloadButton(outputId="dGT", label="Download Table")),
                              br(),
                              tableOutput(outputId="user.choice"))),
            tabPanel("Differential Expression Results",
                     fixedRow(align="center",
                              br(),
                              column(6, downloadButton(outputId="dplotFC",label="Download Fold-Change Plot")),
                              column(6, downloadButton(outputId="dplotNLOP", label="Download P-Value Plot"))),
                     hr(),
                     fluidRow(align="center",
                              column(6, plotOutput(outputId="heatmapMATFC",width="100%", height="800px")),
                              column(6, plotOutput(outputId="heatmapMATNLOP",width="100%", height="800px")))),
            tabPanel("Selected Gene Results",
                     br(),
                     fluidRow(align="center",
                              tableOutput(outputId="tableforgraph"), align="center")))))
)

