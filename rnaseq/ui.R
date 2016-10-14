library(shiny)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  br(),
  
  # Application title
  #titlePanel(h2("Lung Cell Transcriptome Explorer", align="center")),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      p("Expression results for individual gene transcripts using publicly available RNA-Seq datasets."),
      br(),
      
      selectInput("tissue", label = "Tissue:",
                  choices = list("Airway Smooth Muscle 1" = "ASM1", "Airway Smooth Muscle 2" = "ASM2", "Airway Epithelium" = "AE"), 
                  selected="AE"),
      
      selectInput("gene", "Gene (official gene symbol):", multiple=FALSE, choices=NULL),
      
      br(),
      
      p("Abundance of each transcript is expressed in Transcripts per Million Reads Mapped (TPM) as computed by",
        a("kallisto.", href="https://pachterlab.github.io/kallisto/"),
        "Differential expression results were obtained with",
        a("sleuth.", href="https://pachterlab.github.io/sleuth/"),
        "The table includes transcripts for which at least 47% of samples had 5 or more reads.",
        "The plot includes transcripts in the table with average TPM > 1 across all conditions."),
      
      br(),
      #textInput("debugcode", "Debug:", ""),     
      br(),
      br(),
      br(),
      img(src="bigorb.png", height=42, width=48),
      "Created with ",
      a("RStudio's shiny", href="http://www.rstudio.com/shiny")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(p(""),
              #h5("Results For Each Experimental Condition", align="center"),
              #h3(textOutput("gene", container=span)),
              plotOutput("GeneBoxPlot"),
              br(),br(),br(),br(),br(),
              uiOutput("studyText"),
              downloadButton('downloadPic', 'Download Figure'),
              br(),
              br(),
              h5("Corresponding Differential Expression Results", align="center"),
              dataTableOutput("diffResults"),
              br()
    )
  )
))
